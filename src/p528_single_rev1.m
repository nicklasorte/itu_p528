function [res, rtn] = p528_single_rev1(d_km, h_1_m, h_2_m, f_mhz, T_pol, p)
%P528_SINGLE_REV1 One ITU-R P.528-5 basic transmission loss calculation (scalar inputs).
%   [res, rtn] = P528_SINGLE_REV1(d_km, h_1_m, h_2_m, f_mhz, T_pol, p)
%   follows Annex 2, sec. 3 step by step and returns the result struct
%   (d_km, A_db, A_fs_db, A_a_db, theta_h1_rad, propagation_mode, warnings)
%   and the reference implementation's return code. Use P528_REV1 for
%   array inputs and error handling.
%
%   Port of src/p528/P528.cpp and ValidateInputs.cpp (github.com/NTIA/p528,
%   commit 74870ab), modified November 2026 for MATLAB: see README.
c = p528_constants_rev1();
res = struct('d_km', 0, 'A_db', 0, 'A_fs_db', 0, 'A_a_db', 0, 'theta_h1_rad', 0, ...
    'propagation_mode', c.MODE_NOT_SET, 'warnings', c.WARNING_NO_WARNINGS);

[err, res.warnings] = validate_inputs(d_km, h_1_m, h_2_m, f_mhz, T_pol, p, res.warnings, c);
if err ~= c.SUCCESS
    if err == c.ERROR_HEIGHT_AND_DISTANCE
        rtn = c.SUCCESS;            % reference returns zero loss with SUCCESS
    else
        rtn = err;
    end
    return
end

% Step 1: terminal geometry
t1 = p528_terminal_geometry_rev1(f_mhz, h_1_m / 1000);
t2 = p528_terminal_geometry_rev1(f_mhz, h_2_m / 1000);

% Step 2
path.d_ML_km = t1.d_r_km + t2.d_r_km;                                   % Eqn 3-1
path.d_0_km = 0;
path.d_d_km = 0;

% Step 3: smooth-Earth diffraction line
d_3_km = path.d_ML_km + 0.5 * (c.a_e_km^2 / f_mhz)^c.third;           % Eqn 3-2
d_4_km = path.d_ML_km + 1.5 * (c.a_e_km^2 / f_mhz)^c.third;           % Eqn 3-3
A_3_db = p528_smooth_earth_diffraction_rev1(t1.d_r_km, t2.d_r_km, f_mhz, d_3_km, T_pol);
A_4_db = p528_smooth_earth_diffraction_rev1(t1.d_r_km, t2.d_r_km, f_mhz, d_4_km, T_pol);
M_d = (A_4_db - A_3_db) / (d_4_km - d_3_km);                           % Eqn 3-4
A_d0 = A_4_db - M_d * d_4_km;                                          % Eqn 3-5
A_dML_db = (M_d * path.d_ML_km) + A_d0;                                % Eqn 3-6
path.d_d_km = -(A_d0 / M_d);                                           % Eqn 3-7

% Step 4: line-of-sight region
if path.d_ML_km - d_km > 0.001
    res.propagation_mode = c.MODE_LOS;
    [~, ~, res] = p528_line_of_sight_rev1(path, t1, t2, f_mhz, -A_dML_db, p, d_km, T_pol, res);
    rtn = return_code(res.warnings, c);
    return
end

% K_LOS from just inside the horizon
[path, ~, res, K_LOS] = p528_line_of_sight_rev1(path, t1, t2, f_mhz, -A_dML_db, p, path.d_ML_km - 1, T_pol, res);

% Step 6: transhorizon search
[M_d, A_d0, d_crx_km, CASE, res.warnings] = p528_transhorizon_search_rev1(path, t1, t2, f_mhz, A_dML_db, M_d, A_d0, res.warnings);

% Step 7: terrain attenuation
A_d_db = M_d * d_km + A_d0;                                            % Eqn 3-14
tropo = p528_troposcatter_rev1(t1, t2, d_km, f_mhz);
if d_km < d_crx_km
    A_T_db = A_d_db;
    res.propagation_mode = c.MODE_DIFFRACTION;
elseif CASE == c.CASE_1
    if tropo.A_s_db <= A_d_db
        A_T_db = tropo.A_s_db;
        res.propagation_mode = c.MODE_SCATTERING;
    else
        A_T_db = A_d_db;
        res.propagation_mode = c.MODE_DIFFRACTION;
    end
else
    A_T_db = tropo.A_s_db;
    res.propagation_mode = c.MODE_SCATTERING;
end

% Variability (f_theta_h is unity for transhorizon paths)
f_theta_h = 1;
v = p528_variability_rev1("long_term", t1.d_r_km, t2.d_r_km, d_km, f_mhz, p, f_theta_h, -A_T_db);
Y_e_db = v(1);
v = p528_variability_rev1("long_term", t1.d_r_km, t2.d_r_km, d_km, f_mhz, 50, f_theta_h, -A_T_db);
Y_e_50_db = v(1);
ANGLE = 0.02617993878;                                                 % 1.5 deg
if tropo.theta_s >= ANGLE
    K_t_db = 20;
elseif tropo.theta_s <= 0.0
    K_t_db = K_LOS;
else
    K_t_db = (tropo.theta_s * (20.0 - K_LOS) / ANGLE) + K_LOS;
end
Y_pi_db = p528_variability_rev1("nakagami_rice", K_t_db, p);
Y_total_db = p528_variability_rev1("combine", Y_e_50_db, Y_e_db, 0.0, Y_pi_db, p);

% Atmospheric absorption and free-space loss for the transhorizon path
rv = p676_slant_path_rev1(f_mhz / 1000, 0, tropo.h_v_km, pi/2);
res.A_a_db = t1.A_a_db + t2.A_a_db + 2 * rv.A_gas_db;                 % Eqn 3-17
r_fs_km = t1.a_km + t2.a_km + 2 * rv.a_km;                            % Eqn 3-18
res.A_fs_db = 20.0 * log10(f_mhz) + 20.0 * log10(r_fs_km) + 32.45;    % Eqn 3-19

res.d_km = d_km;
res.A_db = res.A_fs_db + res.A_a_db + A_T_db - Y_total_db;            % Eqn 3-20
res.theta_h1_rad = -t1.theta_rad;
rtn = return_code(res.warnings, c);
end

function rtn = return_code(warnings, c)
if warnings == c.WARNING_NO_WARNINGS
    rtn = c.SUCCESS;
else
    rtn = c.SUCCESS_WITH_WARNINGS;
end
end

function [err, warnings] = validate_inputs(d_km, h_1_m, h_2_m, f_mhz, T_pol, p, warnings, c)
err = c.SUCCESS;
if d_km < 0, err = c.ERROR_VALIDATION_D_KM; return; end
if h_1_m < 1.5 || h_1_m > 80000, err = c.ERROR_VALIDATION_H_1; return; end
if h_2_m < 1.5 || h_2_m > 80000, err = c.ERROR_VALIDATION_H_2; return; end
if h_1_m > 20000 && h_1_m <= 80000, warnings = bitor(warnings, c.WARNING_HEIGHT_LIMIT_H_1); end
if h_2_m > 20000 && h_2_m <= 80000, warnings = bitor(warnings, c.WARNING_HEIGHT_LIMIT_H_2); end
if h_1_m > h_2_m, err = c.ERROR_VALIDATION_TERM_GEO; return; end
if f_mhz < 100, err = c.ERROR_VALIDATION_F_MHZ_LOW; return; end
if f_mhz > 30000, err = c.ERROR_VALIDATION_F_MHZ_HIGH; return; end
if T_pol ~= c.POL_HORIZONTAL && T_pol ~= c.POL_VERTICAL, err = c.ERROR_VALIDATION_POLARIZATION; return; end
if p < 1, err = c.ERROR_VALIDATION_PERCENT_LOW; return; end
if p > 99, err = c.ERROR_VALIDATION_PERCENT_HIGH; return; end
if h_1_m == h_2_m && d_km == 0, err = c.ERROR_HEIGHT_AND_DISTANCE; return; end
end
