function [path, prm, res, K_LOS] = p528_line_of_sight_rev1(path, t1, t2, f_mhz, A_dML_db, p, d_km, T_pol, res)
%P528_LINE_OF_SIGHT_REV1 Line-of-sight loss with variability, ITU-R P.528-5 Annex 2 secs. 5-8 and 13.
%   [path, prm, res, K_LOS] = P528_LINE_OF_SIGHT_REV1(path, t1, t2, f_mhz,
%   A_dML_db, p, d_km, T_pol, res) fills res.A_db, res.A_fs_db, res.A_a_db,
%   res.d_km, and res.theta_h1_rad, sets path.d_0_km, and returns the
%   Nakagami-Rice K for the line-of-sight region.
%
%   Port of src/p528/LineOfSight.cpp (github.com/NTIA/p528).
lambda_km = 0.2997925 / f_mhz;                          % Eqn 6-1
terminate = lambda_km / 1e6;

psi_limit = find_psi_at_delta_r(lambda_km / 2, t1, t2, terminate);
d_y6_km = find_distance_at_delta_r(lambda_km / 6, t1, t2, terminate);

% d_0 heuristic [Eqns 8-2, 8-3]
if t1.d_r_km >= path.d_d_km || path.d_d_km >= path.d_ML_km
    if t1.d_r_km > d_y6_km || d_y6_km > path.d_ML_km
        path.d_0_km = t1.d_r_km;
    else
        path.d_0_km = d_y6_km;
    end
elseif path.d_d_km < d_y6_km && d_y6_km < path.d_ML_km
    path.d_0_km = d_y6_km;
else
    path.d_0_km = path.d_d_km;
end

% Tune d_0 forward 1 m at a time without leaving the LOS region
d_temp_km = path.d_0_km;
while true
    psi = find_psi_at_distance(d_temp_km, t1, t2);
    los = p528_ray_optics_rev1(t1, t2, psi);
    if los.d_km >= path.d_0_km || (d_temp_km + 0.001) >= path.d_ML_km
        path.d_0_km = los.d_km;
        break
    end
    d_temp_km = d_temp_km + 0.001;
end

% Loss at d_0
psi_d0 = find_psi_at_distance(path.d_0_km, t1, t2);
prm = p528_ray_optics_rev1(t1, t2, psi_d0);
prm = p528_get_path_loss_rev1(psi_d0, path, f_mhz, psi_limit, A_dML_db, 0, T_pol, prm);
A_LOS_d0 = prm.A_LOS_db;

% Loss at the requested distance
psi = find_psi_at_distance(d_km, t1, t2);
prm = p528_ray_optics_rev1(t1, t2, psi);
[prm, R_Tg] = p528_get_path_loss_rev1(psi, path, f_mhz, psi_limit, A_dML_db, A_LOS_d0, T_pol, prm);

% Atmospheric absorption
slant = p676_slant_path_rev1(f_mhz / 1000, t1.h_r_km, t2.h_r_km, pi/2 - prm.theta_h1_rad);
res.A_a_db = slant.A_gas_db;

% Free-space loss [Eqn 6-4]
res.A_fs_db = 20.0 * log10(prm.r_0_km) + 20.0 * log10(f_mhz) + 32.45;

% Variability [Eqn 13-1]
if prm.theta_h1_rad <= 0.0
    f_theta_h = 1.0;
elseif prm.theta_h1_rad >= 1.0
    f_theta_h = 0.0;
else
    f_theta_h = max(0.5 - (1 / pi) * (atan(20.0 * log10(32.0 * prm.theta_h1_rad))), 0);
end
v = p528_variability_rev1("long_term", t1.d_r_km, t2.d_r_km, d_km, f_mhz, p, f_theta_h, prm.A_LOS_db);
Y_e_db = v(1);
v = p528_variability_rev1("long_term", t1.d_r_km, t2.d_r_km, d_km, f_mhz, 50, f_theta_h, prm.A_LOS_db);
Y_e_50_db = v(1); A_Y = v(2);

if A_Y <= 0.0                                          % Eqn 13-2
    F_AY = 1.0;
elseif A_Y >= 9.0
    F_AY = 0.1;
else
    F_AY = (1.1 + (0.9 * cos((A_Y / 9.0) * pi))) / 2.0;
end
if prm.delta_r_km >= (lambda_km / 2.0)                 % Eqn 13-3
    F_delta_r = 1.0;
elseif prm.delta_r_km <= lambda_km / 6.0
    F_delta_r = 0.1;
else
    F_delta_r = 0.5 * (1.1 - (0.9 * cos(((3.0 * pi) / lambda_km) * (prm.delta_r_km - (lambda_km / 6.0)))));
end
R_s = R_Tg * F_delta_r * F_AY;                         % Eqn 13-4
Y_pi_99_db = 10.0 * log10(f_mhz * slant.a_km^3) - 84.26;   % Eqn 13-5
K_t = p528_variability_rev1("find_k_99", Y_pi_99_db);
W_a = 10.0^(K_t / 10.0);                               % Eqn 13-6
W_R = R_s^2 + 0.01^2;                                  % Eqn 13-7
W = W_R + W_a;                                         % Eqn 13-8
if W <= 0.0                                            % Eqn 13-9
    K_LOS = -40.0;
else
    K_LOS = max(10.0 * log10(W), -40.0);
end
Y_pi_db = p528_variability_rev1("nakagami_rice", K_LOS, p);
Y_total_db = -p528_variability_rev1("combine", Y_e_50_db, Y_e_db, 0.0, Y_pi_db, p);

res.d_km = prm.d_km;
res.A_db = res.A_fs_db + res.A_a_db - prm.A_LOS_db + Y_total_db;
res.theta_h1_rad = prm.theta_h1_rad;
end

function psi = find_psi_at_distance(d_km, t1, t2)
if d_km == 0
    psi = pi/2;
    return
end
psi = pi/2;
delta_psi = -pi/4;
while true
    psi = psi + delta_psi;
    tmp = p528_ray_optics_rev1(t1, t2, psi);
    d_psi_km = tmp.d_km;
    if d_psi_km > d_km
        delta_psi = abs(delta_psi) / 2;
    else
        delta_psi = -abs(delta_psi) / 2;
    end
    if ~(abs(d_km - d_psi_km) > 1e-3 && abs(delta_psi) > 1e-12)
        break
    end
end
end

function psi = find_psi_at_delta_r(delta_r_km, t1, t2, terminate)
psi = pi/2;
delta_psi = -pi/4;
for iter = 1:2000
    psi = psi + delta_psi;
    tmp = p528_ray_optics_rev1(t1, t2, psi);
    if tmp.delta_r_km > delta_r_km
        delta_psi = -abs(delta_psi) / 2;
    else
        delta_psi = abs(delta_psi) / 2;
    end
    if ~(abs(tmp.delta_r_km - delta_r_km) > terminate)
        break
    end
end
check_converged(tmp.delta_r_km, delta_r_km, terminate);   % the C++ loop is unbounded; fail loud instead of hanging
end

function d = find_distance_at_delta_r(delta_r_km, t1, t2, terminate)
psi = pi/2;
delta_psi = -pi/4;
for iter = 1:2000
    psi = psi + delta_psi;
    tmp = p528_ray_optics_rev1(t1, t2, psi);
    if tmp.delta_r_km > delta_r_km
        delta_psi = -abs(delta_psi) / 2;
    else
        delta_psi = abs(delta_psi) / 2;
    end
    if ~(abs(tmp.delta_r_km - delta_r_km) > terminate)
        break
    end
end
check_converged(tmp.delta_r_km, delta_r_km, terminate);
d = tmp.d_km;
end

function check_converged(x, target, tol)
if abs(x - target) > tol
    error('p528:noConvergence', 'delta_r search did not converge in 2000 steps (target %.3g km).', target);
end
end
