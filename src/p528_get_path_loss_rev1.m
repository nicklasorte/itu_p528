function [prm, R_Tg] = p528_get_path_loss_rev1(psi_rad, path, f_mhz, psi_limit, A_dML_db, A_d_0_db, T_pol, prm)
%P528_GET_PATH_LOSS_REV1 Line-of-sight loss for grazing angle psi, ITU-R P.528-5 Annex 2 sec. 8.
%   [prm, R_Tg] = P528_GET_PATH_LOSS_REV1(psi_rad, path, f_mhz, psi_limit,
%   A_dML_db, A_d_0_db, T_pol, prm) sets prm.A_LOS_db and returns the total
%   ground reflection coefficient R_Tg.
%
%   Port of src/p528/GetPathLoss.cpp (github.com/NTIA/p528).
[R_g, phi_g] = p528_reflection_coefficients_rev1(psi_rad, f_mhz, T_pol);
if tan(psi_rad) >= 0.1
    D_v = 1.0;
else
    r_1 = prm.D_km(1) / cos(psi_rad);                  % Eqn 8-3
    r_2 = prm.D_km(2) / cos(psi_rad);
    R_r = (r_1 * r_2) / prm.r_12_km;                   % Eqn 8-4
    term_1 = (2 * R_r * (1 + sin(psi_rad)^2)) / (prm.a_a_km * sin(psi_rad));
    term_2 = (2 * R_r / prm.a_a_km)^2;
    D_v = (1.0 + term_1 + term_2)^(-0.5);              % Eqn 8-5
end
F_r = min(prm.r_0_km / prm.r_12_km, 1);                % Eqn 8-6
R_Tg = R_g * D_v * F_r;                                % Eqn 8-7
if prm.d_km > path.d_0_km
    prm.A_LOS_db = ((prm.d_km - path.d_0_km) * (A_dML_db - A_d_0_db) / (path.d_ML_km - path.d_0_km)) + A_d_0_db;   % Eqn 8-1
else
    lambda_km = 0.2997925 / f_mhz;                     % Eqn 8-2
    if psi_rad > psi_limit
        prm.A_LOS_db = 0;                              % Step 8-2: ignore the phase lag
    else
        phi_Tg = (2 * pi * prm.delta_r_km / lambda_km) + phi_g;            % Eqn 8-8
        cplx = complex(R_Tg * cos(phi_Tg), -R_Tg * sin(phi_Tg));            % Eqn 8-9
        W_RL = min(abs(1.0 + cplx), 1.0);              % Eqn 8-10
        W_R0 = W_RL^2;                                 % Eqn 8-11
        prm.A_LOS_db = 10.0 * log10(W_R0);             % Eqn 8-12
    end
end
end
