function tr = p528_troposcatter_rev1(t1, t2, d_km, f_mhz)
%P528_TROPOSCATTER_REV1 Troposcatter loss, ITU-R P.528-5 Annex 2 sec. 11.
%   tr = P528_TROPOSCATTER_REV1(t1, t2, d_km, f_mhz) returns a struct with
%   d_s_km, d_z_km, h_v_km, theta_s, theta_A, A_s_db for terminals t1, t2
%   (from P528_TERMINAL_GEOMETRY_REV1).
%
%   Port of src/p528/Troposcatter.cpp (github.com/NTIA/p528).
c = p528_constants_rev1();
a_e = c.a_e_km;
N_s = c.N_s;
tr.d_s_km = d_km - t1.d_r_km - t2.d_r_km;                          % Eqn 11-2
if tr.d_s_km <= 0.0
    tr.d_z_km = 0.0; tr.A_s_db = 0.0; tr.d_s_km = 0.0;
    tr.h_v_km = 0.0; tr.theta_s = 0.0; tr.theta_A = 0.0;
    return
end
tr.d_z_km = 0.5 * tr.d_s_km;                                       % Eqn 11-6
A_m = 1 / c.a_0_km;                                                % Eqn 11-7
dN = A_m - (1.0 / a_e);                                            % Eqn 11-8
gamma_e_km = (N_s * 1e-6) / dN;                                    % Eqn 11-9
z_a_km = 1.0 / (2 * a_e) * (tr.d_z_km / 2)^2;                      % Eqn 11-10
z_b_km = 1.0 / (2 * a_e) * tr.d_z_km^2;                            % Eqn 11-11
Q_o = A_m - dN;                                                    % Eqn 11-12
Q_a = A_m - dN / exp(min(35.0, z_a_km / gamma_e_km));              % Eqn 11-13
Q_b = A_m - dN / exp(min(35.0, z_b_km / gamma_e_km));
Z_a_km = (7.0 * Q_o + 6.0 * Q_a - Q_b) * (tr.d_z_km^2 / 96.0);     % Eqn 11-14
Z_b_km = (Q_o + 2.0 * Q_a) * (tr.d_z_km^2 / 6.0);                  % Eqn 11-15
Q_A = A_m - dN / exp(min(35.0, Z_a_km / gamma_e_km));              % Eqn 11-16
Q_B = A_m - dN / exp(min(35.0, Z_b_km / gamma_e_km));
tr.h_v_km = (Q_o + 2.0 * Q_A) * (tr.d_z_km^2 / 6.0);               % Eqn 11-17
tr.theta_A = (Q_o + 4.0 * Q_A + Q_B) * tr.d_z_km / 6.0;            % Eqn 11-18
tr.theta_s = 2 * tr.theta_A;                                       % Eqn 11-19

% Scattering efficiency
epsilon_1 = 5.67e-6 * N_s^2 - 0.00232 * N_s + 0.031;               % Eqn 11-20
epsilon_2 = 0.0002 * N_s^2 - 0.06 * N_s + 6.6;                     % Eqn 11-21
gamma = 0.1424 * (1.0 + epsilon_1 / exp(min(35.0, (tr.h_v_km / 4.0)^6)));      % Eqn 11-22
S_e_db = 83.1 - epsilon_2 / (1.0 + 0.07716 * tr.h_v_km^2) ...
    + 20 * log10((0.1424 / gamma)^2 * exp(gamma * tr.h_v_km));     % Eqn 11-23

% Scattering volume
X_A1 = t1.h_e_km^2 + 4.0 * (a_e + t1.h_e_km) * a_e * sin(t1.d_r_km / (a_e * 2))^2;  % Eqn 11-24
X_A2 = t2.h_e_km^2 + 4.0 * (a_e + t2.h_e_km) * a_e * sin(t2.d_r_km / (a_e * 2))^2;
ell_1 = sqrt(X_A1) + tr.d_z_km;                                    % Eqn 11-25
ell_2 = sqrt(X_A2) + tr.d_z_km;
ell = ell_1 + ell_2;                                               % Eqn 11-26
s = (ell_1 - ell_2) / ell;                                         % Eqn 11-27
eta = gamma * tr.theta_s * ell / 2;                                % Eqn 11-28
kappa = f_mhz / 0.0477;                                            % Eqn 11-29
rho_1 = 2.0 * kappa * tr.theta_s * t1.h_e_km;                      % Eqn 11-30
rho_2 = 2.0 * kappa * tr.theta_s * t2.h_e_km;
SQRT2 = sqrt(2);
A = (1 - s^2)^2;                                                   % Eqn 11-36
X_v1 = (1 + s)^2 * eta;                                            % Eqn 11-32
X_v2 = (1 - s)^2 * eta;                                            % Eqn 11-33
q_1 = X_v1^2 + rho_1^2;                                            % Eqn 11-34
q_2 = X_v2^2 + rho_2^2;                                            % Eqn 11-35
B_s = 6 + 8 * s^2 ...
    + 8 * (1.0 - s) * X_v1^2 * rho_1^2 / q_1^2 ...
    + 8 * (1.0 + s) * X_v2^2 * rho_2^2 / q_2^2 ...
    + 2 * (1.0 - s^2) * (1 + 2 * X_v1^2 / q_1) * (1 + 2 * X_v2^2 / q_2);         % Eqn 11-37
C_s = 12 * ((rho_1 + SQRT2) / rho_1)^2 * ((rho_2 + SQRT2) / rho_2)^2 ...
    * (rho_1 + rho_2) / (rho_1 + rho_2 + 2 * SQRT2);               % Eqn 11-38
temp = (A * eta^2 + B_s * eta) * q_1 * q_2 / (rho_1^2 * rho_2^2);
S_v_db = 10 * log10(temp + C_s);

tr.A_s_db = S_e_db + S_v_db + 10.0 * log10(kappa * tr.theta_s^3 / ell);
end
