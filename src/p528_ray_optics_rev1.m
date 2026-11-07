function prm = p528_ray_optics_rev1(t1, t2, psi)
%P528_RAY_OPTICS_REV1 Two-ray geometry for grazing angle psi, ITU-R P.528-5 Annex 2 sec. 7.
%   prm = P528_RAY_OPTICS_REV1(t1, t2, psi) returns the line-of-sight
%   parameter struct: z_km, theta, D_km (1x2), d_km, r_0_km, r_12_km,
%   delta_r_km, theta_h1_rad, theta_h2_rad, a_a_km, A_LOS_db.
%
%   Port of src/p528/RayOptics.cpp (github.com/NTIA/p528).
a_0 = 6371.0; a_e = 9257.0;
z = (a_0 / a_e) - 1;                                   % Eqn 7-1
k_a = 1 / (1 + z * cos(psi));                          % Eqn 7-2
prm.a_a_km = a_0 * k_a;                                % Eqn 7-3
dh_a1 = t1.delta_h_km * (prm.a_a_km - a_0) / (a_e - a_0);   % Eqn 7-4
dh_a2 = t2.delta_h_km * (prm.a_a_km - a_0) / (a_e - a_0);
H = [t1.h_r_km - dh_a1, t2.h_r_km - dh_a2];            % Eqn 7-5
prm.z_km = prm.a_a_km + H;                             % Eqn 7-6
prm.theta = acos(prm.a_a_km * cos(psi) ./ prm.z_km) - psi;   % Eqn 7-7
prm.D_km = prm.z_km .* sin(prm.theta);                 % Eqn 7-8
if psi > 1.56                                          % Eqn 7-9
    Hp = H;
else
    Hp = prm.D_km * tan(psi);
end
delta_z = abs(prm.z_km(1) - prm.z_km(2));              % Eqn 7-10
prm.d_km = max(prm.a_a_km * (prm.theta(1) + prm.theta(2)), 0);   % Eqn 7-11
alpha = atan((Hp(2) - Hp(1)) / (prm.D_km(1) + prm.D_km(2)));      % Eqn 7-12
prm.r_0_km = max(delta_z, (prm.D_km(1) + prm.D_km(2)) / cos(alpha));   % Eqn 7-13
prm.r_12_km = (prm.D_km(1) + prm.D_km(2)) / cos(psi);  % Eqn 7-14
prm.delta_r_km = 4.0 * Hp(1) * Hp(2) / (prm.r_0_km + prm.r_12_km);   % Eqn 7-15
prm.theta_h1_rad = alpha - prm.theta(1);               % Eqn 7-16
prm.theta_h2_rad = -(alpha + prm.theta(2));            % Eqn 7-17
prm.A_LOS_db = 0;
end
