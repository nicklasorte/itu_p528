function t = p528_terminal_geometry_rev1(f_mhz, h_r_km)
%P528_TERMINAL_GEOMETRY_REV1 Terminal geometry, ITU-R P.528-5 Annex 2 sec. 4.
%   t = P528_TERMINAL_GEOMETRY_REV1(f_mhz, h_r_km) ray-traces a horizontal
%   ray from the ground to the terminal height and returns the terminal
%   struct: h_r_km, theta_rad, A_a_db, a_km, d_r_km, phi_rad, h_e_km,
%   delta_h_km.
%
%   Port of src/p528/TerminalGeometry.cpp (github.com/NTIA/p528).
c = p528_constants_rev1();
theta_tx_rad = 0;
r = p676_slant_path_rev1(f_mhz / 1000, 0, h_r_km, pi/2 - theta_tx_rad);
t.h_r_km = h_r_km;
t.theta_rad = pi/2 - r.angle_rad;
t.A_a_db = r.A_gas_db;
t.a_km = r.a_km;
central_angle = (pi/2 - r.angle_rad) - theta_tx_rad + r.bending_rad;   % [Thayer, Eqn 2], rearranged
t.d_r_km = c.a_0_km * central_angle;
t.phi_rad = t.d_r_km / c.a_e_km;                                       % Eqn 4-1
t.h_e_km = (c.a_e_km / cos(t.phi_rad)) - c.a_e_km;                     % Eqn 4-2
t.delta_h_km = t.h_r_km - t.h_e_km;                                    % Eqn 4-3
end
