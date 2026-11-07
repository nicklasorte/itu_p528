function [R_g, phi_g] = p528_reflection_coefficients_rev1(psi_rad, f_mhz, T_pol)
%P528_REFLECTION_COEFFICIENTS_REV1 Smooth-Earth reflection coefficient, ITU-R P.528-5 Annex 2 sec. 9.
%   [R_g, phi_g] = P528_REFLECTION_COEFFICIENTS_REV1(psi_rad, f_mhz, T_pol)
%   returns magnitude and phase for grazing angle psi_rad.
%
%   Port of src/p528/ReflectionCoefficients.cpp (github.com/NTIA/p528).
c = p528_constants_rev1();
if psi_rad <= 0.0
    sin_psi = 0.0; cos_psi = 1.0;
elseif psi_rad >= pi/2
    sin_psi = 1.0; cos_psi = 0.0;
else
    sin_psi = sin(psi_rad); cos_psi = cos(psi_rad);
end
er = c.epsilon_r;
X = (18000.0 * c.sigma) / f_mhz;                       % Eqn 9-1
Y = er - cos_psi^2;                                    % Eqn 9-2
T = sqrt(Y^2 + X^2) + Y;                               % Eqn 9-3
P = sqrt(T * 0.5);                                     % Eqn 9-4
Q = X / (2.0 * P);                                     % Eqn 9-5
if T_pol == c.POL_HORIZONTAL
    B = 1.0 / (P^2 + Q^2);                             % Eqn 9-6
    A = (2.0 * P) / (P^2 + Q^2);                       % Eqn 9-7
else
    B = (er^2 + X^2) / (P^2 + Q^2);
    A = (2.0 * (P * er + Q * X)) / (P^2 + Q^2);
end
R_g = sqrt((1.0 + (B * sin_psi^2) - (A * sin_psi)) / (1.0 + (B * sin_psi^2) + (A * sin_psi)));  % Eqn 9-8
if T_pol == c.POL_HORIZONTAL
    alpha = atan2(-Q, sin_psi - P);                    % Eqn 9-9
    beta = atan2(Q, sin_psi + P);                      % Eqn 9-10
else
    alpha = atan2((er * sin_psi) - Q, er * sin_psi - P);
    beta = atan2((X * sin_psi) + Q, er * sin_psi + P);
end
phi_g = alpha - beta;                                  % Eqn 9-11
end
