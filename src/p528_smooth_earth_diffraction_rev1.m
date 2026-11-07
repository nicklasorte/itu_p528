function A_db = p528_smooth_earth_diffraction_rev1(d_1_km, d_2_km, f_mhz, d_0_km, T_pol)
%P528_SMOOTH_EARTH_DIFFRACTION_REV1 Smooth-Earth diffraction loss, ITU-R P.528-5 Annex 2 sec. 10.
%   A_db = P528_SMOOTH_EARTH_DIFFRACTION_REV1(d_1_km, d_2_km, f_mhz, d_0_km, T_pol)
%   uses Vogler (1964) with the IF-77 height-function blend.
%
%   Port of src/p528/SmoothEarthDiffraction.cpp (github.com/NTIA/p528).
c = p528_constants_rev1();
s = 18000 * c.sigma / f_mhz;
if T_pol == c.POL_HORIZONTAL
    K = 0.01778 * f_mhz^(-c.third) * ((c.epsilon_r - 1)^2 + s^2)^(-0.25);
else
    K = 0.01778 * f_mhz^(-c.third) * ((c.epsilon_r^2 + s^2) / ((c.epsilon_r - 1)^2 + s^2)^0.5)^0.5;
end
B_0 = 1.607;
x_0_km = (B_0 - K) * f_mhz^c.third * d_0_km;      % [Vogler 1964, Eqn 2], C_0 = 1
x_1_km = (B_0 - K) * f_mhz^c.third * d_1_km;
x_2_km = (B_0 - K) * f_mhz^c.third * d_2_km;
G_x_db = distance_function(x_0_km);
F_x1_db = height_function(x_1_km, K);
F_x2_db = height_function(x_2_km, K);
A_db = G_x_db - F_x1_db - F_x2_db - 20.0;          % [Vogler 1964, Eqn 1], C_1 = 20
end

function G = distance_function(x)
G = 0.05751 * x - 10.0 * log10(x);                 % [Vogler 1964, Eqn 13]
end

function F = height_function(x, K)
y_db = 40.0 * log10(x) - 117.0;                    % [FAA-ES-83-3, Eqn 73]
G_db = distance_function(x);
if x <= 200.0
    x_t = 450 / -(log10(K)^3);                     % Eqn 109
    if x >= x_t                                    % Eqn 110
        if abs(y_db) < 117
            F = y_db;
        else
            F = -117;
        end
    else
        F = 20 * log10(K) - 15 + (0.000025 * x^2 / K);
    end
elseif x > 2000.0
    F = G_db;                                      % [Vogler 1964] F ~= G for large x
else
    W = 0.0134 * x * exp(-0.005 * x);              % [FAA-ES-83-3, Eqn 72]
    F = W * y_db + (1.0 - W) * G_db;               % [FAA-ES-83-3, Eqn 75]
end
end
