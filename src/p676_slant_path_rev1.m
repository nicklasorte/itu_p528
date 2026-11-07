function res = p676_slant_path_rev1(f_ghz, h_1_km, h_2_km, beta_1_rad)
%P676_SLANT_PATH_REV1 ITU-R P.676 Annex 1 slant-path gaseous attenuation by ray tracing.
%   res = P676_SLANT_PATH_REV1(f_ghz, h_1_km, h_2_km, beta_1_rad) traces a
%   ray from height h_1_km to h_2_km starting at incident angle beta_1_rad
%   (measured from zenith) through the P.835 global reference atmosphere.
%   For beta_1_rad > pi/2 (negative elevation) it finds the grazing height
%   and traces both legs (P.676 Annex 1, sec. 2.2.2).
%
%   res fields: A_gas_db, bending_rad, a_km (ray length), angle_rad
%   (incident angle at the top), delta_L_km (excess path length).
%
%   Port of src/p676/SlantPath.cpp and src/p676/RayTrace.cpp from the NTIA
%   reference implementation (github.com/NTIA/p528).
arguments
    f_ghz (1,1) double
    h_1_km (1,1) double
    h_2_km (1,1) double
    beta_1_rad (1,1) double
end
a_0_km = 6371.0;
if beta_1_rad > pi/2
    [T, p, e] = p835_atmosphere_rev1(h_1_km);
    n_1 = refractive_index(p, T, e);
    h_G_km = h_1_km;
    delta = h_1_km / 2;
    diff = 100;
    for iter = 1:200
        if iter == 200
            % The C++ loop is unbounded and hangs when no grazing point exists above
            % the surface (ray steeper than the horizon dip); fail loud instead.
            error('p676:noGrazingPoint', ...
                'No grazing height found for h_1 = %.6g km, beta_1 = %.10g rad.', h_1_km, beta_1_rad);
        end
        if diff > 0
            h_G_km = h_G_km - delta;
        else
            h_G_km = h_G_km + delta;
        end
        delta = delta / 2;
        [T, p, e] = p835_atmosphere_rev1(h_G_km);
        n_G = refractive_index(p, T, e);
        diff = n_G * (a_0_km + h_G_km) - n_1 * (a_0_km + h_1_km) * sin(beta_1_rad);
        if ~(abs(diff) > 0.001), break; end
    end
    r1 = ray_trace(f_ghz, h_G_km, h_1_km, pi/2);
    r2 = ray_trace(f_ghz, h_G_km, h_2_km, pi/2);
    res.A_gas_db = r1.A_gas_db + r2.A_gas_db;
    res.bending_rad = r1.bending_rad + r2.bending_rad;
    res.a_km = r1.a_km + r2.a_km;
    res.angle_rad = r2.angle_rad;
    res.delta_L_km = r1.delta_L_km + r2.delta_L_km;
else
    res = ray_trace(f_ghz, h_1_km, h_2_km, beta_1_rad);
end
end

function res = ray_trace(f_ghz, h_1_km, h_2_km, beta_1_rad)
% RayTrace.cpp, vectorized over layers. Layer properties do not depend on
% the ray, and Snell's law (Eqn 19b) uses only the bottom layer, so every
% layer can be evaluated at once; the sums match the C++ loop order of terms.
a_0_km = 6371.0;
e1 = exp(1/100);
i_lower = floor(100 * log(1e4 * h_1_km * (e1 - 1) + 1) + 1);          % Eqn 16(a)
i_upper = ceil(100 * log(1e4 * h_2_km * (e1 - 1) + 1) + 1);           % Eqn 16(b)
m = ((exp(2/100) - e1) / (exp(i_upper/100) - exp(i_lower/100))) * (h_2_km - h_1_km);   % Eqn 16(c)

res = struct('A_gas_db', 0, 'bending_rad', 0, 'a_km', 0, 'angle_rad', beta_1_rad, 'delta_L_km', 0);
if i_upper - 1 < i_lower
    return
end
i = (i_lower:i_upper)';                                               % layers i_lower..i_upper
h = h_1_km + m * ((exp((i - 1)/100) - exp((i_lower - 1)/100)) / (e1 - 1));
delta = m * exp((i - 1)/100);                                         % LayerThickness, Eqn 14
[T, p, e] = p835_atmosphere_rev1(h + delta/2);
[gamma, n] = p676_specific_attenuation_rev1(f_ghz, T, e, p);
r = a_0_km + h;

k = 1:(numel(i) - 1);                                                 % loop i = i_lower .. i_upper-1
n_1 = n(1); r_1 = r(1);
beta_i = asin(min(1, (n_1 * r_1) ./ (n(k) .* r(k)) * sin(beta_1_rad)));       % Eqn 19b
alpha_i = asin(min(1, (n_1 * r_1) ./ (n(k) .* r(k + 1)) * sin(beta_1_rad)));  % Eqn 18a
a_i = -r(k) .* cos(beta_i) + sqrt(r(k).^2 .* cos(beta_i).^2 + 2 * r(k) .* delta(k) + delta(k).^2);  % Eqn 17
arg = n(k) ./ n(k + 1) .* sin(alpha_i);
arg(abs(arg) > 1) = NaN;                                              % C++ asin returns NaN outside [-1,1]
beta_ii = asin(arg);

res.a_km = sum_seq(a_i);
res.A_gas_db = sum_seq(a_i .* gamma(k));
res.delta_L_km = sum_seq(a_i .* (n(k) - 1));
kb = 1:(numel(k) - 1);                                                % Eqn 22a: summation stops at i_max - 1
res.bending_rad = sum_seq(beta_ii(kb) - alpha_i(kb));
res.angle_rad = alpha_i(end);
end

function s = sum_seq(x)
% Sequential left-to-right sum, matching the C++ accumulation order.
s = 0;
for j = 1:numel(x)
    s = s + x(j);
end
end

function n = refractive_index(p, T, e)
N = 77.6 * p / T + 72 * e / T + 3.75e5 * e / T^2;
n = 1 + N * 1e-6;
end
