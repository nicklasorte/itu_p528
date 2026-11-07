function [T_kelvin, p_hPa, e_hPa] = p835_atmosphere_rev1(h_km)
%P835_ATMOSPHERE_REV1 ITU-R P.835 mean annual global reference atmosphere, as used by P.528-5.
%   [T_kelvin, p_hPa, e_hPa] = P835_ATMOSPHERE_REV1(h_km) returns temperature,
%   dry pressure, and water vapour partial pressure at geometric height
%   h_km (any size). Water vapour uses rho_0 = 7.5 g/m^3 with a 2 km scale
%   height, floored at a mixing ratio of 2e-6 (GlobalWetPressure in the
%   NTIA reference implementation).
%
%   Port of src/p835/MeanAnnualGlobalReferenceAtmosphere.cpp,
%   src/p835/Conversions.cpp, and src/p676/GlobalWetPressure.cpp from the
%   NTIA reference implementation (github.com/NTIA/p528). Out-of-range
%   heights return the reference implementation's error codes (-1 below
%   0 km, -2 above 100 km) in place of values, as the C++ does.
arguments
    h_km double
end
T_kelvin = global_temperature(h_km);
p_hPa = global_pressure(h_km);
rho_0 = 7.5;                                            % RHO_0__M_KG
rho = rho_0 * exp(-h_km / 2);                            % GlobalWaterVapourDensity, 2 km scale height
rho(h_km < 0) = -1; rho(h_km > 100) = -2;
rho = max(rho, 2e-6 * 216.7 * p_hPa ./ T_kelvin);
e_hPa = (rho .* T_kelvin) / 216.7;                       % WaterVapourDensityToPressure
end

function T = global_temperature(h)
T = zeros(size(h));
low = h < 86;
hp = (6356.766 * h) ./ (6356.766 + h);                   % ConvertToGeopotentialHeight
T(low) = temp_regime1(hp(low));
T(~low) = temp_regime2(h(~low));
T(h < 0) = -1;
T(h > 100) = -2;
end

function T = temp_regime1(hp)
T = -2 * ones(size(hp));
T(hp < 0) = -1;
m = hp >= 0 & hp <= 11;           T(m) = 288.15 - 6.5 * hp(m);
m = hp > 11 & hp <= 20;           T(m) = 216.65;
m = hp > 20 & hp <= 32;           T(m) = 216.65 + (hp(m) - 20);
m = hp > 32 & hp <= 47;           T(m) = 228.65 + 2.8 * (hp(m) - 32);
m = hp > 47 & hp <= 51;           T(m) = 270.65;
m = hp > 51 & hp <= 71;           T(m) = 270.65 - 2.8 * (hp(m) - 51);
m = hp > 71 & hp <= 84.852;       T(m) = 214.65 - 2.0 * (hp(m) - 71);
end

function T = temp_regime2(h)
T = -2 * ones(size(h));
T(h < 86) = -1;
m = h >= 86 & h <= 91;            T(m) = 186.8673;
m = h > 91 & h <= 100;            T(m) = 263.1905 - 76.3232 * sqrt(1 - ((h(m) - 91) / 19.9429).^2);
end

function p = global_pressure(h)
p = zeros(size(h));
low = h < 86;
hp = (6356.766 * h) ./ (6356.766 + h);
p(low) = press_regime1(hp(low));
hh = h(~low);
a = [95.571899, -4.011801, 6.424731e-2, -4.789660e-4, 1.340543e-6];
p(~low) = exp(a(1) + a(2) * hh + a(3) * hh.^2 + a(4) * hh.^3 + a(5) * hh.^4);
p(h < 0) = -1;
p(h > 100) = -2;
end

function p = press_regime1(hp)
p = -2 * ones(size(hp));
p(hp < 0) = -1;
m = hp >= 0 & hp <= 11;      p(m) = 1013.25 * (288.15 ./ (288.15 - 6.5 * hp(m))).^(-34.1632 / 6.5);
m = hp > 11 & hp <= 20;      p(m) = 226.3226 * exp(-34.1632 * (hp(m) - 11) / 216.65);
m = hp > 20 & hp <= 32;      p(m) = 54.74980 * (216.65 ./ (216.65 + (hp(m) - 20))).^34.1632;
m = hp > 32 & hp <= 47;      p(m) = 8.680422 * (228.65 ./ (228.65 + 2.8 * (hp(m) - 32))).^(34.1632 / 2.8);
m = hp > 47 & hp <= 51;      p(m) = 1.109106 * exp(-34.1632 * (hp(m) - 47) / 270.65);
m = hp > 51 & hp <= 71;      p(m) = 0.6694167 * (270.65 ./ (270.65 - 2.8 * (hp(m) - 51))).^(-34.1632 / 2.8);
m = hp > 71 & hp <= 84.852;  p(m) = 0.03956649 * (214.65 ./ (214.65 - 2.0 * (hp(m) - 71))).^(-34.1632 / 2.0);
end
