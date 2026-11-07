function [gamma_db_km, n] = p676_specific_attenuation_rev1(f_ghz, T_kelvin, e_hPa, p_hPa)
%P676_SPECIFIC_ATTENUATION_REV1 ITU-R P.676 Annex 1 gaseous specific attenuation and refractive index.
%   [gamma_db_km, n] = P676_SPECIFIC_ATTENUATION_REV1(f_ghz, T_kelvin, e_hPa, p_hPa)
%   returns specific attenuation (dB/km) and refractive index for column
%   vectors of temperature, water vapour pressure, and dry pressure at one
%   frequency.
%
%   Port of SpecificAttenuation.cpp, Refractivity.cpp, LineShapeFactor.cpp,
%   NonresonantDebyeAttenuation.cpp, and RefractiveIndex.cpp from the NTIA
%   reference implementation (github.com/NTIA/p528). Vectorized over
%   heights (rows) and spectral lines (columns); the C++ loops over lines.
arguments
    f_ghz (1,1) double
    T_kelvin (:,1) double
    e_hPa (:,1) double
    p_hPa (:,1) double
end
[ox, wv] = p676_line_data_rev1();
theta = 300 ./ T_kelvin;

% Oxygen, line-by-line (implicit expansion: heights x lines)
S_o = ox.a_1 * 1e-7 .* p_hPa .* theta.^3 .* exp(ox.a_2 .* (1 - theta));
df_o = ox.a_3 * 1e-4 .* (p_hPa .* theta.^(0.8 - ox.a_4) + 1.1 * e_hPa .* theta);
df_o = sqrt(df_o.^2 + 2.25e-6);
dl_o = (ox.a_5 + ox.a_6 .* theta) * 1e-4 .* (p_hPa + e_hPa) .* theta.^0.8;
N_o = sum(S_o .* line_shape(f_ghz, ox.f_0, df_o, dl_o), 2);

% Non-resonant Debye term
d = 5.6e-4 * (p_hPa + e_hPa) .* theta.^0.8;
frac_1 = 6.14e-5 ./ (d .* (1 + (f_ghz ./ d).^2));
frac_2 = (1.4e-12 * p_hPa .* theta.^1.5) / (1 + 1.9e-5 * f_ghz^1.5);
N_D = f_ghz * p_hPa .* theta.^2 .* (frac_1 + frac_2);

% Water vapour
S_w = 0.1 * wv.b_1 .* e_hPa .* theta.^3.5 .* exp(wv.b_2 .* (1 - theta));
df_w = 1e-4 * wv.b_3 .* (p_hPa .* theta.^wv.b_4 + wv.b_5 .* e_hPa .* theta.^wv.b_6);
term1 = 0.217 * df_w.^2 + (2.1316e-12 * wv.f_0.^2 ./ theta);
df_w = 0.535 * df_w + sqrt(term1);
N_w = sum(S_w .* line_shape(f_ghz, wv.f_0, df_w, 0), 2);

gamma_db_km = 0.1820 * f_ghz * (N_o + N_D) + 0.1820 * f_ghz * N_w;

N_dry = 77.6 * p_hPa ./ T_kelvin;
N_wet = 72 * e_hPa ./ T_kelvin + 3.75e5 * e_hPa ./ T_kelvin.^2;
n = 1 + (N_dry + N_wet) * 1e-6;
end

function F = line_shape(f, f_i, df, delta)
term1 = f ./ f_i;
term2 = (df - delta .* (f_i - f)) ./ ((f_i - f).^2 + df.^2);
term3 = (df - delta .* (f_i + f)) ./ ((f_i + f).^2 + df.^2);
F = term1 .* (term2 + term3);
end
