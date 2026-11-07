function out = p528_variability_rev1(fn, varargin)
%P528_VARIABILITY_REV1 Variability helpers of ITU-R P.528-5 Annex 2 secs. 13-14.
%   out = P528_VARIABILITY_REV1(fn, ...) dispatches to one helper:
%
%   "long_term"   [Y_e_db, A_Y] = (d_r1_km, d_r2_km, d_km, f_mhz, p, f_theta_h, A_T)
%                 LongTermVariability.cpp; returns a 1x2 vector [Y_e_db, A_Y].
%   "nakagami_rice"   Y = (K, p)              NakagamiRice.cpp
%   "find_k_99"       K = (Y_pi_99_db)        FindKForYpiAt99Percent.cpp
%   "combine"         C = (A_M, A_p, B_M, B_p, p)   CombineDistributions.cpp
%   "inv_ccdf"        Q = (q)                 InverseComplementaryCumulativeDistributionFunction.cpp
%
%   Ports from the NTIA reference implementation (github.com/NTIA/p528).
%   NakagamiRice reproduces the reference code exactly, including the
%   expression NakagamiRiceCurves[d_K - 1][d_p] - 1 in the second
%   interpolation, which subtracts 1 dB where a column index (d_p - 1)
%   appears to be intended. See README, "Known differences".
switch fn
    case "long_term"
        [Y_e, A_Y] = long_term(varargin{:});
        out = [Y_e, A_Y];
    case "nakagami_rice"
        out = nakagami_rice(varargin{:});
    case "find_k_99"
        out = find_k_99(varargin{:});
    case "combine"
        out = combine(varargin{:});
    case "inv_ccdf"
        out = inv_ccdf(varargin{:});
    otherwise
        error('p528_variability_rev1:fn', 'Unknown helper "%s".', fn);
end
end

function [Y_e_db, A_Y] = long_term(d_r1_km, d_r2_km, d_km, f_mhz, p, f_theta_h, A_T)
c = p528_constants_rev1();
d_qs_km = 65.0 * (100.0 / f_mhz)^c.third;          % Eqn 14-1
d_Lq_km = d_r1_km + d_r2_km;                       % Eqn 14-2
d_q_km = d_Lq_km + d_qs_km;                        % Eqn 14-3
if d_km <= d_q_km                                  % Eqn 14-4
    d_e_km = (130.0 * d_km) / d_q_km;
else
    d_e_km = 130.0 + d_km - d_q_km;
end
if f_mhz > 1600.0                                  % Eqns 14-5, 14-6
    g_10 = 1.05; g_90 = 1.05;
else
    g_10 = (0.21 * sin(5.22 * log10(f_mhz / 200.0))) + 1.28;
    g_90 = (0.18 * sin(5.22 * log10(f_mhz / 200.0))) + 1.23;
end
% Tech Note 101, Vol 2: Table III.4 row A*, Table III.3 row A*, Table III.5 continental temperate
c_1 = [2.93e-4, 5.25e-4, 1.59e-5];
c_2 = [3.78e-8, 1.57e-6, 1.56e-11];
c_3 = [1.02e-7, 4.70e-7, 2.77e-8];
n_1 = [2.00, 1.97, 2.32];
n_2 = [2.88, 2.31, 4.08];
n_3 = [3.15, 2.90, 3.25];
f_inf = [3.2, 5.4, 0.0];
f_m = [8.2, 10.0, 3.9];
Z = zeros(1, 3);                                   % [Y_0(90) Y_0(10) V(50)]
for i = 1:3
    f_2 = f_inf(i) + ((f_m(i) - f_inf(i)) * exp(-c_2(i) * d_e_km^n_2(i)));
    Z(i) = (c_1(i) * d_e_km^n_1(i) - f_2) * exp(-c_3(i) * d_e_km^n_3(i)) + f_2;
end
if p == 50
    Y_p_db = Z(3);
elseif p > 50
    z_90 = inv_ccdf(90.0 / 100.0);
    z_p = inv_ccdf(p / 100.0);
    c_p = z_p / z_90;
    Y_p_db = c_p * (-Z(1) * g_90) + Z(3);
else
    if p >= 10
        z_10 = inv_ccdf(10.0 / 100.0);
        z_p = inv_ccdf(p / 100.0);
        c_p = z_p / z_10;
    else
        ps = [1, 2, 5, 10];                        % [15], Table 10, p. 34, Climate 6
        c_ps = [1.9507, 1.7166, 1.3265, 1.0000];
        dist = upper_bound(c.P, p);
        c_p = lin_interp(ps(dist), c_ps(dist), ps(dist + 1), c_ps(dist + 1), p);
    end
    Y_p_db = c_p * (Z(2) * g_10) + Z(3);
end
Y_10_db = (Z(2) * g_10) + Z(3);                    % Eqn 14-20
Y_eI_db = f_theta_h * Y_p_db;                      % Eqn 14-21
Y_eI_10_db = f_theta_h * Y_10_db;                  % Eqn 14-22
A_YI = (A_T + Y_eI_10_db) - 3.0;                   % Eqn 14-23
A_Y = max(A_YI, 0);                                % Eqn 14-24
Y_e_db = Y_eI_db - A_Y;                            % Eqn 14-25
if p < 10                                          % [Gierhart 1970] free-space cap
    c_Y = [-5.0, -4.5, -3.7, 0.0];
    dist = upper_bound(c.P, p);
    c_Yi = lin_interp(c.P(dist), c_Y(dist), c.P(dist + 1), c_Y(dist + 1), p);
    Y_e_db = Y_e_db + A_T;
    if Y_e_db > -c_Yi
        Y_e_db = -c_Yi;
    end
    Y_e_db = Y_e_db - A_T;
end
end

function Y = nakagami_rice(K, p)
% Indices below are C++ 0-based distances; arrays are read at (index + 1).
c = p528_constants_rev1();
d_K = lower_bound(c.K, K);
d_p = lower_bound(c.P, p);
C = c.NakagamiRiceCurves;
nK = numel(c.K);
if d_K == 0                                        % K <= -40
    if d_p == 0
        Y = C(1, 1);
    else
        Y = lin_interp(c.P(d_p + 1), C(1, d_p + 1), c.P(d_p), C(1, d_p), p);
    end
elseif d_K == nK                                   % K > 20
    if d_p == 0
        Y = C(d_K, 1);
    else
        Y = lin_interp(c.P(d_p + 1), C(d_K, d_p + 1), c.P(d_p), C(d_K, d_p), p);
    end
else
    if d_p == 0
        Y = lin_interp(c.K(d_K + 1), C(d_K + 1, 1), c.K(d_K), C(d_K, 1), K);
    else
        v1 = lin_interp(c.K(d_K + 1), C(d_K + 1, d_p + 1), c.K(d_K), C(d_K, d_p + 1), K);
        % Reference code: NakagamiRiceCurves[d_K - 1][d_p] - 1 (reproduced as written)
        v2 = lin_interp(c.K(d_K + 1), C(d_K + 1, d_p), c.K(d_K), C(d_K, d_p + 1) - 1, K);
        Y = lin_interp(c.P(d_p + 1), v1, c.P(d_p), v2, p);
    end
end
end

function K = find_k_99(Y_pi_99_db)
c = p528_constants_rev1();
col = c.NakagamiRiceCurves(:, c.Y_pi_99_index);
if Y_pi_99_db < col(1)
    K = c.K(1);
    return
end
for i = 1:numel(c.K)
    if Y_pi_99_db - col(i) < 0
        K = (c.K(i) * (Y_pi_99_db - col(i - 1)) - c.K(i - 1) * (Y_pi_99_db - col(i))) / (col(i) - col(i - 1));
        return
    end
end
K = c.K(end);
end

function C_out = combine(A_M, A_p, B_M, B_p, p)
C_M = A_M + B_M;
Y_3 = sqrt((A_p - A_M)^2 + (B_p - B_M)^2);
if p < 50
    C_out = C_M + Y_3;
else
    C_out = C_M - Y_3;
end
end

function Q_q = inv_ccdf(q)
% Abramowitz and Stegun 26.2.23 rational approximation, as in the reference code.
C_0 = 2.515516; C_1 = 0.802853; C_2 = 0.010328;
D_1 = 1.432788; D_2 = 0.189269; D_3 = 0.001308;
x = q;
if q > 0.5
    x = 1.0 - x;
end
T_x = sqrt(-2.0 * log(x));
zeta_x = ((C_2 * T_x + C_1) * T_x + C_0) / (((D_3 * T_x + D_2) * T_x + D_1) * T_x + 1.0);
Q_q = T_x - zeta_x;
if q > 0.5
    Q_q = -Q_q;
end
end

function y = lin_interp(x1, y1, x2, y2, x)
y = (y1 * (x2 - x) + y2 * (x - x1)) / (x2 - x1);
end

function d = lower_bound(v, x)
% std::lower_bound distance: count of elements < x.
d = find(v >= x, 1) - 1;
if isempty(d), d = numel(v); end
end

function d = upper_bound(v, x)
% std::upper_bound distance: count of elements <= x.
d = find(v > x, 1) - 1;
if isempty(d), d = numel(v); end
end
