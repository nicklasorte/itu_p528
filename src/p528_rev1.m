function [A_db, res, rtn] = p528_rev1(d_km, h_1_m, h_2_m, f_mhz, T_pol, p, opts)
%P528_REV1 ITU-R P.528-5 basic transmission loss for aeronautical paths (MATLAB port of the U.S. reference).
%   [A_db, res, rtn] = P528_REV1(d_km, h_1_m, h_2_m, f_mhz, T_pol, p)
%   returns the basic transmission loss A_db (dB) not exceeded for p% of
%   time, following Annex 2 (step-by-step method) of Recommendation
%   ITU-R P.528-5.
%
%   Inputs (scalars or arrays of a common size; scalars expand):
%     d_km    great-circle distance between terminals, km (>= 0)
%     h_1_m   height of the low terminal, m (1.5-20000; warnings to 80000)
%     h_2_m   height of the high terminal, m (h_1_m <= h_2_m)
%     f_mhz   frequency, MHz (100-30000)
%     T_pol   polarization: 0 horizontal, 1 vertical
%     p       time percentage (1-99)
%
%   Name-value options:
%     error_mode     "error" (default) raises p528:invalidInput for any
%                    invalid input; "code" returns NaN in A_db and the
%                    reference error code in rtn.
%     parallel_flag  0 (default) serial; 1 uses parfor over elements.
%
%   Outputs:
%     A_db  basic transmission loss, dB (size of the expanded inputs)
%     res   struct array with fields d_km, A_db, A_fs_db, A_a_db,
%           theta_h1_rad, propagation_mode (1 LOS, 2 diffraction,
%           3 troposcatter), and warnings (bit flags)
%     rtn   return codes of the reference implementation (0 success,
%           11 success with warnings, 1-10 validation errors)
%
%   Example (README test case 1):
%     A = p528_rev1(15, 10, 1000, 500, 0, 50)   % 110.0 dB
%
%   This is a derivative of the NTIA U.S. reference implementation of
%   ITU-R P.528-5 (github.com/NTIA/p528, commit 74870ab), which is not
%   subject to copyright in the United States. See README for changes.
arguments
    d_km double
    h_1_m double
    h_2_m double
    f_mhz double
    T_pol double
    p double
    opts.error_mode (1,1) string {mustBeMember(opts.error_mode, ["error", "code"])} = "error"
    opts.parallel_flag (1,1) double {mustBeMember(opts.parallel_flag, [0 1])} = 0
end
[d_km, h_1_m, h_2_m, f_mhz, T_pol, p] = expand_inputs(d_km, h_1_m, h_2_m, f_mhz, T_pol, p);
n = numel(d_km);
A_db = NaN(size(d_km));
rtn = zeros(size(d_km));
template = struct('d_km', 0, 'A_db', 0, 'A_fs_db', 0, 'A_a_db', 0, 'theta_h1_rad', 0, ...
    'propagation_mode', 0, 'warnings', 0);
res = repmat(template, size(d_km));

if opts.parallel_flag == 1
    res_lin = repmat(template, 1, n);
    rtn_lin = zeros(1, n);
    parfor k = 1:n
        [res_lin(k), rtn_lin(k)] = p528_single_rev1(d_km(k), h_1_m(k), h_2_m(k), f_mhz(k), T_pol(k), p(k));
    end
    res(:) = res_lin;
    rtn(:) = rtn_lin;
else
    for k = 1:n
        [res(k), rtn(k)] = p528_single_rev1(d_km(k), h_1_m(k), h_2_m(k), f_mhz(k), T_pol(k), p(k));
    end
end

ok = rtn == 0 | rtn == 11;
A_db(ok) = [res(ok).A_db];
if opts.error_mode == "error" && any(~ok)
    k = find(~ok, 1);
    error('p528:invalidInput', ...
        'Element %d: reference error code %d (%s). Inputs d=%.15g km, h1=%.15g m, h2=%.15g m, f=%.15g MHz, pol=%d, p=%.15g.', ...
        k, rtn(k), error_text(rtn(k)), d_km(k), h_1_m(k), h_2_m(k), f_mhz(k), T_pol(k), p(k));
end
end

function varargout = expand_inputs(varargin)
sz = [1 1];
for k = 1:nargin
    if ~isscalar(varargin{k})
        if isequal(sz, [1 1])
            sz = size(varargin{k});
        elseif ~isequal(size(varargin{k}), sz)
            error('p528:sizeMismatch', 'Non-scalar inputs must all be the same size.');
        end
    end
end
varargout = cell(1, nargin);
for k = 1:nargin
    if isscalar(varargin{k})
        varargout{k} = repmat(varargin{k}, sz);
    else
        varargout{k} = varargin{k};
    end
end
end

function s = error_text(code)
names = ["d_km < 0", "h_1 out of range", "h_2 out of range", "h_1 > h_2", ...
    "f_mhz < 100", "f_mhz > 30000", "p < 1", "p > 99", "invalid polarization", ...
    "h_1 == h_2 at zero distance"];
if code >= 1 && code <= numel(names)
    s = names(code);
else
    s = "unknown";
end
end
