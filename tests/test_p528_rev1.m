function tests = test_p528_rev1
%TEST_P528_REV1 Validates the MATLAB port against the NTIA C++ reference implementation.
%   Run from the repository root: results = runtests('tests/test_p528_rev1.m')
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'));
tc.TestData.ref = readtable(fullfile(root, 'tests', 'reference_ntia_cpp.csv'));
end

function test_readme_examples(tc)
% Example values published in the NTIA README (to 0.1 dB).
A = p528_rev1([15 100 1500 30], [10 100 15 8], [1000 15000 10000 20000], ...
    [500 3600 5700 22000], [0 0 0 1], [50 90 10 50]);
verifyEqual(tc, round(A, 1), [110.0 151.6 293.4 151.1]);
end

function test_against_cpp_reference(tc)
% 800 cases from the compiled NTIA C++ code (tools/p528_reference_driver.cpp).
r = tc.TestData.ref;
[A, res, rtn] = p528_rev1(r.d_km, r.h_1_m, r.h_2_m, r.f_mhz, r.T_pol, r.p);
verifyEqual(tc, rtn, r.rtn);
verifyEqual(tc, [res.propagation_mode]', r.propagation_mode);
verifyEqual(tc, [res.warnings]', r.warnings);
verifyEqual(tc, A, r.A_db, 'AbsTol', 1e-6);
verifyEqual(tc, [res.A_fs_db]', r.A_fs_db, 'AbsTol', 1e-6);
verifyEqual(tc, [res.A_a_db]', r.A_a_db, 'AbsTol', 1e-6);
verifyEqual(tc, [res.d_km]', r.res_d_km, 'AbsTol', 1e-6);
end

function test_invalid_inputs(tc)
verifyError(tc, @() p528_rev1(-1, 10, 100, 1000, 0, 50), 'p528:invalidInput');
[A, ~, rtn] = p528_rev1([-1 10 10 10 10 10 10], [10 1 10 200 10 10 10], [100 100 100 100 100 100 100], ...
    [1000 1000 99 1000 1000 31000 1000], [0 0 0 0 2 0 0], [50 50 50 50 50 50 0.5], 'error_mode', "code");
verifyTrue(tc, all(isnan(A)));
verifyEqual(tc, rtn, [1 2 5 4 9 6 7]);
end

function test_zero_distance_equal_heights(tc)
% Reference returns zero loss with SUCCESS.
[A, ~, rtn] = p528_rev1(0, 100, 100, 1000, 0, 50);
verifyEqual(tc, A, 0);
verifyEqual(tc, rtn, 0);
end

function test_free_space_floor(tc)
% Short, steep LOS paths at 50% stay within a few dB of free space.
[A, res] = p528_rev1(0.5, 10, 121.92, 8160, 1, 50);
verifyLessThan(tc, abs(A - res.A_fs_db), 3);
end

function test_parallel_matches_serial(tc)
pool = gcp('nocreate');
tc.assumeNotEmpty(pool, 'No parallel pool open; parallel check filtered.');
d = [1 10 60 120];
A0 = p528_rev1(d, 10, 121.92, 8160, 1, 50);
A1 = p528_rev1(d, 10, 121.92, 8160, 1, 50, 'parallel_flag', 1);
verifyEqual(tc, A1, A0);
end

function test_transhorizon_reference(tc)
% 300 cases from the NTIA C++ code aimed at the diffraction and troposcatter
% regions near the radio horizon, the h > 20 km warnings, equal heights, and
% a fine time-percentage sweep (tests/reference_cases_transhorizon_input.txt).
root = fileparts(fileparts(mfilename('fullpath')));
r = readtable(fullfile(root, 'tests', 'reference_ntia_cpp_transhorizon.csv'));
[A, res, rtn] = p528_rev1(r.d_km, r.h_1_m, r.h_2_m, r.f_mhz, r.T_pol, r.p);
verifyEqual(tc, rtn, r.rtn);
verifyEqual(tc, [res.propagation_mode]', r.propagation_mode);
verifyEqual(tc, [res.warnings]', r.warnings);
verifyEqual(tc, A, r.A_db, 'AbsTol', 1e-6);
verifyEqual(tc, [res.A_a_db]', r.A_a_db, 'AbsTol', 1e-6);
% Every propagation mode and both warning types are covered.
verifyTrue(tc, all(ismember([1 2 3], [res.propagation_mode])));
verifyTrue(tc, any(bitand([res.warnings], 1)) || any([res.warnings] == 4));
end

function test_components_against_cpp(tc)
% P.835 atmosphere, P.676 specific attenuation and slant path vs the C++
% (tools/component_reference_driver.cpp -> tests/reference_components_cpp.txt).
root = fileparts(fileparts(mfilename('fullpath')));
lines = readlines(fullfile(root, 'tests', 'reference_components_cpp.txt'));
lines = lines(strlength(lines) > 0);
for k = 1:numel(lines)
    parts = split(lines(k));
    v = str2double(parts(2:end));
    switch parts(1)
        case "atm"
            [T, p, e] = p835_atmosphere_rev1(v(1));
            verifyEqual(tc, [T p e], v(2:4)', 'RelTol', 1e-12);
        case "gam"
            [T, p, e] = p835_atmosphere_rev1(0);
            g = p676_specific_attenuation_rev1(v(1), T, e, p);
            verifyEqual(tc, g, v(2), 'RelTol', 1e-10);
        case "slt"
            s = p676_slant_path_rev1(v(1), v(2), v(3), v(4));
            verifyEqual(tc, [s.A_gas_db s.bending_rad s.a_km s.angle_rad s.delta_L_km], v(5:9)', 'RelTol', 1e-9);
    end
end
end

function test_physical_properties(tc)
% Loss grows with distance, falls with time percentage, and never drops
% unrealistically far below free space on short LOS paths.
d = [0.5 1 2 5 10 20 40 80 150 300];
A = p528_rev1(d, 10, 1000, 3000, 0, 50);
verifyTrue(tc, all(diff(A) > 0), 'Loss must increase with distance.');
pp = [1 5 10 50 90 99];
Ap = p528_rev1(200, 15, 10000, 1200, 0, pp);
verifyTrue(tc, all(diff(Ap) >= 0), 'Loss must not decrease as p increases.');
[A1, r1] = p528_rev1(1, 10, 121.92, 8160, 1, 1);
verifyGreaterThan(tc, A1, r1.A_fs_db - 10);
end

function test_free_space_formula(tc)
% The reported free-space loss uses 20log10(f) + 20log10(r) + 32.45.
[~, res] = p528_rev1(5, 100, 1000, 2000, 1, 50);
verifyEqual(tc, res.A_fs_db, 20*log10(2000) + 20*log10(sqrt((5)^2 + 0.9^2)) + 32.45, 'AbsTol', 0.2);
end

function test_unreachable_grazing_point_errors(tc)
% The C++ hangs here; the port raises an error.
verifyError(tc, @() p676_slant_path_rev1(8.16, 0.12192, 0.5, 1.6), 'p676:noGrazingPoint');
end

function test_variability_helpers(tc)
% Inverse CCDF: Q(0.1) = 1.2816 (normal table); symmetric about 0.5.
q = p528_variability_rev1("inv_ccdf", 0.1);
verifyEqual(tc, q, 1.2816, 'AbsTol', 5e-4);
verifyEqual(tc, p528_variability_rev1("inv_ccdf", 0.9), -q, 'AbsTol', 1e-12);
% Nakagami-Rice at the 50% point is zero for every K.
for K = [-40 -10 0 20]
    verifyEqual(tc, p528_variability_rev1("nakagami_rice", K, 50), 0, 'AbsTol', 1e-12);
end
% find_k_99 inverts the 99% column at table points.
c = p528_constants_rev1();
verifyEqual(tc, p528_variability_rev1("find_k_99", c.NakagamiRiceCurves(10, 17) + 1e-9), c.K(10), 'AbsTol', 1e-6);
end
