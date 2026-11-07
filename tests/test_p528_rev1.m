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
