function [M_d, A_d0, d_crx_km, CASE, warnings] = p528_transhorizon_search_rev1(path, t1, t2, f_mhz, A_dML_db, M_d, A_d0, warnings)
%P528_TRANSHORIZON_SEARCH_REV1 Diffraction/troposcatter crossover search, ITU-R P.528-5 Annex 2 sec. 3, Step 6.
%   [M_d, A_d0, d_crx_km, CASE, warnings] = P528_TRANSHORIZON_SEARCH_REV1(path,
%   t1, t2, f_mhz, A_dML_db, M_d, A_d0, warnings) steps out from the
%   line-of-sight limit in 1 km increments for up to 100 km.
%
%   Port of src/p528/TranshorizonSearch.cpp (github.com/NTIA/p528). As in
%   the C++, the no-crossover default returns CONST_MODE_DIFFRACTION, whose
%   value equals CASE_1.
c = p528_constants_rev1();
CASE = c.CONST_MODE_SEARCH; %#ok<NASGU> mirrors the C++ initialization
k = 0;
d_search = [path.d_ML_km + 3, path.d_ML_km + 2];   % d', d" (Eqns 3-8, 3-9)
A_s = [0, 0];
SEARCH_LIMIT = 100;
for i_search = 1:SEARCH_LIMIT
    A_s(2) = A_s(1);
    tr = p528_troposcatter_rev1(t1, t2, d_search(1), f_mhz);   % Step 6.2
    A_s(1) = tr.A_s_db;
    if tr.A_s_db < 20.0          % not within the valid part of the model
        d_search(2) = d_search(1);
        d_search(1) = d_search(1) + 1;
        continue
    end
    k = k + 1;
    if k <= 1                    % need two points for a line
        d_search(2) = d_search(1);
        d_search(1) = d_search(1) + 1;
        continue
    end
    M_s = (A_s(1) - A_s(2)) / (d_search(1) - d_search(2));    % Eqn 3-10
    if M_s <= M_d
        d_crx_km = d_search(1);
        A_d_db = M_d * d_search(2) + A_d0;                     % Eqn 3-11
        if A_s(2) >= A_d_db
            CASE = c.CASE_1;
        else
            M_d = (A_s(2) - A_dML_db) / (d_search(2) - path.d_ML_km);   % Eqn 3-12
            A_d0 = A_s(2) - (M_d * d_search(2));                         % Eqn 3-13
            CASE = c.CASE_2;
        end
        return
    end
    d_search(2) = d_search(1);
    d_search(1) = d_search(1) + 1;
end
CASE = c.CONST_MODE_DIFFRACTION;
d_crx_km = d_search(2);
warnings = bitor(warnings, c.WARNING_DFRAC_TROPO_REGION);
end
