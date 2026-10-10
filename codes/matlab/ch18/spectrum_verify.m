function report = spectrum_verify(lam1, lam2, tol)
%% spectrum_verify - Drift-with-N diagnostic for spectral eigenvalue problems.
%
% Chapter 18: Linear Spectral Eigenproblems --- shared utility.
%
% Implements the drift-with-N diagnostic of Boyd (2000, Ch. 7.5). Given two
% real sorted spectra lam1 (coarse) and lam2 (fine), classifies each coarse
% mode as trusted or suspect by measuring its movement between resolutions
% scaled by the local intermodal separation.
%
% Quantities:
%   sigma_j         (Boyd 7.19) intermodal separation of lam1
%   delta_ord_j     (Boyd 7.20) scaled ordinal drift |lam1(j) - lam2(j)| / sigma_j
%   delta_nst_j     (Boyd 7.21) scaled nearest drift min_k |lam1(j) - lam2(k)| / sigma_j
%
% Input:
%   lam1 - real column/row vector, coarse-resolution spectrum
%   lam2 - real column/row vector, fine-resolution spectrum
%   tol  - trust tolerance (optional; default 1e-3)
%
% Output (struct):
%   report.lam1, report.lam2         - sorted input spectra
%   report.sigma                     - intermodal separation
%   report.delta_ordinal             - ordinal drift
%   report.delta_nearest             - nearest drift
%   report.trusted                   - logical vector: mode j trusted iff max(drifts) < tol
%   report.tol                       - tolerance used
%   report.n_trusted                 - sum(trusted)
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

    if nargin < 3
        tol = 1e-3;
    end

    lam1 = sort(real(lam1(:)));
    lam2 = sort(real(lam2(:)));

    sigma = intermodal_separation(lam1);
    d_ord = ordinal_drift(lam1, lam2, sigma);
    d_nst = nearest_drift(lam1, lam2, sigma);

    trusted = (d_ord < tol) & (d_nst < tol) & isfinite(d_ord) & isfinite(d_nst);

    report = struct();
    report.lam1          = lam1;
    report.lam2          = lam2;
    report.sigma         = sigma;
    report.delta_ordinal = d_ord;
    report.delta_nearest = d_nst;
    report.trusted       = trusted;
    report.tol           = tol;
    report.n_trusted     = sum(trusted);
end


function sigma = intermodal_separation(lam)
%% Intermodal separation sigma_j per Boyd Eq. 7.19.
    n = numel(lam);
    sigma = zeros(n, 1);
    if n == 1
        sigma(1) = abs(lam(1));
        return;
    end
    sigma(1) = abs(lam(1) - lam(2));
    if n > 2
        diffs = abs(diff(lam));          % length n-1
        sigma(2:end-1) = 0.5 * (diffs(1:end-1) + diffs(2:end));
    end
    sigma(end) = abs(lam(end) - lam(end-1));
    tiny = sigma < 1e-14;
    sigma(tiny) = max(abs(lam(tiny)), 1e-14);
end


function d = ordinal_drift(lam1, lam2, sigma)
%% delta_j,ordinal = |lam1(j) - lam2(j)| / sigma_j   (Boyd Eq. 7.20).
    n1 = numel(lam1);
    n2 = numel(lam2);
    d = inf(n1, 1);
    m = min(n1, n2);
    d(1:m) = abs(lam1(1:m) - lam2(1:m)) ./ sigma(1:m);
end


function d = nearest_drift(lam1, lam2, sigma)
%% delta_j,nearest = min_k |lam1(j) - lam2(k)| / sigma_j   (Boyd Eq. 7.21).
    n1 = numel(lam1);
    d = inf(n1, 1);
    if isempty(lam2)
        return;
    end
    % |lam1(j) - lam2(k)| as an (n1 x n2) matrix
    diffs = abs(lam1(:) - lam2(:)');
    mn = min(diffs, [], 2);
    d = mn ./ sigma;
end
