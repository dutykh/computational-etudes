function [D, x] = cheb(N)
% CHEB  Chebyshev differentiation matrix of size (N+1) x (N+1).
%   [D, x] = cheb(N) returns the differentiation matrix D and
%   the Chebyshev nodes x = cos(pi*(0:N)/N)'.
%
%   Reference: Trefethen, Spectral Methods in MATLAB, SIAM 2000.

    if N == 0, D = 0; x = 1; return, end
    x = cos(pi*(0:N)/N)';
    c = [2; ones(N-1,1); 2] .* (-1).^(0:N)';
    X = repmat(x, 1, N+1);
    dX = X - X' + eye(N+1);
    D = (c * (1./c)') ./ dX;
    D = D - diag(sum(D,2));
end
