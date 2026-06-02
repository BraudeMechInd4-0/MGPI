function [coeffs_all, powers_all] = gegenbauer_precompute(N, alpha)
    % Precompute coefficients and powers for gegenbauer_eval
    % coeffs_all: cell array, coeffs_all{n+1} is 1 x (floor(n/2)+1)
    % powers_all: cell array, powers_all{n+1} is 1 x (floor(n/2)+1)
    % Uses gammaln for numerical stability (avoids overflow for large arguments)
    % Guard: alpha=0 is the Chebyshev limit and must be routed to CGL branch before calling this
    if alpha == 0
        error('gegenbauer_precompute: alpha=0 is the Chebyshev limit, use the CGL branch instead');
    end
    coeffs_all = cell(N+1, 1);
    powers_all = cell(N+1, 1);
    for n = 0:N
        k = 0:floor(n/2);
        coeffs_all{n+1} = (-1).^k .* exp(gammaln(n-k+alpha) - gammaln(alpha) - gammaln(k+1) - gammaln(n-2*k+1));
        powers_all{n+1} = n - 2*k;
    end
end