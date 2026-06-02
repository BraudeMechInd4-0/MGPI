function C = gegenbauer_eval(coeffs_all, powers_all, tau)
    % Evaluate Gegenbauer matrix at tau points given precomputed coeffs
    % C is (N+1) x length(tau)
    % coeffs_all, powers_all: output of gegenbauer_precompute
    tau = tau(:)';
    N = length(coeffs_all) - 1;
    C = zeros(N+1, length(tau));
    for n = 0:N
        C(n+1,:) = coeffs_all{n+1} * (2*tau).^(powers_all{n+1}');
    end
    C = C';
end