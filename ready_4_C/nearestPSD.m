function P_corrected = nearestPSD(P)
    P = (P + P')/2;
    [V, D] = eig(P);           % Ensure symmetry and decompose
    d = diag(D);
    if any(d < 1e-10) % if eigenvalue is negative or close to negative
        fprintf('Rats! Covariance not PSD.\n');
        d(d < 1e-10) = 1e-8;               % Clamp negative/near-zero eigenvalues
        P_corrected = zeros(size(P));
        for i = 1:length(d)
            vi = V(:,i);
            P_corrected = P_corrected + d(i) * (vi*vi');              % Reconstruct matrix
        end
    else
        fprintf('Yay! Covariance is PSD.\n');
        P_corrected = P;
    end
end