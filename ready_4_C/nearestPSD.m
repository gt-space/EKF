function P_corrected = nearestPSD(P)
    P = (P + P')/2;
    [V, D] = eig(P);           % Ensure symmetry and decompose
    if sum(D < 1e-10) > 0 % if eigenvalue is negative or close to negative
        D(D < 1e-10) = 1e-10;               % Clamp negative/near-zero eigenvalues
        P_corrected = V * D * V';              % Reconstruct matrix
    else
        P_corrected = P;
    end
end