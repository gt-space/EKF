function dPdt = RiccatiHybrid(t, P, A, Qtil)
    % Integrates the differential Riccati equation for the minimal error covariance.
    P = reshape(P,size(Qtil));
    dPdt = A * P + P * A' + Qtil;
    dPdt = dPdt(:);
end