% converts radius vector in the ECEF frame to geodetic coordinates
% sources:
% Jekeli, Christopher (August, 2016). Geometric Reference Systems in Geodesy. Ohio State University.
% data from TOPEX/Poseidon (1992) following IERS recommendations
% algorithm from Fukushima, Toshio (December, 1999). Fast transform from geocentric to geodetic coordinates. Journal of Geodesy. National Astronomical Observatory.
function geo_lla = improved_ecef2lla(r)
    a = 6378136.3; % semi-major axis in m
    f = (298.257)^-1; % flattening
    b = a * (1.0 - f); % semi-minor axis
    e = sqrt(1.0 - (b^2 / a^2)); % first eccentricity
    epsilon = 10^-12; % allowable error in the Newton-Raphson solver
    max_iter = 100; % maximum allowable iterations before truncation
    x = r(1);
    y = r(2);
    p = sqrt(x^2 + y^2);
    z = r(3);
    z_norm = abs(z);
    if p == 0
        if z == 0
            phi = 0.0;
            lambda = 0.0;
            h_geo = -b;
        else
            if z > 0
                phi = pi / 2.0;
            else
                phi = -pi / 2.0;
            end
            lambda = 0.0;
            h_geo = z_norm - b;
        end
    elseif z == 0
        phi = 0.0;
        lambda = atan2(y, x);
        h_geo = p - a;
    else
        % in the non-trivial case, we are essentially solving Bowring's formula, but parameterized into a more efficient quartic equation
        % this saves computation time by avoiding unnecessary expensive calls to transcendental functions
        % to avoid the issue of the Newton-Raphson solver converging on the incorrect root, we must first initialize the iteration variable, t, for a number of different circumstances: 
        c = a * e^2;
        e_prime = sqrt(1 - e^2);
        z_prime = e_prime * z_norm;
        u = 2.0 * (z_prime - c);
        v = 2.0 * (z_prime + c);
        t_M = (c - z_prime) / p;
        f_M = p * t_M^4 + u * t_M^3 + v * t_M - p;
        t_1 = (p - c + z_prime) / (p - c + 2.0 * z_prime);
        t_0 = p / (z_prime + c);
        if t_M <= 0
            t = t_1;
        elseif t_M >= 1 
            t = t_0;
        elseif f_M >= 0
            t = t_0;
        else
            t = t_1;
        end
        for i = 1:max_iter
            delta_t = (p - (p * t^4 + u * t^3 + v * t)) / (4.0 * p * t^3 + 3.0 * u * t^2 + v);
            prev_t = t;
            t = t + delta_t;
            TE = abs(t - prev_t); % truncation error
            if TE <= epsilon
                break
            end
        end
        if r(3) > 0
            phi = atan2(1.0 - t^2, 2.0 * e_prime * t);
        else
            phi = -atan2(1.0 - t^2, 2.0 * e_prime * t);
        end
        h_geo = ((2.0 * p * e_prime * t) + (z_norm * (1.0 - t^2)) - (a * e_prime * (1.0 + t^2))) / (sqrt((1.0 + t^2)^2 - (4.0 * e^2 * t^2)));
        lambda = atan2(y, x);
        if lambda < 0
            lambda = 2.0 * pi + lambda;
        end
    end
    geo_lla = [phi; lambda; h_geo];
end