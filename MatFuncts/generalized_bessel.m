function J_n = generalized_bessel(n, x, y)
    % Fourier component of x*sin⁡(theta)+y*cos(theta)
    % Its relation with conventional Bessel function of the first kind j_n(x), Fourier of
    % x*sin(theta) is
    % J_n(a, 0) = j_n(a)
    % J_n(0, a) = i^n * j_n(a)
    % J_n(r*cosphi, r*sinphi) = e^{-i*n*phi) * j_n(r)
    if ~isequal(size(x), size(y))
        error('The x and y have to be the same size');
    end
    [phi,rho] = cart2pol(x,y);
    J_n = exp(-1i*n*phi) .* besselj(n, rho);
    
    % Old version for numerical integration
    % sz = size(x);
    % x = x(:);
    % y = y(:);
    % J_n = zeros(size(x), 'like', complex(1));
    % 
    % for k = 1:numel(x)
    %     integrand = @(theta) exp(1i * (n*theta - x(k)*sin(theta) - y(k)*cos(theta)));
    %     J_n(k) = (1/(2*pi)) * integral(integrand, -pi, pi, 'ArrayValued', true);
    % end
    % J_n = reshape(J_n, sz);
    
end