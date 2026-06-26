function phi_ave = FloPhase_ave(A0,R,wc,alatt,nw,nR,q)
% This function perform time average for the phase factor
% The phi_ave is time average of e^[i(A0 * r) + iqt]
% Here <0m | r | Rn > = R - tau_m + tau_n
% A0 is in size of (1,3)
% wr is in size of (nR,nw,nw,3)
% We now apply Bessel function directly to evaluate the Fourier components
% The old version is totally commented here

% T = 2*pi/omega;
% phi_ave = zeros(nR,nw,nw);
% T_weight = 1/NT;
% t_line = linspace(0, T, NT+1);
% t_line = t_line(1:NT);
% wcdiff1 = wc(:,1)*ones(1,nw)-ones(nw,1)*wc(:,1)';
% wcdiff2 = wc(:,2)*ones(1,nw)-ones(nw,1)*wc(:,2)';
% wcdiff3 = wc(:,3)*ones(1,nw)-ones(nw,1)*wc(:,3)';
% for iR = 1:nR
%     aR = R(iR,:)*alatt;
%     aR1 = aR(1) - wcdiff1;
%     aR2 = aR(2) - wcdiff2;
%     aR3 = aR(3) - wcdiff3;
%     t_ave = zeros(nw,nw);
% for it = 1:NT
%     tau =  t_line(it);
%     A_F = A0 .* real(exp(1i*omega*tau)*eta);
%     AdotR = aR1 * A_F(1) + aR2 * A_F(2) + aR3 * A_F(3);
%     phi = AdotR + q * omega * tau;
%     t_ave = t_ave + exp(1i*phi)*T_weight;
% end
% phi_ave(iR,:,:) = t_ave;
% end

% this is the new version that uses generalized bessel function

phi_ave = zeros(nR,nw,nw);
wcdiff1 = wc(:,1)*ones(1,nw)-ones(nw,1)*wc(:,1)';
wcdiff2 = wc(:,2)*ones(1,nw)-ones(nw,1)*wc(:,2)';
% wcdiff3 = wc(:,3)*ones(1,nw)-ones(nw,1)*wc(:,3)';
for iR = 1:nR
    aR = R(iR,:)*alatt;
    aR1 = aR(1) - wcdiff1;
    aR2 = aR(2) - wcdiff2;
    % aR3 = aR(3) - wcdiff3;

    % --- Elliptical / circular / linear unified treatment ---
    % A0 = [Ax, Ay, 0]   (real amplitudes)
    % A(t) = (Ax * sin wt, Ay * cos wt)
    if ~isreal(A0(1)) && ~isreal(A0(2))
        error('The current version only supports real number in light amplitude')
    end
    if abs(A0(3)) > 1e-6
        warning('We usually assume light propagating along z, not it seems not')
    end
    Ax = A0(1);
    Ay = A0(2);
    % Az = real(A0(3));
    AxR = Ax * aR1;   % sin(wt) coefficient
    AyR = Ay * aR2;   % cos(wt) coefficient
    phi_ave(iR,:,:) = generalized_bessel(q, AxR, AyR);
end

end

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