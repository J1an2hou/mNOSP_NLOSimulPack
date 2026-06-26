function deltafunct = delta_funct(x,eps0)
deltafunct = eps0./(x.^2+eps0^2)/pi;
%deltafunct = 1/eps0/sqrt(pi)*exp(-x.^2/eps0^2);
end