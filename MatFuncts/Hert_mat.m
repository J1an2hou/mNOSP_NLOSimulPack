function A = Hert_mat(A)
% Use with great care!
% This function hermitianize a complex matrix that is slightly deviated
% from Hermitian
B = A';
D = (B - A)/2;
if norm(D) > 0.1
    warning('Your input matrix is largely deviated from Hermitian')
end
A = A + D;
A = roundn(A,ceil(log10(norm(D))));
if ~ishermitian(A)
    error('Hermitianization failed')
end

end