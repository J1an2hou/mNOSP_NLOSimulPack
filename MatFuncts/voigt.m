function [a_out, b_out] = voigt(a_in,b_in)
% To convert two indice from 1-3 into a single Voigt index  from 1-6
if nargin == 2 && nargout == 1
    if a_in == b_in
        a_out = a_in;
    else
        a_out = 9-a_in-b_in;
    end

elseif nargin == 1 && nargout == 2
if a_in == 1
    a_out = 1;b_out = 1;
elseif a_in == 2
    a_out = 2;b_out = 2;
elseif a_in == 3
    a_out = 3;b_out = 3;
elseif a_in == 4
    a_out = 2;b_out = 3;
elseif a_in == 5
    a_out = 1;b_out = 3;
elseif a_in == 6
    a_out = 1;b_out = 2;
else
    error('Wrong number in inverse Voigt')
end

else
    error('Wrong use of input and output index numbers')
end
% End of function
end