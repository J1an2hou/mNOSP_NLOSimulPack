function [a,b] = inv_Voigt1D(c)
if c == 1
    a = 2; b = 3;
elseif c == 2
    a = 3; b = 1;
elseif c == 3
    a = 1; b = 2;
else
    warning('Wrong direction index')
end

end