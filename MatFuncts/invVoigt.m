function [a,b] = invVoigt(v)
if v == 1
    a = 1;b = 1;
elseif v == 2
    a = 2;b = 2;
elseif v == 3
    a = 3;b = 3;
elseif v == 4
    a = 2;b = 3;
elseif v == 5
    a = 1;b = 3;
elseif v == 6
    a = 1;b = 2;
else
    error('Wrong number in inverse Voigt')
end
end