function x = mod(x,y)
% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at

if isa(x, 'myAD')
    if isa(y, 'myAD')
        n = floor(x.values./y.values);
    else
        n = floor(x.values./y);
    end
else
    n = floor(x./y.values);
end
x = x - n.*y;
