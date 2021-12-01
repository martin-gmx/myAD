function out = myClancy2002Rudy(t, x)
global Vclamp steadyNinfo

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

OpenState = 1;

if 2 == steadyNinfo
% SS: 4(+open=5) of 9 states
% 3 params not identifiable per se
p = [0.1027/3.802; 12; 0.25/3.802; 150;
     0.1917; 15.3; 3.7933e-7; 7.7;
     0.0084; 9.178; 29.68; 100; 9.5e4; 50];

    out = {'myClancy2002Rudy',OpenState,p};
    return;
end
p = x(end-13:end);
x = x(1:end-14);

% states: 1.O 2.C1 3.C2 4.C3 5.IC3 6.IC2 7.IF 8.IM1 (9.IM2)
a = 1./(p(1).*exp(-V./p(2))+p(3).*exp(-V./p(4)));
b = p(5)./2.5.*exp(-V./p(6));
c1 = p(7).*exp(-V./p(8));
d1 = p(9);
f1 = p(10).*exp(V./p(11));

alphaM = [3*b; 2*b; b; d1; 3*a; 2*a; f1./p(12); f1./p(13)];
betaM = [a; 2*a; 3*a; c1; b; 2*b; c1; c1./p(14)];

    n = 9;
    A(1:n+1:n^2) = [-alphaM; 0] + [0; -betaM];
    A(2:n+1:n^2) = alphaM;
    A(n+1:n+1:n^2) = betaM;
    A = setA(A,n,1,7,f1,f1.*betaM(1).*c1./(d1.*alphaM(1)));
    A = setA(A,n,2,7,d1,c1);
    A = setA(A,n,3,6,d1,c1);
% A = diag(alphaM,-1) + diag(betaM,+1); % -diag([alphaM; 0]+[0; betaM])
% A(2,7) = c1; A(3,6) = c1;
% A(7,2) = d1; A(6,3) = d1;
% A(7,1) = f1;
% A(1,7) = A(7,1).*A(1,2).*A(2,7)./(A(7,2).*A(2,1));
% A = A - diag(sum(A));

if 1 == steadyNinfo % steady states
    gammaM = alphaM./betaM;
    s = [1; cumprod(gammaM)];
    s = normToOne(s);
    if isa(s, 'myAD')
        condA = cond(reshape(getvalue(A),n,n));
    else
        condA = cond(reshape(A,n,n));
    end
    out = {s(1:end-1), condA};
elseif 3 == steadyNinfo % get rate constants and matrix A
    trans = [alphaM; betaM];
    transStr = {'CO'; 'OI'; 'OC'; 'IO'};
    A = reshape(A,n,n);
    out = {A, trans, transStr};
else % Remove last state
    y = [x; 1-sum(x)];
    if isa(y, 'myAD')
        dy = zeros(y);
    else
        dy = y.*0;
        A = A(:);
    end
    for i=1:n
        dy = dy + A((1:n)+(i-1)*n).*y(i*ones(n,1));
    end
%     dy = A*y;
    out = dy(1:end-1);
end
