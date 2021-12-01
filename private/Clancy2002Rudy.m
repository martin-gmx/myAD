function out = Clancy2002Rudy(t, x)
global Vclamp steadyNinfo

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

OpenState = 1;

if 2 == steadyNinfo
% SS: 4(+open=5) of 9 states
% 3 params not identifiable per se
p = [3.802; 0.1027; 150; 17; 0.20; 15; 0.23; 12; 0.25;
     20.3; 0.1917; 0.20; 5; 0.22; 10; 3.7933e-7; 7.7;
     0.0084; 0.00002; 9.178; 29.68; 100; 9.5e4; 50];

    out = {'Clancy2002Rudy',OpenState,p};
    return;
end
p = x(end-23:end);
x = x(1:end-24);

% states: 1.O 2.C1 3.C2 4.C3 5.IC3 6.IC2 7.IF 8.IM1 (9.IM2)
a1 = p(1)./(p(2).*exp(-V./p(4))+p(5).*exp(-V./p(3)));
a2 = p(1)./(p(2).*exp(-V./p(6))+p(7).*exp(-V./p(3)));
a3 = p(1)./(p(2).*exp(-V./p(8))+p(9).*exp(-V./p(3)));
b1 = p(11).*exp(-V./p(10));
b2 = p(12).*exp(-(V-p(13))./p(10));
b3 = p(14).*exp(-(V-p(15))./p(10));
c1 = p(16).*exp(-V./p(17));
d1 = p(18)+p(19).*V;
f1 = p(20).*exp(V./p(21));
% a1 = 3.802./(0.1027.*exp(-V./17.0)+0.20.*exp(-V./150));
% a2 = 3.802./(0.1027.*exp(-V./15.0)+0.23.*exp(-V./150));
% a3 = 3.802./(0.1027.*exp(-V./12.0)+0.25.*exp(-V./150));
% b1 = 0.1917.*exp(-V./20.3);
% b2 = 0.20.*exp(-(V-5)./20.3);
% b3 = 0.22.*exp(-(V-10)./20.3);
% c1 = 3.7933e-7.*exp(-V./7.7);
% d1 = 0.0084+0.00002.*V;
% f1 = 9.178.*exp(V./29.68);

alphaM = [b3; b2; b1; d1; a1; a2; f1./p(22); f1./p(23)];
betaM = [a3; a2; a1; c1; b1; b2; c1; c1./p(24)];

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
