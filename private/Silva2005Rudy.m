function out = Silva2005Rudy(t, x)
global Vclamp steadyNinfo

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

OpenState = [10;11];

if 2 == steadyNinfo
% SS: 3-5? of 17
%%% NO MICROSCOPIC REVERSABILITY
    p = [3.98e-4; 3.61e-1; 5.74e-5; 9.23e-2; 3.41e-3; 8.68e-1; 1.10e-3; 3.3e-1;
         6.47e-3; 1.25e-2; 4.81e-1; 6.33e-3; 1.27; 4.91e-3; 6.79e-1];

    out = {'Silva2005Rudy',OpenState,p};
    return;
end
p = x(end-14:end);
x = x(1:end-15);

F = 96485;
R = 8314;
T = 310;
VFonRT = V*F/(R*T);

alpha = p(1)*exp(p(2)*VFonRT);
beta = p(3)*exp(-p(4)*VFonRT);
gamma = p(5)*exp(p(6)*VFonRT);
delta = p(7)*exp(-p(8)*VFonRT);
phi = p(9);
eta = p(10)*exp(-p(11)*VFonRT);
psi = p(12)*exp(p(13)*VFonRT);
omega = p(14)*exp(-p(15)*VFonRT);
% alpha = 3.98e-4*exp(3.61e-1*VFonRT);
% beta = 5.74e-5*exp(-9.23e-2*VFonRT);
% gamma = 3.41e-3*exp(8.68e-1*VFonRT);
% delta = 1.10e-3*exp(-3.3e-1*VFonRT);
% phi = 6.47e-3;
% eta = 1.25e-2*exp(-4.81e-1*VFonRT);
% psi = 6.33e-3*exp(1.27*VFonRT);
% omega = 4.91e-3*exp(-6.79e-1*VFonRT);

% C1 C2 C3 C4 C5 C9 C12 C14 C15 O1 O2 || C6 C7 C8 C11 C13 || C10
alphaM = [4*alpha; 3*alpha; 2*alpha; alpha; 4*gamma; 3*gamma; 2*gamma; gamma; phi; psi];
betaM = [beta; 2*beta; 3*beta; 4*beta; delta; 2*delta; 3*delta; 4*delta; eta; omega];

alphaS = [3*alpha; 2*alpha; 2*gamma; gamma];
betaS = [2*beta; 3*beta; 2*delta; 3*delta];

    n = 17;
    A(1:n+1:n^2) = -[alphaM; 0; alphaS; 0; 0] - [0; betaM; 0; betaS; 0];
    A(2:n+1:n^2) = [alphaM; 0; alphaS; 0];
    A(n+1:n+1:n^2) = [betaM; 0; betaS; 0];
    for i = 1:3
        A = setA(A,n,1+i,11+i,gamma,delta);
        A = setA(A,n,5+i,13+i,4*beta,alpha);
    end
    A = setA(A,n,13,17,gamma,2*delta);
    A = setA(A,n,15,17,3*beta,2*alpha);
% A = diag([alphaM; 0; alphaS; 0],-1) + diag([betaM; 0; betaS; 0],+1);
% for i=1:3
%     A(1+i,11+i) = delta; % C6->C2
%     A(11+i,1+i) = gamma;
%     A(5+i,13+i) = alpha; % C8->C9
%     A(13+i,5+i) = 4*beta;
% end
% A(13,17) = 2*delta; % C10->C7
% A(17,13) = gamma;
% A(15,17) = 2*alpha; % C10->C11
% A(17,15) = 3*beta;
% 
% A = A - diag(sum(A));
   
if 1 == steadyNinfo % steady states
    if isa(A, 'myAD')
        out = {A([]), []};
    else
        A = reshape(A,n,n);
        B = A(2:end,2:end);
        b = A(2:end,1);
        B = B - b(:,ones(length(A)-1,1));
        steadyStatesB = B\(-b);
        s = [1-sum(steadyStatesB); steadyStatesB];
        out = {s(1:end-1), cond(reshape(A,n,n))};
    end
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
