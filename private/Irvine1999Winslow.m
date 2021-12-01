function out = Irvine1999Winslow(t, x)
global Vclamp steadyNinfo

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

OpenState = [6,7];

	F=96.5;     			% Faraday's constant (C./mmol)
	T=310;  				% absolute temperature (K)
	R=8.315;				% ideal gas constant (J./[mol.*K])
    k = 1.381e-23;          % Boltzmann's (J./K)
    h = 6.626e-34;          % Planck's (J./s)

if 2 == steadyNinfo
% SS: 4(+open=5) states of 13
% p = [224.114; 708.146; 529.952; 229.205; 39.295; 1.510; -578.317; -130.639; 70.078; 225.175; 338.915; 193.265; 786.217; 0.00711;
%     0; -0.9701; 1.5703; -1.3266; 0.6625; 0; 0; -3.5596; 0; 0; 1.5717; -1.3281; 0; 0;
%     1.4]; % 29 components
% p(1:14) = p(1:14)/R;
% p(15:28) = p(15:28)*F/R/T;
p0 = [224.114; 708.146; 529.952; 229.205; 39.295; 1.510; -578.317; -130.639; 70.078; 225.175; 338.915; 193.265; 786.217; 0.00711;
    -0.9701; 1.5703; -1.3266; 0.6625; -3.5596; 1.5717; -1.3281;
    1.4]; % 29 components
p0(1:14) = p0(1:14)/R;
p0(15:21) = p0(15:21)*F/R/T;
    out = {'Irvine1999Winslow',OpenState,p0};
    return;
end
p = x(end-21:end);
x = x(1:end-22);

% 1C0, 2C1, 3C2, 4C3, 5C4, 6O1, 7O2, 8I0, 9I1, 10I2, 11I3, 12I4, 13II
% 14 rate constants
dHoRT = [116900, 263870, 200240, 127970, 62385, 79035, -99967, 62555, 79183, 123020, 150333, 121900, 293270, 57533]'/R/T;
dSoR = p(1:14);
z = p(15:21);
a = p(end);
logkTh = log(k*T/h);
% lam = k.*T./h .* exp(-dH./R./T + dS./R + z.*F.*V./R./T);
%     lam(i) = (k.*T./h) .* exp( (dS(i) + (z(i).*F.*V-dH(i))./T)./R);
alpha = exp(logkTh+ dSoR(1) -dHoRT(1));
beta = exp(logkTh+ dSoR(2) + z(1).*V -dHoRT(2));
gamma = exp(logkTh+ dSoR(3) + z(2).*V -dHoRT(3));
delta = exp(logkTh+ dSoR(4) + z(3).*V -dHoRT(4));
o_n = exp(logkTh+ dSoR(5) + z(4).*V -dHoRT(5));
o_f = exp(logkTh+ dSoR(6) -dHoRT(6));
gamgam = exp(logkTh+ dSoR(7) -dHoRT(7));
deldel = exp(logkTh+ dSoR(8) + z(5).*V -dHoRT(8));
epsilon = exp(logkTh+ dSoR(9) -dHoRT(9));
omega = exp(logkTh+ dSoR(10) -dHoRT(10));
eta = exp(logkTh+ dSoR(11) + z(6).*V -dHoRT(11));
nu = exp(logkTh+ dSoR(12) + z(7).*V -dHoRT(12));
c_n = exp(logkTh+ dSoR(13) -dHoRT(13));
c_f = exp(logkTh+ dSoR(14) -dHoRT(14));
% dS = [224.114, 708.146, 529.952, 229.205, 39.295, 1.510, -578.317, -130.639, 70.078, 225.175, 338.915, 193.265, 786.217, 0.00711]';
% z = [0, -0.9701, 1.5703, -1.3266, 0.6625, 0, 0, -3.5596, 0, 0, 1.5717, -1.3281, 0, 0]';
% a = 1.4;
% lambda = k.*T./h .* exp(-dH./R./T + dS./R + z.*F.*V./R./T); temp = mat2cell(lambda, ones(length(lambda), 1));
% [alpha, beta, gamma, delta, o_n, o_f, gamgam, deldel, epsilon, omega, eta, nu, c_n, c_f] = temp{:};

act = [4.*alpha; 3.*alpha; 2.*alpha; alpha; gamma; epsilon; 0; 4.*alpha.*a; 3.*alpha.*a; 2.*alpha.*a; alpha.*a; gamgam];
deact = [beta; 2.*beta; 3.*beta; 4.*beta; delta; omega; 0; beta./a; 2.*beta./a; 3.*beta./a; 4.*beta./a; deldel];
inact = [c_n; c_n.*a; c_n.*a^2; c_n.*a^3; c_n.*a^4; o_n];
react = [c_f; c_f./a; c_f./a^2; c_f./a^3; c_f./a^4; o_f];

    n = 13;
    A(1:n+1:n^2) = [-act; 0] + [0; -deact] + [-inact; zeros(n-6,1)] + [zeros(n-6,1); -react];
    A(2:n+1:n^2) = act;
    A(n+1:n+1:n^2) = deact;
    A(8:n+1:7*n) = inact;
    A(7*n+1:n+1:n^2) = react;
    A = setA(A,n,5,7,eta,nu);
% A = diag(act, -1) + diag(deact, +1) + diag(inact, -7) + diag(react, +7);
% A(7,5) = eta; A(5,7) = nu;
% A = A - diag(sum(A));

if 1 == steadyNinfo % steady states OO-> C0 -> CI0 -> II
    alphaM = [deact(6:-1:1); inact(1); act(8:12)];
    betaM = [act(6:-1:1); react(1); deact(8:12)];
    gammaM = alphaM./betaM;
    s = [1; cumprod(gammaM)];
    s = normToOne(s);
    if isa(s, 'myAD')
        condA = cond(reshape(getvalue(A),n,n));
    else
        condA = cond(reshape(A,n,n));
    end
    out = {s([7:-1:1 8:12]), condA};
elseif 3 == steadyNinfo % get rate constants and matrix A
    trans = [alpha; beta; gamma; delta; o_n; o_f; gamgam; deldel; epsilon; omega; eta; nu; c_n; c_f];
    transStr = {'alpha'; 'beta'; 'gamma'; 'delta'; 'o_n'; 'o_f'; 'gamgam'; 'deldel'; 'epsilon'; 'omega'; 'eta'; 'nu'; 'c_n'; 'c_f'};
    A = reshape(A,n,n);
    out = {A, trans, transStr};
else
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
