function out = myFaber2007Rudy2nd(t, x)
global Vclamp CassClamp steadyNinfo

V = interp1(Vclamp(:,1),Vclamp(:,2),t);
Cass = interp1(CassClamp(:,1),CassClamp(:,2),t);

OpenState = 8;

if 2 == steadyNinfo
% SS: 2./4(+open)+HH-CaInact of 14
% omegafs=phis
p = [0.925; 30; 0.39; 40; 0.245; 10;
     0.005; 40; 0.02; 0; 0.03; 0;
     0.035; 0; 0.0011; 0; 4; 1; 0.01];

    out = {'myFaber2007Rudy2nd',OpenState,p};
    return;
end
p = x(end-18:end);
x = x(1:end-19);

alpha = p(1).*exp(V./p(2));
beta = p(3).*exp(-V./p(4));
gammaf = p(5).*exp(V./p(6));
gammas = p(7).*exp(-V./p(8));
phif = p(9);
phis = p(11);
lambdaf = p(13);
lambdas = p(15);
omegaf = 4.*beta.*lambdaf.*gammaf./alpha./phif;
omegas = 4.*beta.*lambdas.*gammas./alpha./phis;
omegasf = lambdas.*phif./lambdaf;
omegafs = phis;

delta = p(17)./(p(18)+1./Cass);
theta = p(19);

% CCCCOIsIv
alphaM = [4.*alpha; 3.*alpha; 2.*alpha; alpha; phis; omegasf];
betaM = [beta; 2.*beta; 3.*beta; 4.*beta; lambdas; omegafs];

% include C3-Is, C3-If, O-If
alphaI = [gammas; gammaf; phif];
betaI = [omegas; omegaf; lambdaf];

    n = 7;
    A(1:n+1:n^2) = - [alphaM; 0] - [0; betaM];
    A(2:n+1:n^2) = alphaM;
    A(n+1:n+1:n^2) = betaM;
    idxI = [4, 6; 4, 7; 5, 7];
    for i = 1:3
        A = setA(A,n,idxI(i,1),idxI(i,2),alphaI(i),betaI(i));
    end
    
    inf_CaIn = theta./(delta+theta);
    tau_CaIn = 1./(delta+theta);

if 1 == steadyNinfo % steady states
    gammaM = alphaM./betaM;
    s = [1; cumprod(gammaM)];
    s = normToOne(s);
    if isa(s, 'myAD')
        condA = cond(reshape(getvalue(A),n,n));
    else
        condA = NaN; %cond(reshape(A,n,n));
    end
    out = {[inf_CaIn; s(1:end-1); inf_CaIn.*s(5)], condA};
elseif 3 == steadyNinfo % get rate constants and matrix A
    trans = [gammas; omegas; gammaf; omegaf; phif; lambdaf; phis; lambdas; omegafs; omegasf];
    transStr = {'gammas'; 'omegas'; 'gammaf'; 'omegaf'; 'phif'; 'lambdaf'; 'phis'; 'lambdas'; 'omegafs'; 'omegasf'};
    A = reshape(A,n,n);
    out = {A, trans, transStr};
else
    CaIn = x(1);
    x = x(2:end-1);
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
    dCaIn = (inf_CaIn-CaIn)./tau_CaIn;
    out = [dCaIn; dy(1:end-1); dCaIn.*y(5) + CaIn.*dy(5)];
end
