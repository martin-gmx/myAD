function out = LiuRasmusson_Ito_model_MM(t,x)
global Vclamp steadyNinfo
% No-block

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

OpenState = 4;

if 2 == steadyNinfo
    p = [0.03577; 3.8102; 0.06237; 4.5941;
         28.3; 0.7075];
    out = {'LR_MM',OpenState,p};
    return;
end
p = x(end-5:end);
x = x(1:end-6);

alpha_a = exp(p(1)*V+p(2))*1e-3;
beta_a = exp(-p(3)*V+p(4))*1e-3;
Kf = p(5)*1e-3; Kb = p(6)*1e-3;

pp = 360; % 1; % 360;
% Block dynamics
% Kon = 4.5; Koff = 3.4;

%% Markov-model C1-C2-C3-O-I1-I2-I3
alphaM = [3*alpha_a; 2*alpha_a; alpha_a; Kf; 3*beta_a; 2*beta_a./pp];
betaM = [beta_a; 2*beta_a; 3*beta_a; Kb; alpha_a; 2*alpha_a];

    n = 7;
    A(1:n+1:n^2) = [-alphaM; 0] + [0; -betaM];
    A(2:n+1:n^2) = alphaM;
    A(n+1:n+1:n^2) = betaM;
% C2 - I3 (2 - 7)
    A = setA(A,n,2,7,Kf,Kb*pp);
% C3 - I2 (3 - 6)
    A = setA(A,n,3,6,Kf,Kb);

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
