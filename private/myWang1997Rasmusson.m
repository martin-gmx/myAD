function out = myWang1997Rasmusson(t, x)
global Vclamp steadyNinfo

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

T = 310;
Ko = 2.0;

OpenState = 2;

if 2 == steadyNinfo
% SS: 3 of 5
p = [0.038198*1.13; 0.01176/22; 105; 0.047002/500; 0.0631/1.7;
%      0.013733; 0.038198; 6.89e-5; 0.04178;
     0.090821; 0.023391/1.8; 0.006497; 0.03268];

    out = {'myWang1997Rasmusson',OpenState,p};
    return;
end
p = x(end-8:end);
x = x(1:end-9);

C0toO =  exp(-p(1).*V.*T./296 - p(2).*(V-p(3)).^2.*T./296);
OtoC0 =  p(4).*exp(-p(5).*V.*T./296);

% Ko = 2 mM
OtoI    = p(6).*exp( p(7).*V.*T./296);
ItoO    = p(8).*exp(-p(9).*V.*T./296);

    ItoO = ItoO .* (2.0./Ko)^0.4;

    % C0, C1, C2, O, I
    alpha = [C0toO; OtoI];
    beta = [OtoC0; ItoO];

    n = 3;
    A(1:n+1:n^2) = [-alpha; 0] + [0; -beta];
    A(2:n+1:n^2) = alpha;
    A(n+1:n+1:n^2) = beta;
%     A = diag([-alpha; 0])+diag([0; -beta])+diag(alpha,-1)+diag(beta,+1);
   
if 1 == steadyNinfo % steady states
    gammaM = alpha./beta;
    s = [1; cumprod(gammaM)];
    s = normToOne(s);
    if isa(s, 'myAD')
        condA = cond(reshape(getvalue(A),n,n));
    else
        condA = cond(reshape(A,n,n));
    end
    out = {s(1:end-1), condA};
elseif 3 == steadyNinfo % get rate constants and matrix A
    trans = [alpha; beta];
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
