function out = Wang1997Rasmusson(t, x)
global Vclamp steadyNinfo

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

T = 310;
Ko = 2.0;

OpenState = 4;

if 2 == steadyNinfo
% SS: 3 of 5
p = [0.022348; 0.01176; 0.047002; 0.0631; 0.023761; 0.036778;
     0.013733; 0.038198; 6.89e-5; 0.04178;
     0.090821; 0.023391; 0.006497; 0.03268];

    out = {'Wang1997Rasmusson',OpenState,p};
    return;
end
p = x(end-13:end);
x = x(1:end-14);

C0toC1  = p(1).*exp( p(2).*V.*T./296);
C1toC0  = p(3).*exp(-p(4).*V.*T./296);
C1toC2  = p(5);
C2toC1  = p(6);
C2toO   = p(7).*exp( p(8).*V.*T./296);
OtoC2   = p(9).*exp(-p(10).*V.*T./296);
% Ko = 2 mM
OtoI    = p(11).*exp( p(12).*V.*T./296);
ItoO    = p(13).*exp(-p(14).*V.*T./296);
%     C0toC1  = 0.022348 .* exp( 0.01176.*V.*T./296);
%     C1toC0  = 0.047002 .* exp(-0.0631 .*V.*T./296);
%     C1toC2  = 0.023761;
%     C2toC1  = 0.036778;
%     C2toO   = 0.013733 .* exp( 0.038198.*V.*T./296);
%     OtoC2   = 6.89e-5.* exp(-0.04178 .*V.*T./296);
%     % Ko = 2 mM
%     OtoI    = 0.090821 .* exp( 0.023391.*V.*T./296);
%     ItoO    = 0.006497 .* exp(-0.03268.*V.*T./296);

    ItoO = ItoO .* (2.0./Ko)^0.4;

    % C0, C1, C2, O, I
    alpha = [C0toC1; C1toC2; C2toO; OtoI];
    beta = [C1toC0; C2toC1; OtoC2; ItoO];

    n = 5;
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
