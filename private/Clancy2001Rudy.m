function out = Clancy2001Rudy(t, x)
global Vclamp steadyNinfo

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

Ko = 2.0;

OpenState = 4;

if 2 == steadyNinfo
% SS: all 5 states
% p(2) twice, C1toI = C1toO. 3 params not identifiable per se
p = [55.5e-3; 0.05547153; 12; 2.357e-3; 0.036588; 2.172; 1.077;
     65.5e-3; 36; 2.9357e-3; 0.02158; 0.656; 0.000942; 0.439;
     0.02352; 25];

    out = {'Clancy2001Rudy',OpenState,p};
    return;
end
p = x(end-15:end);
x = x(1:end-16);

    C3toC2 = p(1) .* exp(p(2).*(V-p(3)));
    C2toC3 = p(4) .* exp(-p(5).*(V));
    C2toC1 = p(6);
    C1toC2 = p(7);
    C1toO  = p(8) .* exp(p(2).*(V-p(9)));
    OtoC1  = p(10) .* exp(-p(11).*(V));
    OtoI   = p(12) .* exp(p(13).*(V)) .* (5.4./Ko)^0.3;
    ItoO   = p(14) .* exp(-p(15).*(V+p(16))) .* 5.4./Ko;
%     C3toC2 = 55.5e-3 .* exp(0.05547153.*(V-12));
%     C2toC3 = 2.357e-3 .* exp(-0.036588.*(V));
%     C2toC1 = 2.172;
%     C1toC2 = 1.077;
%     C1toO  = 65.5e-3 .* exp(0.05547153.*(V-36));
%     OtoC1  = 2.9357e-3 .* exp(-0.02158.*(V));
%     OtoI   = 0.656 .* exp(0.000942.*(V)) .* (5.4./Ko)^0.3;
%     ItoO   = 0.439 .* exp(-0.02352.*(V+25)) .* 5.4../Ko;
    C1toI = C1toO;
    ItoC1  = (ItoO .* OtoC1 .* C1toO)./(C1toO .* OtoI);

    % C3, C2, C1, O, I - and ItoC1
    alpha = [C3toC2; C2toC1; C1toO; OtoI];
    beta = [C2toC3; C1toC2; OtoC1; ItoO];

    n = 5;
    A(1:n+1:n^2) = [-alpha; 0] + [0; -beta];
    A(2:n+1:n^2) = alpha;
    A(n+1:n+1:n^2) = beta;
    A = setA(A,n,3,5,C1toI,ItoC1);
%     A = diag([-alpha; 0])+diag([0; -beta])+diag(alpha,-1)+diag(beta,+1);
%     xs = 3; ys = 5;
%     A(xs,xs) = A(xs,xs) -C1toI;
%     A(ys,ys) = A(ys,ys) -ItoC1;
%     A(ys,xs) = A(ys,xs) +C1toI;
%     A(xs,ys) = A(xs,ys) +ItoC1;

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
