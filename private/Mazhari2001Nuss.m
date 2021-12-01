function out = Mazhari2001Nuss(t, x)
global Vclamp steadyNinfo

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

OpenState = 4;

if 2 == steadyNinfo
% SS: 3 states of 5
p = [0.017147641733086d0; 0.03304608038835d0;
	0.03969328381141d0; -0.04306054163980d0;
	0.02057448605977d0; 0.02617412715118d0;
	0.00134366604423d0; -0.02691385498399d0;
	0.10666316491288d0; 0.00568908859717d0;
	0.00646393910049d0; -0.04536642959543d0;
	0.00008039374403d0; 0.00000069808924d0;
    0.02608362043337d0; 0.14832978132145d0];
    
    out = {'Mazhari2001Nuss',OpenState,p};
    return;
end
p = x(end-15:end);
x = x(1:end-16);

% HERG+hKCNE2 at T=296
C1toC2 = p(1).*exp(p(2).*V);
C2toC1 = p(3).*exp(p(4).*V);
C2toC3 = p(15);
C3toC2 = p(16);
C3toO = p(5).*exp(p(6).*V);
OtoC3  = p(7).*exp(p(8).*V);
OtoI   = p(9).*exp(p(10).*V);
ItoO   = p(11).*exp(p(12).*V);
C3toI  = p(13).*exp(p(14).*V);
ItoC3  = (OtoC3.*ItoO.*C3toI)./(C3toO.*OtoI);

% C1, C2, C3, O, I - with C3toI
alpha = [C1toC2; C2toC3; C3toO; OtoI];
beta = [C2toC1; C3toC2; OtoC3; ItoO];
    
    n = 5;
    A(1:n+1:n^2) = [-alpha; 0] + [0; -beta];
    A(2:n+1:n^2) = alpha;
    A(n+1:n+1:n^2) = beta;
    A = setA(A,n,3,5,C3toI,ItoC3);
% A = diag([-alpha; 0])+diag([0; -beta])+diag(alpha,-1)+diag(beta,+1);
% xs = 3; ys = 5;
% A(xs,xs) = A(xs,xs) -C3toI;
% A(ys,ys) = A(ys,ys) -ItoC3;
% A(ys,xs) = A(ys,xs) +C3toI;
% A(xs,ys) = A(xs,ys) +ItoC3;
    
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
