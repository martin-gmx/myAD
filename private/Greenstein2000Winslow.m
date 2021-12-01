function out = Greenstein2000Winslow(t, x)
global Vclamp steadyNinfo

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

OpenState = 1;

if 2 == steadyNinfo
% SS: 4(+open=5) of 10
p = [0.5437; 0.02898; 0.08019; 0.04684;
    0.04984; 3.7302e-4; 8.1948e-4; 5.374e-8;
    1.8936; 14.225; 158.574; 142.937;
    6.7735; 15.621; 28.753; 524.576];

    out = {'Greenstein2000Winslow',OpenState,p};
    return;
end
p = x(end-15:end);
x = x(1:end-16);

alpha_a = p(1).*exp(p(2).*V);
beta_a = p(3).*exp(-p(4).*V);
alpha_i = p(5).*exp(-p(6).*V);
beta_i = p(7).*exp(p(8).*V);
% alpha_a0 = 0.5437; a_a = 0.02898;
% beta_a0 = 0.08019; b_a = 0.04684;
% alpha_i0 = 0.04984; a_i = 3.7302e-4;
% beta_i0 = 8.1948e-4; b_i = 5.374e-8;
% f1 = 1.8936; f2 = 14.225;
% f3 = 158.574; f4 = 142.937;
% b1 = 6.7735; b2 = 15.621;
% b3 = 28.753; b4 = 524.576;
% alpha_a = alpha_a0.*exp(a_a.*V);
% beta_a = beta_a0.*exp(-b_a.*V);
% alpha_i = alpha_i0.*exp(-a_i.*V);
% beta_i = beta_i0.*exp(b_i.*V);

% O C3 C2 C1 C0 CI0 CI1 CI2 CI3 OI
f = p(9:12);
b = p(13:16);
alpha = [4.*beta_a; 3.*beta_a; 2.*beta_a; beta_a; beta_i; b(1).*4.*alpha_a; b(2).*3.*alpha_a./b(1); b(3).*2.*alpha_a./b(2); b(4).*alpha_a./b(3)];
beta = [alpha_a; 2.*alpha_a; 3.*alpha_a; 4.*alpha_a; alpha_i; beta_a./f(1); f(1).*2.*beta_a./f(2); f(2).*3.*beta_a./f(3); f(3).*4.*beta_a./f(4)];


    n = 10;
    A(1:n+1:n^2) = [-alpha; 0] + [0; -beta];
    A(2:n+1:n^2) = alpha;
    A(n+1:n+1:n^2) = beta;
    for i = 1:4
        A = setA(A,n,i,11-i,beta_i.*f(5-i),alpha_i./b(5-i));
    end
% A = diag(alpha,-1) + diag(beta,+1);
% for i = 1:4
%     xs = i; ys = 11-i;
%     A(ys,xs) = beta_i.*f(5-i);
%     A(xs,ys) = alpha_i./b(5-i);
% end
% A = A - diag(sum(A));
   
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
