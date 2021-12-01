function out = Priebe1998Beuckelmann(t, x)
global Vclamp steadyNinfo

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

OpenState = 4;

if 2 == steadyNinfo
p = [0.13; 10.66; 11.1; 0.3; 2.525e-7; 0.1; 32;
    0.135; 80; 6.8; 3.56; 0.079; 3.1e5; 0.35;
    1.2714e5; 0.2444; 3.474e-5; 0.04391; 37.78; 0.311; 79.23;
    0.1212; 0.01052; 0.1378; 40.14;
    0.32; 47.13; 0.1; 0.08; 11];

    out = {'Priebe1998Beuckelmann',OpenState,p};
    return;
end
p = x(end-29:end);
x = x(1:end-30);

if V >= -40
    a_h = 0; a_j = 0;
    b_h = 1./(p(1).*(1+exp(-(V+p(2))./p(3))));
    b_j = p(4).*exp(-p(5).*V)./(1+exp(-p(6).*(V+p(7))));
else
	a_h = p(8).*exp(-(p(9)+V)./p(10));
    b_h = p(11).*exp(p(12).*V)+p(13).*exp(p(14).*V);
    a_j = (-p(15).*exp(p(16).*V)-p(17).*exp(-p(18).*V)).*(V+p(19))./(1+exp(p(20).*(V+p(21))));
    b_j = p(22).*exp(-p(23).*V)./(1+exp(-p(24).*(V+p(25))));
end
a_m = p(26).*(V+p(27))./(1-exp(-p(28).*(V+p(27))));
b_m = p(29).*exp(-V./p(30));
% if V >= -40
%     a_h = 0; a_j = 0;
%     b_h = 1./(0.13.*(1+exp(-(V+10.66)./11.1)));
%     b_j = 0.3.*exp(-2.535e-7.*V)./(1+exp(-0.1.*(V+32)));
% else
% 	a_h = 0.135.*exp(-(80+V)./6.8);
%     b_h = 3.56.*exp(0.079.*V)+3.1e5.*exp(0.35.*V);
%     a_j = (-1.2714e5.*exp(0.2444.*V)-3.474e-5.*exp(-0.04391.*V)).*(V+37.78)./(1+exp(0.311.*(V+79.23)));
%     b_j = 0.1212.*exp(-0.01052.*V)./(1+exp(-0.1378.*(V+40.14)));
% end
% a_m = 0.32.*(V+47.13)./(1-exp(-0.1.*(V+47.13)));
% b_m = 0.08.*exp(-V./11);

alpha = [a_h; a_j; a_m];
beta = [b_h; b_j; b_m];

% varargout = {23.*x(:,3).^3.*x(:,2).*x(:,1).*(V-ENa)};
if 1 == steadyNinfo % steady states
%     out = alpha./(alpha+beta);
    s = alpha./(alpha+beta);
    out = {[s; s(1).*s(2).*s(3).^3],1};
else
%     out = alpha.*(1-x(1:3))-beta.*x(1:3);
    dx = alpha.*(1-x(1:3))-beta.*x(1:3);
    out = [dx; dx(1).*x(2).*x(3).^3 + x(1).*dx(2).*x(3).^3 + x(1).*x(2).*3.*x(3).^2.*dx(3)];
end
