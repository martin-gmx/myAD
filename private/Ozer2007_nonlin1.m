function out = Ozer2007_nonlin1(t,x)
global Vclamp steadyNinfo
% No-block

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

OpenState = 1;

if 2 == steadyNinfo
    p = [0.03518; 0.1167; 60; 0.0008729;
        0.2102; 0.004187; 0.0001557];

    out = {'Ozer07 nonlin1',OpenState,p};
    return;
end
p = x(end-6:end);
x = x(1:end-7);

alpha_a = p(1).*exp(-( p(2).*(V+p(3)) + p(4).*(V+p(3)).^2 ));
beta_a  = p(5).*exp(-( p(6).*(V+p(3)) + p(7).*(V+p(3)).^2 ));

% Block dynamics
% Kon = 4.5; Koff = 3.4;

%% HH-model
a_tau = 1./(alpha_a+beta_a);
a_inf = alpha_a.*a_tau;

if 1 == steadyNinfo % steady states
    out = {a_inf, 1};
else
    dx = (a_inf-x(1))./a_tau;
    out = dx;
end
