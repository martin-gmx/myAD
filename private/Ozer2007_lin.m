function out = Ozer2007_lin(t,x)
global Vclamp steadyNinfo
% No-block

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

OpenState = 1;

if 2 == steadyNinfo
    p = [0.01383; 0.004403; 60; 0.2234; 0.01255];

    out = {'Ozer07 lin',OpenState,p};
    return;
end
p = x(end-4:end);
x = x(1:end-5);

alpha_a = p(1).*exp(p(2).*(V+p(3)));
beta_a  = p(4).*exp(-p(5).*(V+p(3)));

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
