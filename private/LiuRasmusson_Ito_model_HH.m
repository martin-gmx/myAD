function out = LiuRasmusson_Ito_model_HH(t,x)
global Vclamp steadyNinfo
% No-block

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

OpenState = 3;

if 2 == steadyNinfo
    p = [0.03577; 3.8102; 0.06237; 4.5941;
        1.9; 13.5; 11.3; 0.051335; 0.067083];

    out = {'LR_HH',OpenState,p};
    return;
end
p = x(end-8:end);
x = x(1:end-9);

alpha_a = 1e-3.*exp(p(1)*V+p(2));
beta_a = 1e-3.*exp(-p(3)*V+p(4));

beta_i = 1e-3.*p(5).*exp((V+p(6))./p(7))./(1+p(8).*exp((V+p(6))./p(7)));
alpha_i = 1e-3.*p(5).*exp(-(V+p(6))./p(7))./(1+p(9).*exp(-(V+p(6))./p(7)));

% Block dynamics
% Kon = 4.5; Koff = 3.4;

%% HH-model
a_tau = 1./(alpha_a+beta_a);
a_inf = alpha_a.*a_tau;

i_tau = 1./(alpha_i+beta_i);
i_inf = alpha_i.*i_tau;

if 1 == steadyNinfo % steady states
    out = {[a_inf; i_inf; a_inf.^3.*i_inf],1};
else
    dx = [(a_inf-x(1))./a_tau; (i_inf-x(2))./i_tau];
    out = [dx; 3.*x(1).^2.*x(2).*dx(1) + x(1).^3.*dx(2)];
end
