function out = Ozer2007_nonlin2(t,x)
global Vclamp steadyNinfo
% No-block

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

OpenState = 1;

if 2 == steadyNinfo
    p = [2.3; -24.39; 16.98; 1.414; 1.555; 29.07;
        0.3092; 80.84; 18.59; 0.1755; 38.78; 44.33];

    out = {'Ozer07 nonlin2',OpenState,p};
    return;
end
p = x(end-11:end);
x = x(1:end-12);

alpha_a = p(1).*exp(-( (V+p(2))./p(3) ).^2) + p(4).*exp(-( (V+p(5))./p(6) ).^2);
beta_a  = p(7).*exp(-( (V+p(8))./p(9) ).^2) + p(10).*exp(-( (V+p(11))./p(12) ).^2);

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
