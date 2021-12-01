function out = tenTusscher2004Panfilov(t, x)
global Vclamp steadyNinfo

V = interp1(Vclamp(:,1),Vclamp(:,2),t);

OpenState = 3;

if 2 == steadyNinfo
    p = [26; 7; 88; 24; 450; 45; 10; 6; 30; 11.5; 3; 60; 20; 1.12];

    out = {'tenTusscher2004Panfilov',OpenState,p};
    return;
end
p = x(end-13:end);
x = x(1:end-14);

Xr1_INF=1./(1.+exp((-p(1)-V)./p(2)));
Xr2_INF=1./(1.+exp((V-(-p(3)))./p(4)));
        
axr1=p(5)./(1.+exp((-p(6)-V)./p(7)));
bxr1=p(8)./(1.+exp((V-(-p(9)))./p(10)));
TAU_Xr1=axr1.*bxr1;

axr2=p(11)./(1.+exp((-p(12)-V)./p(13)));
bxr2=p(14)./(1.+exp((V-p(12))./p(13)));
TAU_Xr2=axr2.*bxr2;
% Xr1_INF=1./(1.+exp((-26.-V)./7.));
% Xr2_INF=1./(1.+exp((V-(-88.))./24.));
%         
% axr1=450./(1.+exp((-45.-V)./10.));
% bxr1=6./(1.+exp((V-(-30.))./11.5));
% TAU_Xr1=axr1.*bxr1;
% 
% axr2=3./(1.+exp((-60.-V)./20.));
% bxr2=1.12./(1.+exp((V-60.)./20.));
% TAU_Xr2=axr2.*bxr2;

if 1 == steadyNinfo % steady states
    out = {[Xr1_INF; Xr2_INF; Xr1_INF.*Xr2_INF],1};
%     out = [Xr1_INF; Xr2_INF];
else
%     out = [(Xr1_INF-x(1))./TAU_Xr1; (Xr2_INF-x(2))./TAU_Xr2];
    dx = [(Xr1_INF-x(1))./TAU_Xr1; (Xr2_INF-x(2))./TAU_Xr2];
    out = [dx; dx(1).*x(2) + x(1).*dx(2)];
end
