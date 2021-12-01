function paramID(chooseModel, chooseVclamp, calcSens, useSavedResults)

if nargin < 4
    chooseModel = 6; % 1..16 for all used models
    chooseVclamp = 2; % 0 short, 1 long, 2 optimal
    calcSens = 3; %% 0..solve, 1..basic sens, 3..subset selection
    useSavedResults = 0;
end
calcCompTime = 0;

tRange   = [];%[0 3000; 6000 14500];% 28900 29500];
selectParam = [];%[1:4 6 8:12 14 16:18 20:24];%[27 1 30 14 3 28 2 9 21 12];

%% Do not change anything beneath
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doICaL = 0;
switch chooseModel
case 1 % INa
    modelfunc = @Priebe1998Beuckelmann; E = 65;
case 2
    modelfunc = @Clancy2002Rudy; E = 65;
case 3
    modelfunc = @myClancy2002Rudy; E = 65;
case 4 % ICaL
    modelfunc = @Faber2007Rudy; E = 150; doICaL = 1;
case 5
    modelfunc = @myFaber2007Rudy; E = 150; doICaL = 1;
case 6 % IKr
    modelfunc = @tenTusscher2004Panfilov; E = -86;
case 7
    modelfunc = @Mazhari2001Nuss; E = -86;
case 8
    modelfunc = @Wang1997Rasmusson; E = -86;
case 9
    modelfunc = @Clancy2001Rudy; E = -86;
case 10 % IKs
    modelfunc = @Silva2005Rudy; E = -86;
case 11 % Ito
    modelfunc = @Greenstein2000Winslow; E = -86;
case 12
    modelfunc = @LiuRasmusson_Ito_model_HH; E = -86;
case 13
    modelfunc = @LiuRasmusson_Ito_model_MM; E = -86;
case 14
    modelfunc = @myFaber2007Rudy2nd; E = 150; doICaL = 1;
case 15 % IKr
    modelfunc = @myWang1997Rasmusson; E = -86;
case 16 % INa
    modelfunc = @FinkINa; E = 65;

case 18 % Cav with sodium
    modelfunc = @Ozer2007_lin; E = 65;
case 19
    modelfunc = @Ozer2007_nonlin1; E = 65;
case 20
    modelfunc = @Ozer2007_nonlin2; E = 65;
case 23
    modelfunc = @Irvine1999Winslow; E = 65;
end

if calcCompTime; doICaL = 0; calcSens = 0; end

%% Prepare V-clamp
global Vclamp CassClamp steadyNinfo

[Vclamp00, addToName] = getVclamp(chooseVclamp);

if (doICaL) % 3 Cass-steps
    Vclamp = [Vclamp00; Vclamp00(2:end,:); Vclamp00(2:end,:)];
else
    Vclamp = Vclamp00;
end
Vclamp(:,1) = cumsum(Vclamp(:,1));
Vclamp = [0 Vclamp(1,2); Vclamp];

if isempty(tRange); tRange = [0 Vclamp(end,1)]; end
myeps = 1e-10;
CassClamp = [Vclamp(1,1)  1e-3; tRange(2)/3-myeps 1e-3; tRange(2)/3 5e-3; tRange(2)*2/3-myeps 5e-3; tRange(2)*2/3 15e-3; tRange(2) 15e-3];

%% Get information on the model
steadyNinfo = 2; % Info-mode
out = feval(modelfunc, 0, []);
modelName = out{1};
openStates = out{2};
param = out{3};
disp(modelName);

%% Get initial steady state
steadyNinfo = 1; % Steady-state mode
out = feval(modelfunc, 0, param);
x0 = out{1};
% condNr = out{2};
% whos x0
% return

%% Do sensitivity analysis?
if (calcSens)
    idxCell = {[1 4],[2 5],[3 6]};
    if chooseModel == 23 % Irvine1999Winslow
        Vclamp = []; % Too stiff for sensitivity analysis
    end

    if ~exist('Sens','dir'); mkdir('Sens'); end
    combineStates = 1; %% 0 .. normal generalized sensitivity for each state,
                       %% 1 .. combine information of TcCO2 and TcO2 (B*x)
    sameFig = chooseModel;
    nrStates = length(x0);
    nrParams = length(param);
    stateStr = num2cell([repmat('S_',nrStates,1) num2str((1:nrStates)',2)],2);
    parStr = num2cell([repmat('p_',nrParams,1) num2str((1:nrParams)',2)],2);
    name = ['Sens/' modelName addToName];
    functionDE = modelfunc;
    functionUpdateStates = [];
    tspan = Vclamp([1 end],1);
    lag = [];
    auxilFun = [];
    MaxStep = 10;
    B = zeros(1,nrStates);
    B(1,openStates) = 1;
    DEname = '';

    steadyNinfo = 1; % Steady-state mode
    parAD = myAD(param);
    out = feval(modelfunc, 0, parAD);
    xICs = out{1};
    ICs = getderivs(xICs);

    steadyNinfo = 0; % ODE-mode
    run_sensitivities(functionDE, x0, param, stateStr, parStr, tspan, functionUpdateStates, lag, MaxStep, calcSens, useSavedResults, sameFig, [name DEname], idxCell, B, combineStates, selectParam, tRange, 1, auxilFun, ICs(:));
    return;
end

%% Simulate V-clamp protocol
tic
options = odeset('MaxStep',10);
steadyNinfo = 0; % ODE-mode
[t,x] = ode15s(@runODE, Vclamp([1 end],1), x0, options, param, modelfunc);
toc
if calcCompTime; return; end
V = interp1(Vclamp(:,1),Vclamp(:,2),t);

fAdd = chooseModel*10;

figure(1+fAdd);
plot(t,x); ylabel('States'); xlabel('t (ms)');

figure(2+fAdd);
% subplot(2,1,1); 
plot(t,sum(x(:,openStates),2)); ylabel('Open Prob'); ylim([0 1]);
xlabel('t (ms)');
 
figure(6+fAdd);
% subplot(2,1,2);
current = [t, sum(x(:,openStates),2).*(V-E)];
plot(t,current(:,2)); ylabel('Current'); xlabel('t (ms)');
if chooseVclamp == 2; save(['Data/' modelName '-Vopt-current.dat'],'-ascii','current'); end

figure(3+fAdd);
plot(t,V); ylabel('V (mV)'); xlabel('t (ms)'); title('V clamp protocol');

figure(4+fAdd);
steadyNinfo = 1; %% Calculating steady state
nSS = 50;
Vclamp = [0 -100; 1 80]; CassClamp = [0 5e-3; 1 5e-3];
tSS = linspace(Vclamp(1,1),Vclamp(end,1),nSS);
Vss = interp1(Vclamp(:,1),Vclamp(:,2),tSS);
for i = 1:nSS
    out = feval(modelfunc, tSS(i), param);
    x0 = out{1};
    condNr(i) = out{2};
    currSS(i) = sum(x0(openStates));
    if condNr(i) == 1 % HH
        allSS(:,i) = x0;
    else % MM
        allSS(:,i) = [x0; 1-sum(x0)];
    end
end
plot(Vss, currSS.*(Vss-E)); ylabel('Current'); xlabel('V (mV)'); title('I-V steady state curve');

figure(5+fAdd); clf;
% subplot(2,1,1)
plot(Vss, allSS); ylabel('Steady states');
% subplot(2,1,2)
% plot(Vss, condNr); ylabel(['Condition Nr';' of Matrix  ']);
xlabel('V (mV)');

function dx = runODE(t, x, p, modelfunc)

dx = feval(modelfunc, t, [x; p]);
