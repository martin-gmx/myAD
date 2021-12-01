function varargout = sensitivities(varargin)
%%% This package computes sensitivities and generalized sensitivities for ODEs
%%% It uses the myAD-package for automatic differentiation
%%%
%%% June/October, 2006, January, 2007, etc.
%%% Martin Fink
%%% martinfink@gmx.at
%
%
% Now you are able to set sens.sameFig to influence the figure numbers of the output
% to compare results of different systems.
%
%% calcSensitivities for ODEs or DDEs
% function sens = sensitivitiesODEs(functionDE, x0, param, stateStr, parStr, data, functionUpdateStates, lag, maxStep, justSolve)
% OR old function call: function sens = sensitivitiesODEs(functionDE, x0, param, stateStr, parStr, data, functionUpdateStates [, lag])
% OR old function call: function sens = sensitivitiesODEs(x0, stateStr, param, parStr, functionDE, data, functionUpdateStates [, lag])
%       x0          .. initial condition
%       stateStr    .. names of the states (cell array)
%       param       .. parameters for sensitivity analysis
%       parStr      .. names of parameters (cell array)
%       functionDE  .. function handle(!) for ODEs/DDEs - not only string!
%                       syntax: dx = function(t,x) or dx = function(t,x,xDelay)
%       data        .. data with time in first column
%                       OR times as column vector
%       functionUpdateStates .. function handle(!) for updating (only the) states at data points
%                       syntax: xNew = updateStatesAtData(n, xOld, data)
%                                   n .. idx of data point
%                                   xOld, xNew .. states and sensitivities
%                                   data .. data
%                       OR [] (empty matrix) if no update necessary
%       lag         .. provide delay-lag for DDEs or [] for ODEs
%       maxStep     .. maximal step length solving the system or [] for none
%       justSolve   .. -1 numerical sensitivity, 0 calc myAD sensitivity,
%                       1 just solve system
%
%       sens        .. output structure containing (all) information
%                       to be used to plot sensitivities and generalized sensitivities
%                       it's a good idea to save that structure
%
%       In all user functions you can use
%                       nState (number of states) and
%                       nParam (number of sensitivity parameters)
%
%
%% plotSensitivities
% function sensitivities(sens)
%       sens        ..  output from sensitivity calculation
%
%% calcGeneralizedSensitivity
% function gsCell = sensitivity(sens, idxCell, k)
%       sens        .. output from sensitivity calculation
%       idxCell     .. parameter choices to calculate gerenalized sensitivities for
%                       syntax: idxCell = {[1 2], [1 4 5], [3]}
%                           gen. sensitivity for 3 cases: two parameters (1 and 2), three param (1, 4, and 5), and one parameter (3)
%       k/B         .. idx of state to calculate and plot generalized sensitivities or
%                         observation matrix B thus measurements are B*x, where x are system states
%       gsCell      .. generalized sensitivities for cases in idxCell (cell array)
%
%% calcParamSubset
% function sensitivity(sens, k)
%       sens        .. output from sensitivity calculation
%       k/B         .. idx of state to calculate and plot generalized sensitivities or
%                         observation matrix B thus measurements are B*x, where x are system states
%

switch nargin
    case 11
        varargout{:} = calcSensitivitiesDEs(varargin{:});
    case {7, 8}
        varargout{:} = oldCalcSensitivitiesDEs(varargin{:});
    case 1
        varargout{:} = plotSensitivities(varargin{:});
    case 2
        varargout{:} = calcParamSubset(varargin{:});
    case 3
        varargout{:} = 0;
        %varargout{:} = calcGeneralizedSensitivity(varargin{:});
    otherwise
        error('No proper function call.');
end


%%
%% Backward compatibility
function sens = oldCalcSensitivitiesDEs(varargin)

    if nargin == 7
        varargin(end+1) = {[]};
    end

    if isa(varargin{1}, 'numeric')
        % OLD: x0, stateStr, param, parStr, functionDE
        % NEW: functionDE, x0, param, stateStr, parStr
        sens = calcSensitivitiesDEs(varargin{[5 1 3 2 4]}, varargin{6:end}, [], 0);
        warning('SENS:depricateInput','The order for the arguments has changed AND there are additional arguments for solving the system.\n         Please check the program sensitivities.m');
    else
        sens = calcSensitivitiesDEs(varargin{:}, [], 0);
        warning('SENS:depricateInput','There are additional arguments for solving the system.\n         Please check the program sensitivities.m');
    end

%%
%% Calculate Sensitivities
function sens = calcSensitivitiesDEs(functionDE, x0, param, stateStr, parStr, data, functionUpdateStates, lag, maxStep, justSolve, ICs)

global nState nParam nDelay time solution

x0 = x0(:);
param = param(:);
nState = length(stateStr);
nParam = length(parStr);
nDelay = length(lag);
s = size(data,1);
if length(x0)~=nState; error('SENS:INPUT','Length of initial states and state string does not match.'); end
if length(param)~=nParam; error('SENS:INPUT','Length of parameter vector and parameter string does not match.'); end
if s < 2; error('SENS:INPUT','Length of column of data/time vector smaller than 2.'); end

if isempty(lag)
    options = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep', maxStep);
%     options = odeset('MaxStep', maxStep);
else
    options = ddeset('MaxStep', maxStep);
end

if justSolve == 0 % use myAD
    funODE = @sensMyADode;
    funDDE = @sensMyADdde;
    if length(lag)>1
        funDDE = @sensMyMatrixADdde;
    end
    % initial condition for sensitivities
    if isempty(ICs)
        p0 = zeros(nParam*nState,1); % zero or
        x0 = [x0;p0];
    else
        x0 = [x0;ICs];
    end
else % solve numerically
    funODE = @solveDE;
    funDDE = @solveDE;
end

for repeat = 0:nParam
    time=data(1,1);
    solution=x0';
    dataIdx=1;

    for n=1:s(1)-1
        if isempty(lag)
            [t,x] = ode15s(funODE,data(n:n+1,1),x0,options,param,functionDE);
        else
            sol = dde23(funDDE,lag,@getDelayHist,data(n:n+1,1),options,param,functionDE);
            t = sol.x'; x = sol.y';
        end
        time=[time(1:end-1);t];
        solution=[solution(1:end-1,:);x];
        dataIdx = [dataIdx; length(time)];
        if ~(isempty(functionUpdateStates))
            solution(end,1:nState) = feval(functionUpdateStates, n+1, solution(end,1:nState), data);
        end
        x0 = solution(end,:)';
    end

    if justSolve >= 0; break; end % several loops only for numerical solution
    if repeat == 0
        origTime = time;
        origSolution = solution;
        origX0 = solution(1,:)';
        origParam = param;
        saveSolution = [];
%         saveSolution = zeros(length(time),nParam*nState);
    else
        solution = interp1(time, solution, origTime);
        saveSolution = [saveSolution, (solution-origSolution)/(param(repeat)-origParam(repeat))];
%         saveSolution(:,(repeat-1)*nState:repeat*nState) = (solution-origSolution)/(param(repeat)-origParam(repeat));
    end
    x0 = origX0; param = origParam;
    if repeat < nParam
        param(repeat+1) = origParam(repeat+1)*1.15; % increase 15%
    end
end

if justSolve < 0 % Rename variables
    time = origTime;
    solution = [origSolution, saveSolution]; clear saveSolution;
end

% add delays to parameters
% TODO
if 0;%~isempty(lag)
%     sol = ode15s(@delaySensDDE,time([1 end]),zeros(nDelay*nState,1),options,lag,param,functionDE);
    sol = dde23(@delaySensDDE,[lag, 2*lag],zeros(nDelay*nState,1),time([1 end]),options,lag,param,functionDE);
    figure(55); clf;
    for i = 1:3
        subplot(3,1,i);
        plot(sol.x, sol.y(i+[0 3],:)); legend('h1','h2'); xlabel(['x' num2str(i)]);
    end
%     figure(56); clf;
%     plot(time, -solution(:,1)+solution(:,2));
end

sens.t = time;
sens.x = solution;
sens.stateStr = stateStr;
sens.parStr = parStr;
sens.param = param;
sens.data = data;
sens.dataIdx = dataIdx;
sens.sameFig = 1;
sens.selectParam = [];

if justSolve > 0
    figure(justSolve*5231); clf;
    sy = ceil(sqrt(nState/2));
    sx = ceil(nState/sy);
    for k = 1:nState
        subplot(sx,sy,k)
        plot(time, solution(:,k));
        ylabel(stateStr{k});
    end
end

%%
%% Solve ODE/DDE
function dx = solveDE(t, x, varargin)

if nargin > 4 % delay
    dx = feval(varargin{end}, t, [x; varargin{end-1}], varargin{1});
else
    dx = feval(varargin{end}, t, [x; varargin{end-1}]);
end

%%
%% History function for the delay
function states = getDelayHist(t, varargin)

global time solution

if t <= time(1)
    states = solution(1,:)';
else
    states = interp1(time, solution, t)';
end

%%
%% Include sensitivity to delays
function  dx = delaySensDDE(t, xsens, xsensDel, lags, param, functionDDE)

global nState nParam nDelay time solution
% TODO

state = interp1(time, solution(:,1:nState), t)';
for i = 1:nDelay
    if t-lags(i) <=0
        stateDel(:,i) = zeros(nState,1);
    else
        stateDel(:,i) = interp1(time, solution(:,1:nState), t-lags(i))';
    end
    if t-2*lags(i) <=0
        state2Del(:,i) = zeros(nState,1);
    else
        state2Del(:,i) = interp1(time, solution(:,1:nState), t-2*lags(i))';
    end
end

gAD=[state;param]; %% [ current states ; parameters ]

%%% Solve derivatives with respect to time
gTime=myMatrixAD(t);
ffTime=feval(functionDDE, gTime, gAD, stateDel);
dfdTime=getderivs(ffTime);

%%% Solve derivatives with respect to states and parameters (=df)
gp=[myAD(state);param];
ff=feval(functionDDE, t, gp, stateDel);
dfdstate=getderivs(ff);

%%% Solve derivatives with respect to delayed states (=dfDelay)
gpDelay=myMatrixAD(stateDel);
ffDelay=feval(functionDDE, t, gAD, gpDelay);
dfdstateDelay=getderivs(ffDelay);

%%% Vector of derivatives for
for j = 1:nDelay
    idxA = (j-1)*nState;
    for k = 1:nState   %% b) the sensitivities
        dx(idxA+k,1) = 0;%-xsens(idxA+k); % -Z1DELAY?
        for i = j;%1:nDelay
%             dx(idxA+k) = dx(idxA+k) - dfdstateDelay(k,(1:nState)+(i-1)*nState)*state((1:nState));
            dx(idxA+k) = dx(idxA+k) - dfdstateDelay(k,(1:nState)+(i-1)*nState)*stateDel((1:nState))';
%             dx(idxA+k) = dx(idxA+k) - dfdstateDelay(k,(1:nState)+(i-1)*nState)*dfdTime((1:nState))';
%             dx(idxA+k) = dx(idxA+k) + dfdstate(k,(1:nState))*xsensDel(idxA+(1:nState),i);
%             dx(idxA+k) = dx(idxA+k) + dfdstateDelay(k,(1:nState)+(i-1)*nState)*xsens(idxA+(1:nState));
%             dx(idxA+k) = dx(idxA+k) - xsens(idxA+k);
        end
    end
end



%%
%% Sensitivity ODE
function dxAll = sensMyADode(t, xAll, param, functionODE)

global nState nParam
gAD=[xAll(1:nState);param];
gp=myAD(gAD,eye(nState+nParam));
FF=feval(functionODE, t, gp);
dxAll=getvalue(FF);
Jacobian=squeeze(getderivs(FF)); %squeeze puts vectors of derivatives into a matrix, i.e. squeezes out the 3rd dimension being 1

for j = 1:nParam
    for k = 1:nState
        idxA = nState + (j-1)*nState;
        dxAll(idxA+k,1) = Jacobian(k,1:nState)*xAll(idxA + (1:nState)) + Jacobian(k,nState+j);
    end
end

%%
%% Sensitivity DDE
function dxAll = sensMyADdde(t, xAll, yDelayAll, param, functionDDE)

global nState nParam
%%% Solve derivatives with respect to states and parameters (=df)
gAD=[xAll(1:nState);param]; %% [ current states ; parameters ]
gp=myAD(gAD);
ff=feval(functionDDE, t, gp, yDelayAll);
df=getderivs(ff);

%%% Solve derivatives with respect to delayed states (=dfDelay)
gpDelay=myAD(yDelayAll(1:nState));
ffDelay=feval(functionDDE, t, gAD, gpDelay);
dfDelay=getderivs(ffDelay);

%%% Vector of derivatives for
dxAll = getvalue(ff);  %% a) the model
for j = 1:nParam
    for k = 1:nState    %% b) the sensitivities
        idxA = nState + (j-1)*nState;
        dxAll(idxA+k) = df(k,1:nState)*xAll(idxA + (1:nState)) + df(k,nState+j) + dfDelay(k,1:nState)*yDelayAll(idxA + (1:nState));
    end
end

%%
%% Sensitivity DDE: Nr delays > 1
function dxAll = sensMyMatrixADdde(t, xAll, yDelayAll, param, functionDDE)

global nState nParam
nDelay = size(yDelayAll,2);

%%% Solve derivatives with respect to states and parameters (=df)
gAD=[xAll(1:nState);param]; %% [ current states ; parameters ]
gp=myAD(gAD);
ff=feval(functionDDE, t, gp, yDelayAll);
df=getderivs(ff);

%%% Solve derivatives with respect to delayed states (=dfDelay)
gpDelay=myMatrixAD(yDelayAll(1:nState,:));
ffDelay=feval(functionDDE, t, gAD, gpDelay);
dfDelay=getderivs(ffDelay);

%%% Vector of derivatives for
dxAll = getvalue(ff);  %% a) the model
for j = 1:nParam
    for k = 1:nState   %% b) the sensitivities
        idxA = nState + (j-1)*nState;
        dxAll(idxA+k) = df(k,1:nState)*xAll(idxA + (1:nState)) + df(k,nState+j);
        for i = 1:nDelay
            dxAll(idxA+k) = dxAll(idxA+k) + dfDelay(k,(1:nState)+(i-1)*nState)*yDelayAll(idxA + (1:nState),i);
        end
    end
end

%%
%% Plot Sensitivities
function varargout = plotSensitivities(sens)

nState = length(sens.stateStr);
nParam = length(sens.parStr);
if isfield(sens, 'sameFig'); sameFig = sens.sameFig; else sameFig = 1; end    
if isfield(sens, 'normSens'); normSens = sens.normSens; else normSens = 1; end    
if (isfield(sens, 'selectParam') && ~isempty(sens.selectParam)); selectParam = sens.selectParam(:); else selectParam = 1:nParam; end    

% figure(sameFig*5331); clf;
% sy = ceil(sqrt(nState/2));
% sx = ceil(nState/sy);
% for k = 1:nState
%     subplot(sx,sy,k)
%     plot(sens.t, sens.x(:,k));
%     ylabel(sens.stateStr{k});
% end
nrOfFigs = ceil(nState/12);
for idx = 1:nrOfFigs
    figure(sameFig*5330+idx); clf;
    sx = 4; sy = 3;
    if idx < nrOfFigs
        nPanels = 12;
    else
        nPanels = mod(nState,12);
    end
    for k = 1:nPanels
        kIdx = 12*(idx-1)+k;
        subplot(sx,sy,k)
        plot(sens.t, sens.x(:,kIdx));
        ylabel(sens.stateStr{kIdx});
    end
end

%%% idx for dx_k/dp_j
%%% idx = nState + (j-1)*nState + k;
j=selectParam(:)';
for k = 1:nState
    figure(sameFig*5431+k); clf; hold on;
    idx = nState + (j-1)*nState + k;
    temp = sens.x(:,idx).*sens.param(ones(length(sens.t),1)*j);
    maxY = max(abs(sens.x(:,k))); dt = diff(sens.t([1 end]));
    GTIS = trapz(sens.t, abs(temp))/(dt*maxY);
    GTIS2 = sqrt(trapz(sens.t, (temp).^2))/(dt*maxY);

    if normSens
        temp = temp./sens.x(:,k*ones(size(j)));
    end
    h = plot(sens.t, temp);
    for ii = 8:length(h)
        set(h(ii),'LineStyle','-.');
    end
    legend(sens.parStr{j},2);
    ylabel([sens.stateStr{k} ' - sensitivity']);
    xlabel('t (xx)');
    
    fprintf('\nSensitivity at t=%3.2f and GTIS and GTIS(2-norm) of %s to the parameters:\n', sens.t(end), sens.stateStr{k});
    for i = 1:length(j)
        fprintf('\n%12s: %12.2e %12.2e %12.2e', sens.parStr{j(i)}, temp(end,i), GTIS(i), GTIS2(i));
    end
    fprintf('\n');
    varargout{k,1} = temp(end,:);
    varargout{k,2} = GTIS(:);
end


%%
%% Generalized Sensitivity
function gsCell = calcGeneralizedSensitivity(sens, idxCell, k)

parStr = sens.parStr;
stateStr = sens.stateStr;
nParam = length(parStr);
nState = length(stateStr);
if max([idxCell{:}]) > nParam; error('SENS:INPUT','Index in idxCell is larger than number of parameters.'); end
if isfield(sens, 'sameFig'); sameFig = sens.sameFig; else sameFig = 1; end
if isfield(sens, 'normSens'); normSens = sens.normSens; else normSens = 1; end    

%%% use data points for generalized sensitivity
%%% OR use simulation points
if (size(sens.data,1)>2)
    t = sens.t(sens.dataIdx);
    x = sens.x(sens.dataIdx,:);
else
    t = sens.t;
    x = sens.x;
end

%%% normalize sensitivities to functions (to parameters unnecessary)
if 0%(normSens)
    j=1:nParam;
    for i = 1:nState
        idx = nState + (j-1)*nState + i;
        if normSens
            x(:,idx) = x(:,idx)./x(:,i*ones(size(j)));
        end
    end
end

dim = length(idxCell);

%%% Plot generalized sensitivities either for one state
%%% OR for two or more states (concatenate)
if length(k) < 2
    figure(sameFig*5431+k); clf;
    figure(sameFig*5433+k); clf;
    sx = ceil(sqrt(dim/2));
    sy = ceil(dim/sx);
    redX = x(:, (1:nParam)*nState + k);
    gsCell = cell(dim,1);
    for i = 1:dim
        [gsCell{i}, incGS] = doCalculationGeneralizedSensitivity(redX(:,idxCell{i}));
    figure(sameFig*5433+k);
        subplot(sx,sy,i);
        plot(t, incGS);
        legend(parStr{idxCell{i}},'Location','SouthEast');
    figure(sameFig*5431+k);
        subplot(sx,sy,i);
        plot(t, gsCell{i}');
        legend(parStr{idxCell{i}},'Location','SouthEast');
    end
    subplot(sx,sy,1);
    xlabel('t (min)');
    ylabel([stateStr{k} ' - gen. sensitivity']);
else
    B = k;
    k = find(sum(abs(k),1));
    figure(sameFig*5531+prod(k)); clf;
    figure(sameFig*6531+prod(k)); clf;
    kdim = size(B,1);
    sx = dim;
    sy = kdim;
    for i = 1:nParam
        measure(:,:,i) = B*x(:, i*nState + (1:nState))';
    end
    Jac = reshape(permute(measure,[2 1 3]), [], nParam);
    
    nGS = length(t)-1;
    for i = 1:dim
        [gsCell{i} incGS] = doCalculationGeneralizedSensitivity(Jac(:,idxCell{i}));
        gsCell{i} = gsCell{i}';
        maxGS = max(max(gsCell{i})); minGS = min(min(gsCell{i}));
        figure(sameFig*6531+prod(k));
        for j = 1:kdim
            subplot(sx,sy,(i-1)*kdim+j);
            plot(t(2:end), incGS(:,(1:nGS) + nGS*(j-1)));
        end
        legend(parStr{idxCell{i}},'Location','SouthEast');
        figure(sameFig*5531+prod(k));
        for j = 1:kdim
            subplot(sx,sy,(i-1)*kdim+j);
            plot(t(2:end), gsCell{i}((1:nGS) + nGS*(j-1),:));
            ylim([minGS-.1, maxGS+.1]);
        end
        subplot(sx,sy,(i-1)*kdim+1);
        legend(parStr{idxCell{i}},'Location','NorthWest');

%         outName = 'Gen. Sens.: ';
%         calcNoutMatrixRanking(gsCell{i}, B, parStr, stateStr, idxCell{i}, outName)
    end

    for j = 1:kdim
        subplot(sx,sy,j);
        addStrName = find(abs(B(j,:)));
        outName = stateStr{addStrName(1)};
        for i = addStrName(2:end)
            outName = [outName ' + ' stateStr{i}];
        end
        ylabel(outName);
    end
end

%%
%% Do Calculation of Sensitivities
function [gs, incgs] = doCalculationGeneralizedSensitivity(sensmat)
%%% sensmat..sensitivities of one state with respect to parameters
%%%          size: data points vs. nParam
s = size(sensmat);
S0=0;
for j=1:s(1)%-1
    S0=S0+sensmat(j,:)'*sensmat(j,:);
end
SI=pinv(S0);
    v1=SI*sensmat(1,:)';
    v2=v1.*sensmat(1,:)';
    incgs(:,1) = v2;
for i=2:s(1)%-1
    v1=SI*sensmat(i,:)';
    v2=v1.*sensmat(i,:)';
    incgs(:,i) = v2;
end
gs = cumsum(incgs,2);

%% Initialize Calculation of Parameter Subset Analysis
%% of Basic Sensitivities
function out = calcParamSubset(sens, B)
%%% System states x
%%% Measurements are B*x

parStr = sens.parStr;
stateStr = sens.stateStr;
nParam = length(parStr);
nState = length(stateStr);
if isfield(sens, 'selectParam') && ~isempty(sens.selectParam); selectParam = sens.selectParam(:); else selectParam = (1:nParam)'; end    
if isfield(sens, 'normSens'); normSens = sens.normSens; else normSens = 1; end

%%% use data points OR use simulation points
if (size(sens.data,1)>2)
    x = sens.x(sens.dataIdx,:);
else
    x = sens.x;
end

%%% normalize sensitivities to functions and parameters
if (normSens)
    j=selectParam(:)';
    for k = 1:nState
        idx = nState + (j-1)*nState + k;
        if normSens
            x(:,idx) = x(:,idx).*sens.param(ones(length(sens.t),1)*j)./(x(:,k*ones(size(j)))+eps);
        end
    end
end

%%% if provide index of 1 state k instead of observation matrix B
if length(B) < 2
    k = B;
    B = zeros(1, nState);
    B(1,k) = 1;
end

for i = 1:nParam
    measure(:,:,i) = B*x(:, i*nState + (1:nState))';
end
Jac = reshape(permute(measure,[2 1 3]), [], nParam);

%%% keep selected
Jac = Jac(:,selectParam);

%%% Calculate and output subset selection
outName = '';
out = calcNoutMatrixRanking(Jac, B, parStr, stateStr, selectParam, outName);

%%
%% Calculate and output subset selection
function out = calcNoutMatrixRanking(Jac, B, parStr, stateStr, selectParam, outName)
parStr = parStr(selectParam);
nSelParam = length(selectParam);

%%% The Calculation
temp = Jac'*Jac;
[Q,R,idx] = qr(temp,0);
% eigVal = sort(eig(temp), 'descend');
[eigVec, eigVal] = eig(temp);
[eigVal, srt] = sort(diag(eigVal), 'descend');
eigVec = eigVec(:, srt);

%%% Calc ranks
% thisRanks = [rank(temp), max(find(eigVal > eps^.25*eigVal(1)))];
thisRanks = [rank(temp), find(eigVal > 1e-6*eigVal(1),1,'last')];

%%% Output name
fillStr = '';
for j = 1:size(B,1)
    if j > 1; fillStr = ' & '; end
    addStrName = find(abs(B(j,:)));
    outName = [outName fillStr stateStr{addStrName(1)}];
    for i = addStrName(2:end)
        outName = [outName ' + ' stateStr{i}];
    end
end

%%% Associating parameters with the eigenvalues
fprintf('\nResult for measuring: %s \n\n', outName);
fprintf('Parameter (Nr) - Eigenvalue \n');
fprintf('====================== \n');
for i=1:nSelParam
    fprintf(' %8s (%2d) - %8.5e \n', parStr{idx(i)}, selectParam(idx(i)), eigVal(i));
end
fprintf('\n');
fprintf('%d/%d/%d. Matrix size/rank of Jacobian/numerical rank.\n', length(temp), thisRanks);
fprintf('Assemble parameters of different magnitudes for identification\nand check with generalized sensitivity.\n\n');

% fprintf('The according eigenvectors are given by:\n');
% disp(eigVec);
% 
% fprintf('The according SVD decomposition (U*S*V'') has the following V:\n');
% [u,s,v]=svd(Jac);
% disp(v);
warning('off','MATLAB:divideByZero');
condNr = eigVal(1)/eigVal(end);
warning('on','MATLAB:divideByZero');
fprintf('The condition number is: %4.4e.\n', condNr);

% calcParamClustering(Jac, thisRanks(2), selectParam);
% glueParamsToTop(Jac, thisRanks(2), selectParam, idx);

out = {length(temp), thisRanks, condNr};

%%
function calcParamClustering(Jac, numRank, selectParam)

nParam = size(Jac,2);

% Correlation is equivalent to cosine of vectors, if vectors are centered
% (i.e., v - mean(v) )
for Idx = 1:nParam
    Jac(:,Idx) = Jac(:,Idx)-mean(Jac(:,Idx));
end

A = zeros(nParam,nParam);
for Idx = 1:nParam % A = matrix of angles between sensitivity vectors
    for Jdx = Idx+1:nParam
        A(Idx,Jdx) = (Jac(:,Idx)'*Jac(:,Jdx))/(norm(Jac(:,Idx))*norm(Jac(:,Jdx)));
    end
end
A = abs(A);
% disp(A);

% toReduce = nParam - numRank; % Groups to build
groups = {};
for i = 1:nParam%toReduce
    [temp, iii] = max(A(:)); % smallest angle / largest cos(alpha)
    if temp <= 0.94; break; end
    Idx = mod(iii, nParam);
    Jdx = (iii-Idx)/nParam+1;
    
%     disp('Next.');
%     groups{:}
%     Idx
%     Jdx
    
    newGroup = 1;
    for j = 1:length(groups)
        if newGroup == 1
            if max(ismember(groups{j}, Idx)) || max(max(ismember(groups{j}, Jdx)))
            % if Idx or Jdx part of group then connect
                if max(ismember(groups{j}, Idx))
                    groups{j}(end+1) = Jdx;
                else
                    groups{j}(end+1) = Idx;
                end
            % set all angles within group to 0
                for k = 1:length(groups{j})-1
                    A(groups{j}(k),groups{j}(end)) = 0;
                    A(groups{j}(end),groups{j}(k)) = 0;
                end
            % now we have found a group, but have to look for connection
            % with other group
                newGroup = 0;
                firstGroupIdx = j;
                lookForPar = groups{j}(end);
            end
        else%if j > firstGroupIdx
            if max(ismember(groups{j}, lookForPar))
            % found second index in other group
            % have to connect groups
            % and clean connections in matrix A
                for clearA = 1:length(groups{firstGroupIdx})
                    for clearB = 1:length(groups{j})
                        A(groups{firstGroupIdx}(clearA),groups{j}(clearB)) = 0;
                        A(groups{j}(clearB),groups{firstGroupIdx}(clearA)) = 0;
                    end
                end
                groups{firstGroupIdx} = unique([groups{firstGroupIdx} groups{j}]);
                groups(j) = []; % delete second group
            % does not count as "reduction"
%                 i = i-1;
                break;
            end
        end
    end
    if newGroup == 1
        groups{end+1} = [Idx, Jdx];
        A(Idx,Jdx) = 0;
        A(Jdx,Idx) = 0;
    end
end
for i = 1:length(groups)
    for j = 1:length(groups{i})
        fprintf('%4d ', selectParam(groups{i}(j)));
    end
    fprintf('\n');
end
if ~isempty(groups); disp(temp); end

%%
function glueParamsToTop(Jac, numRank, selectParam, ordIdx)

nParam = size(Jac,2);
topParam = ordIdx(1:numRank);
bottomParam = ordIdx(numRank+1:end);

% Correlation is equivalent to cosine of vectors, if vectors are centered
% (i.e., v - mean(v) )
for Idx = 1:nParam
    Jac(:,Idx) = Jac(:,Idx)-mean(Jac(:,Idx));
end

A = zeros(numRank,nParam-numRank);
for Idx = 1:numRank % A = matrix of angles between sensitivity vectors
    topJ = Jac(:,topParam(Idx));
    topJnorm = norm(topJ);
    if topJnorm == 0; continue; end
    
    for Jdx = 1:nParam-numRank
        botJ = Jac(:,bottomParam(Jdx));
        botJnorm = norm(botJ);
        if botJnorm == 0; continue; end
        A(Idx,Jdx) = (topJ'*botJ)/(norm(topJ)*norm(botJ));
    end
end
A = abs(A);
% disp(A);

% Glue the "not identifiable" parameters to the top numrank parameters
[a, gIdx] = max(A);

for Idx = 1:numRank
    groups{Idx} = selectParam([topParam(Idx) bottomParam(gIdx==Idx)])';
    for Jdx = 1:length(groups{Idx})
        fprintf('%4d ', groups{Idx}(Jdx));
    end
    fprintf('\n');
end


