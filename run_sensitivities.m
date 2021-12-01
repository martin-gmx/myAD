function varargout = run_sensitivities(functionDE, x0, param, stateStr, parStr, data, functionUpdateStates, lag, maxStep, calcSens, useSavedResults, sameFig, name, idxCell, B, combineStates, selectParam, tRange, normSens, auxilFun, ICs)
%% Function for easy use of the sensitivity analysis package.
%% It assumes that one can save/load the sensitivities in order to use the various post-processing options.

%% Parameters for solution of the IVP
% functionDE .. ODE or DDE function
% x0 .. initial values
% param .. vector of (all) parameters to consider
% stateStr .. cell array for names for x0
% parStr .. cell array for names for param
% data .. either array of data matrix (first column is time) - or column-vector for time(span)
% functionUpdateStates .. to include events one can update the states at the time points specified in data
% lag .. delays for differential equations
% maxStep .. maximum step size for DE-solver

%% Parameters/Options for calculating/post-processing of the sensitivities
% calcSens .. 0..solve, 1..basic sens, 2..generalized sens,
%             3..subset selection, -x..numerical sens
% useSavedResults .. use previously run sensitivity output
% sameFig .. 0: creates new figures, 1-...: indicates a specific set of figures (hopefully without overlap)
% name .. name for saving the sensitivity structure
% idxCell .. cell array for parameter sets for generalized sensitivities
% B .. observation matrix - observables are assumed to be linear combinations of the system states
% combineStates .. 0: show output for all states, 1: just for the observables
% selectParam .. for subset selection ([]: use all parameters, any other array shows output only for selected parameters)
% tRange .. restricts the output to specific time interval
% normSens .. normalize sensitivities (default) can be turned off by setting to 0 (if function y goes through zero)
% auxilFun .. funxtion handle for adding algebraic functions to the sensitivity analysis (e.g., a = x(1)+exp(x(2))./x(3))
% ICs .. initial conditions for sensitivities - otherwise IC=zeros.

% if no 'lag', then no delay differential equations
useODE = isempty(lag);
% set default values for some options
if ~exist('tRange','var'); tRange = []; end
if ~exist('normSens','var'); normSens = 1; end
if ~isempty(tRange);
    if (tRange(1,1) < data(1,1)) || (tRange(end,2) > data(end,1));
        error('Selected time range not within simulation time interval.');
    end
end

% just solve the system
if calcSens == 0
    if useODE
        varargout{:}= sensitivities(functionDE, x0, param, stateStr, parStr, data, functionUpdateStates, [], maxStep, sameFig, []);
    else
        varargout{:}= sensitivities(functionDE, x0, param, stateStr, parStr, data, functionUpdateStates, lag, maxStep, sameFig, []);
    end
    return;
end

% % % % calculate model subspaces of the states
% % % if calcSens == -2
% % %     varargout{:}= calcModelReduction(functionDE, x0, param, stateStr);
% % %     return;
% % % end

%% calcSensitivities

if ~useODE; name = [name 'Delay']; end
if calcSens < 0; name = [name '_NumDiff']; end
    
if ~useSavedResults % derive sensitivities
    %% Speed up Matlab stupid solvers by splitting solution into smaller pieces
    if (size(data,1)==2 && data(2,1)-data(1,1) > 10000)
        resetData = 1;
        oldData = data;
        data = (data(1,1):3000:data(2,1))';
    else
        resetData = 0;
    end

    if calcSens < 0
        solveOpt = -1; % solve sens numerically
    else
        solveOpt = 0;  % solve sens with myAD
    end
    tic
    if useODE
        sens = sensitivities(functionDE, x0, param, stateStr, parStr, data, functionUpdateStates, [], maxStep, solveOpt, ICs);
    else
        sens = sensitivities(functionDE, x0, param, stateStr, parStr, data, functionUpdateStates, lag, maxStep, solveOpt, []);
    end
    toc

    %% Speed up Matlab stupid solvers
    if resetData
        sens.data = oldData;
        sens.dataIdx = length(sens.t);
    end

    save([name '.sens'],'sens','-mat');
else % otherwise load old results
    load([name '.sens'],'-mat');
end
calcSens = abs(calcSens); % remove choice between numDiff and AD

% Restrict data range for analysis to interval tRange
if ~isempty(tRange)
    keep = sens.t >= tRange(1,1) & sens.t<= tRange(1,2);
    for j = 2:size(tRange,1)
        keep = keep | (sens.t >= tRange(j,1) & sens.t<= tRange(j,2));
    end
    sens.t = sens.t(keep);
    sens.x = sens.x(keep,:);
    if length(sens.data)>2
        keep = sens.data >= tRange(1,1) & sens.data<= tRange(1,2);
        for j = 2:size(tRange,1)
            keep = keep | (sens.data >= tRange(j,1) & sens.data<= tRange(j,2));
        end
        sens.data = sens.data(keep);
    else
        sens.data = [tRange(1,1) tRange(end,end)];
    end
end

% (re-)set parameters
sens.sameFig = sameFig;
sens.selectParam = selectParam;
sens.normSens = normSens;

% Add sensitivities for auxiliary variables
if exist('auxilFun','var') && ~isempty(auxilFun)
        sens = includeAuxiliaries(sens, auxilFun);
end

% Check if observability matrix is consistent with number of states.
if length(B)>1 && size(B,2)~=length(sens.stateStr)
    error('RUNSENS:Bmatrix','Size of measurement matrix B is not consistent with number of states.');
end

switch calcSens
case {-1, 1}
    %% plotSensitivities
    varargout{:} = sensitivities(sens);
case 2
    %% calcGeneralizedSensitivity
    if combineStates
        varargout{:} = sensitivities(sens, idxCell, B);
    else
        for k = 1:length(sens.stateStr)
            varargout{:} = sensitivities(sens, idxCell, k);
        end
    end
case 3
    %% calcParamSubset
    if combineStates
        varargout{:}= sensitivities(sens, B);
    else
        for k = 1:length(sens.stateStr)
            varargout{:}= sensitivities(sens, k);
        end
    end
end

%%
%% Function: Add auxiliary functions to state structure
%%
function sens = includeAuxiliaries(sens, auxilFun)

    nState = length(sens.stateStr);
    nParam = length(sens.parStr);
    nt = length(sens.t);
    parAD = myAD([zeros(nState,1); sens.param]);
    parAD = parAD(nState+1:end);

    [temp, auxilStr] = feval(auxilFun, 0, ones(nState+nParam,1), parAD);
    nAux = length(auxilStr);
    for i = 1:nt
        states = [sens.x(i,1:nState)'; sens.param];
        statesAD = myAD(states);
        result = feval(auxilFun, sens.t(i), statesAD, parAD);
        auxSens(i,1:nAux) = getvalue(result)';
        auxDer = getderivs(result);
        for j = 1:nParam
            for k = 1:nAux
                idxA = nAux + (j-1)*nAux;
                idxB = nState + (j-1)*nState;
                auxSens(i, idxA+k) = auxDer(k,1:nState)*sens.x(i, idxB + (1:nState))' + auxDer(k, nState+j);
            end
        end
    end
    % Now to merge with the old sens.x structure
    sens.stateStr = {sens.stateStr{:}, auxilStr{:}};
    % Convert to 3D
    temp = reshape(sens.x, nt, nState, nParam+1);
    temp2 = reshape(auxSens, nt, nAux, nParam+1);
    % Concat both tensors
    temp3 = cat(2, temp, temp2);
    % Convert back to 2D
    sens.x = reshape(temp3, nt, (nState+nAux)*(nParam+1));
