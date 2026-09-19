function [LPproblem, varMap, buildInfo] = CAFBA_buildLP(model, varargin)
% Build the linear-programming problem used by CAFBA without solving it.
%
% This function contains the LP-construction part of the original
% CAFBA_OptimizeCbModel_glpk implementation. Both CAFBA optimization and
% CAFBA-flux variability analysis should call this function so that they
% use exactly the same mathematical problem.
%
% USAGE
%
%   [LPproblem,varMap,buildInfo] = CAFBA_buildLP(model)
%
%   [LPproblem,varMap,buildInfo] = CAFBA_buildLP(model, ...
%       'optSense','max', ...
%       'optMode','flux', ...
%       'optProtGroup',2, ...
%       'constraintSense','U', ...
%       'verbose',0);
%
% INPUT
%
%   model
%       COBRA model with the additional CAFBA fields:
%           model.w
%           model.protGroup
%
% OPTIONAL INPUTS
%
%   'optSense'
%       'max' or 'min'. Default: 'max'
%
%   'optMode'
%       'flux': optimize model.c under the CAFBA allocation constraint.
%       'phi' : optimize one proteome group without the global CAFBA row.
%       Default: 'flux'
%
%   'optProtGroup'
%       Proteome-group index used when optMode='phi'. Default: 2
%
%   'constraintSense'
%       GLPK sense of the global CAFBA allocation constraint:
%           'S' equality
%           'U' upper bound, <=
%           'L' lower bound, >=
%       Default: 'S', matching the original Mori implementation.
%
%   'verbose'
%       GLPK message level. Default: 1
%
% OUTPUT
%
%   LPproblem
%       Complete GLPK LP structure:
%           A, b, c, lb, ub, ctype, vartype, sense, param
%
%   varMap
%       Mapping between extended CAFBA variables and original reaction
%       net fluxes. In particular:
%
%           vNet = varMap.T * xExtended
%
%   buildInfo
%       Indices and dimensions used to construct the extended problem.
%

% -------------------------------------------------------------------------
% Defaults
% -------------------------------------------------------------------------

optSense = 'max';
optMode = 'flux';
optProtGroup = 2;
constraintSense = 'U';
verbose = 1;

% -------------------------------------------------------------------------
% Parse name-value inputs
% -------------------------------------------------------------------------

if mod(length(varargin),2) ~= 0
    error('Optional arguments must be supplied as name-value pairs.');
end

for k = 1:2:length(varargin)
    switch varargin{k}
        case 'optSense'
            optSense = varargin{k+1};
        case 'optMode'
            optMode = varargin{k+1};
        case 'optProtGroup'
            optProtGroup = varargin{k+1};
        case 'constraintSense'
            constraintSense = varargin{k+1};
        case 'verbose'
            verbose = varargin{k+1};
        otherwise
            error('Optional argument not recognized: %s',varargin{k});
    end
end

% -------------------------------------------------------------------------
% Validate model and options
% -------------------------------------------------------------------------

requiredFields = {'S','lb','ub','c','w','protGroup'};

for k = 1:length(requiredFields)
    if ~isfield(model,requiredFields{k})
        error('Invalid CAFBA model: ''%s'' field missing.', ...
            requiredFields{k});
    end
end

[nMets,nRxns] = size(model.S);

if length(model.lb) ~= nRxns || length(model.ub) ~= nRxns || ...
        length(model.c) ~= nRxns || length(model.w) ~= nRxns
    error('model.lb, model.ub, model.c and model.w must match model.S.');
end

if strcmp(optSense,'max')
    lpSense = -1;
elseif strcmp(optSense,'min')
    lpSense = 1;
else
    error('optSense must be ''max'' or ''min''.');
end

if ~(strcmp(optMode,'flux') || strcmp(optMode,'phi'))
    error('optMode must be ''flux'' or ''phi''.');
end

if ~(strcmp(constraintSense,'S') || ...
        strcmp(constraintSense,'U') || ...
        strcmp(constraintSense,'L'))
    error('constraintSense must be ''S'', ''U'' or ''L''.');
end

if strcmp(optMode,'phi')
    if ~isscalar(optProtGroup) || optProtGroup < 1 || ...
            optProtGroup > length(model.protGroup) || ...
            optProtGroup ~= round(optProtGroup)
        error('optProtGroup is outside the valid proteome-group range.');
    end
end

% -------------------------------------------------------------------------
% Right-hand side and metabolite constraint senses
% -------------------------------------------------------------------------

if isfield(model,'b')
    bVec = model.b(:);
else
    bVec = zeros(nMets,1);
end

if length(bVec) ~= nMets
    error('model.b must contain one entry per metabolite constraint.');
end

ctypeVec = repmat('S',nMets,1);

% -------------------------------------------------------------------------
% Classify weighted reactions by their original flux bounds
% -------------------------------------------------------------------------

positiveVec = [];
negativeVec = [];
mixedVec = [];

for r = 1:nRxns

    if model.w(r) ~= 0

        if model.ub(r) >= 0 && model.lb(r) >= 0
            positiveVec(end+1,1) = r; %#ok<AGROW>

        elseif model.ub(r) <= 0 && model.lb(r) < 0
            negativeVec(end+1,1) = r; %#ok<AGROW>

        else
            mixedVec(end+1,1) = r; %#ok<AGROW>
        end
    end
end

nPositive = length(positiveVec);
nNegative = length(negativeVec);
nMixed = length(mixedVec);

nRxnsExt = nRxns + nMixed;

if strcmp(optMode,'flux')
    nMetsExt = nMets + 1;
else
    nMetsExt = nMets;
end

% -------------------------------------------------------------------------
% Initialize extended variables
% -------------------------------------------------------------------------

LPproblem.lb = zeros(nRxnsExt,1);
LPproblem.ub = zeros(nRxnsExt,1);
LPproblem.c = zeros(nRxnsExt,1);
LPproblem.vartype = repmat('C',nRxnsExt,1);

LPproblem.lb(1:nRxns) = model.lb(:);
LPproblem.ub(1:nRxns) = model.ub(:);
LPproblem.c(1:nRxns) = model.c(:);

Smixed = sparse(nMets,nMixed);

for k = 1:nMixed

    r = mixedVec(k);
    reverseVar = nRxns + k;

    % The original variable becomes the non-negative forward variable.
    LPproblem.lb(r) = 0;

    % A new non-negative reverse variable is introduced.
    LPproblem.lb(reverseVar) = 0;
    LPproblem.ub(reverseVar) = -model.lb(r);

    % Preserve any objective coefficient on the net reaction flux.
    LPproblem.c(reverseVar) = -model.c(r);

    % Reverse variable has the opposite stoichiometric column.
    Smixed(:,k) = -model.S(:,r);
end

% -------------------------------------------------------------------------
% Build mapping from extended variables to original net reaction fluxes
% -------------------------------------------------------------------------

varMap.T = sparse(nRxns,nRxnsExt);
varMap.T(:,1:nRxns) = speye(nRxns);

for k = 1:nMixed
    varMap.T(mixedVec(k),nRxns+k) = -1;
end

varMap.nRxns = nRxns;
varMap.nRxnsExt = nRxnsExt;
varMap.forwardVarIndex = (1:nRxns)';
varMap.reverseVarIndex = zeros(nRxns,1);

for k = 1:nMixed
    varMap.reverseVarIndex(mixedVec(k)) = nRxns+k;
end

% -------------------------------------------------------------------------
% Build the CAFBA allocation row
% -------------------------------------------------------------------------

if strcmp(optMode,'flux')

    Sac = sparse(1,nRxnsExt);

    for k = 1:nPositive
        r = positiveVec(k);
        Sac(r) = model.w(r);
    end

    for k = 1:nNegative
        r = negativeVec(k);
        Sac(r) = -model.w(r);
    end

    for k = 1:nMixed
        r = mixedVec(k);
        Sac(r) = model.w(r);
        Sac(nRxns+k) = model.w(r);
    end

    phi0Sum = 0;

    for k = 1:length(model.protGroup)
        phi0Sum = phi0Sum + model.protGroup(k).phi0;
    end

    allocationRHS = 1 - phi0Sum;

    bVec(nMets+1,1) = allocationRHS;
    ctypeVec(nMets+1,1) = constraintSense;

else
    Sac = sparse(0,nRxnsExt);
    phi0Sum = NaN;
    allocationRHS = NaN;
end

% -------------------------------------------------------------------------
% Assemble the complete constraint matrix
% -------------------------------------------------------------------------

LPproblem.A = [[model.S,Smixed];Sac];
LPproblem.b = bVec;
LPproblem.ctype = ctypeVec;

% -------------------------------------------------------------------------
% Build the objective
% -------------------------------------------------------------------------

if strcmp(optMode,'phi')

    LPproblem.c = zeros(nRxnsExt,1);
    groupRxns = model.protGroup(optProtGroup).rxns(:);

    for k = 1:nPositive
        r = positiveVec(k);
        if any(groupRxns == r)
            LPproblem.c(r) = model.w(r);
        end
    end

    for k = 1:nNegative
        r = negativeVec(k);
        if any(groupRxns == r)
            LPproblem.c(r) = -model.w(r);
        end
    end

    for k = 1:nMixed
        r = mixedVec(k);
        if any(groupRxns == r)
            LPproblem.c(r) = model.w(r);
            LPproblem.c(nRxns+k) = model.w(r);
        end
    end
end

LPproblem.sense = lpSense;
LPproblem.param.msglev = verbose;

% -------------------------------------------------------------------------
% Diagnostic information
% -------------------------------------------------------------------------

buildInfo.nMets = nMets;
buildInfo.nMetsExt = nMetsExt;
buildInfo.nRxns = nRxns;
buildInfo.nRxnsExt = nRxnsExt;

buildInfo.positiveVec = positiveVec;
buildInfo.negativeVec = negativeVec;
buildInfo.mixedVec = mixedVec;

buildInfo.nPositive = nPositive;
buildInfo.nNegative = nNegative;
buildInfo.nMixed = nMixed;

buildInfo.phi0Sum = phi0Sum;
buildInfo.allocationRHS = allocationRHS;
buildInfo.constraintSense = constraintSense;
buildInfo.optMode = optMode;
buildInfo.optProtGroup = optProtGroup;
buildInfo.optSense = optSense;

% Preserve the exact primary objective used by the generated LP.
varMap.primaryObjective = LPproblem.c(:);

end
