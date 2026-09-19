
function [lambda,v,phi_E,model_t] = dCAFBA_AADown_Glyc_iJR904_post(model,coreState,coreInput)

% Refactored post-shift CAFBA solver for glycerol + amino acids -> glycerol.
%
%
% Required coreState fields
%   phi_Rb   : ribosomal proteome fraction at the current step
%   ratios   : 14 x 1 vector of AA-group capacity ratios, ordered as
%              coreInput.groupNames
%
% Required coreInput fields
%   Dir              : qualitative reaction directions (0 reversible,
%                      1 irreversible-forward)
%   phi_C            : carbon-sector proteome fraction
%   phi_Q            : housekeeping proteome fraction
%   groupNames       : 14 AA-group names
%   groupRxns        : 14 x 1 cell; reaction-index vector for each group
%   steadyGroupFlux  : 14 x 1 cell; corresponding post-shift steady-state
%                      reaction fluxes for each group

nGroups = numel(coreInput.groupRxns);

if numel(coreState.ratios) ~= nGroups
    error('coreState.ratios must contain one value per AA group.');
end
if numel(coreInput.steadyGroupFlux) ~= nGroups
    error('coreInput.steadyGroupFlux must contain one vector per AA group.');
end
if numel(coreInput.groupNames) ~= nGroups
    error('coreInput.groupNames and coreInput.groupRxns must have the same length.');
end

% Global CAFBA proteome allocation constraint.
model.protGroup(4).phi0 = coreInput.phi_Q + coreInput.phi_C + coreState.phi_Rb;
phi_E = 1 - model.protGroup(4).phi0;

% Apply the same reaction-specific AA pathway capacity bounds as the
for g = 1:nGroups
    rxnIndices = coreInput.groupRxns{g};
    steadyFlux = coreInput.steadyGroupFlux{g};
    ratio = coreState.ratios(g);

    steadyFlux = steadyFlux(:);
    if numel(steadyFlux) ~= numel(rxnIndices)
        error('Group %s: steady-state flux vector length does not match reaction-index length.', ...
            coreInput.groupNames{g});
    end

    capacities = ratio .* abs(steadyFlux);
    model = setGroupBounds(model,rxnIndices,coreInput.Dir,capacities);
end

% Save the exact time-specific model passed to CAFBA (for CAFVA/FVA).
model_t = model;

sol = CAFBA_OptimizeCbModel_glpk(model_t,'constraintSense','U');
v = sol.x;
lambda = sol.f;

end


function model = setGroupBounds(model,rxnIndices,Dir,capacities)
% Apply reaction-specific bounds using exactly the original direction rule.

capacities = capacities(:);

for j = 1:length(rxnIndices)
    r = rxnIndices(j);
    capacity = capacities(j);

    if Dir(r) == 1
        model = changeRxnBounds(model,model.rxns{r},0,'l');
        model = changeRxnBounds(model,model.rxns{r},capacity,'u');
    else
        model = changeRxnBounds(model,model.rxns{r},-capacity,'l');
        model = changeRxnBounds(model,model.rxns{r},capacity,'u');
    end
end

end