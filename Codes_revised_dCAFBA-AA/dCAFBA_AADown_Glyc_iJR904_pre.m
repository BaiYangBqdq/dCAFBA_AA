function [lambda,v,phi_E,model_t] = dCAFBA_AADown_Glyc_iJR904_pre(model,coreState,coreInput)

% Compact pre-shift CAFBA solver. Mathematical logic is unchanged.
%
% Required coreState field
%   phi_Rb
%
% Required coreInput fields
%   phi_C
%   phi_Q

model.protGroup(4).phi0 = coreInput.phi_Q + coreInput.phi_C + coreState.phi_Rb;
phi_E=1-model.protGroup(4).phi0;

% Save the complete time-specific model used by CAFBA.
model_t = model;
sol=CAFBA_OptimizeCbModel_glpk(model_t,'constraintSense','U');
v=sol.x;
lambda=sol.f;
end
