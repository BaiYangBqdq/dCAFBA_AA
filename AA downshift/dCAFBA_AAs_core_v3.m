%% Date 2024-05-15 this code works for AA downshift using metabolic model iJR904;

function [v_R,v_AAs_total,v_C,lambda,v,v_E,v_met,v_arg,v_leu,v_ilv,v_glt,v_trp,v_his,v_ser,v_lys,v_thr,v_phetyr,v_aro,v_cys,v_other,phi_E] = dCAFBA_AAs_core_v3(model,Dir,phi_AAs,phi_Rb,phi_met,phi_arg,phi_leu,phi_ilv,phi_glt,phi_trp,phi_his,phi_ser,phi_lys,phi_thr,phi_phetyr,phi_aro,phi_cys,phi_other,phi_C1,w_E1,glyc_r,phi_Q,phiE_r,arg_r,aro_r,cys_r,glt_r,his_r,ilv_r,leu_r,lys_r,met_r,other_r,phetyr_r,ser_r,thr_r,trp_r,AAs_r)

model.protGroup(4).phi0 = phi_Q+phi_C1+phi_Rb;
phi_E=1-model.protGroup(4).phi0;
sol=CAFBA_OptimizeCbModel_glpk(model);
v=sol.x;
v_BOF=sol.f;
lambda=v_BOF;
v_C=abs(sol.x(glyc_r));
v_R=(0.488+0.281+0.229+0.229+0.087+0.25+0.25+0.582+0.09+0.276+0.428+0.326+0.146+0.176+0.21+0.205+0.241+0.054+0.131+0.402)*v_BOF;
v_E=sum(abs(sol.x(phiE_r)));
v_met=sum(abs(sol.x(met_r)));
v_arg=sum(abs(sol.x(arg_r)));
v_leu=sum(abs(sol.x(leu_r)));
v_ilv=sum(abs(sol.x(ilv_r)));
v_glt=sum(abs(sol.x(glt_r)));
v_trp=sum(abs(sol.x(trp_r)));
v_his=sum(abs(sol.x(his_r)));
v_ser=sum(abs(sol.x(ser_r)));
v_lys=sum(abs(sol.x(lys_r)));
v_thr=sum(abs(sol.x(thr_r)));
v_phetyr=sum(abs(sol.x(phetyr_r)));
v_aro=sum(abs(sol.x(aro_r)));
v_cys=sum(abs(sol.x(cys_r)));
v_other=sum(abs(sol.x(other_r)));

for m=1:length(AAs_r)
    v_AAs_total(m)=sol.x(AAs_r(m));
end 

end
