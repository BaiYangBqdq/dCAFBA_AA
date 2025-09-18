%date 2025-04-03:this code works for AA upshift:pre-shift.

function [v_R,v_C,lambda_FBA,v,v_E,v_met,v_arg,v_leu,v_ilv,v_glt,v_trp,v_his,v_ser,v_lys,v_thr,v_phetyr,v_aro,v_cys,v_other,phi_E] = uCAFBA_AAs_core_pre(model,Dir,phi_AAs,phi_Rb,phi_met,phi_arg,phi_leu,phi_ilv,phi_glt,phi_trp,phi_his,phi_ser,phi_lys,phi_thr,phi_phetyr,phi_aro,phi_cys,phi_other,phi_C1,w_C,w_E,glyc_r,phi_Q,phiE_r,arg_r,aro_r,cys_r,glt_r,his_r,ilv_r,leu_r,lys_r,met_r,other_r,phetyr_r,ser_r,thr_r,trp_r,AAs_r) 

model.protGroup(4).phi0 = phi_Q+phi_C1+phi_Rb;
phi_E=1-model.protGroup(4).phi0;

sol=CAFBA_OptimizeCbModel_glpk(model);
v=sol.x;
v_BOF=sol.f;
lambda_FBA=v_BOF;
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

% for m=1:length(AAs_r)
%     v_AAs_total(m)=sol.x(AAs_r(m));
% end 

% for m1=1:length(met_r)
%     v_met_total(m1)=sol.x(met_r(m1));
% end 
% 
% for m2=1:length(arg_r)
%     v_arg_total(m2)=sol.x(arg_r(m2));
% end 
% 
% for m3=1:length(leu_r)
%     v_leu_total(m3)=sol.x(leu_r(m3));
% end 
% 
% for m4=1:length(ilv_r)
%     v_ilv_total(m4)=sol.x(ilv_r(m4));
% end 
% 
% for m5=1:length(glt_r)
%     v_glt_total(m5)=sol.x(glt_r(m5));
% end 
% 
% for m6=1:length(trp_r)
%     v_trp_total(m6)=sol.x(trp_r(m6));
% end 
% 
% for m7=1:length(his_r)
%     v_his_total(m7)=sol.x(his_r(m7));
% end 
% 
% for m8=1:length(ser_r)
%     v_ser_total(m8)=sol.x(ser_r(m8));
% end 
% 
% for m9=1:length(lys_r)
%     v_lys_total(m9)=sol.x(lys_r(m9));
% end 
% 
% for m10=1:length(thr_r)
%     v_thr_total(m10)=sol.x(thr_r(m10));
% end 
% 
% for m11=1:length(phetyr_r)
%     v_phetyr_total(m11)=sol.x(phetyr_r(m11));
% end 
% 
% for m12=1:length(aro_r)
%     v_aro_total(m12)=sol.x(aro_r(m12));
% end 
% 
% for m13=1:length(cys_r)
%     v_cys_total(m13)=sol.x(cys_r(m13));
% end 
% 
% for m14=1:length(other_r)
%     v_other_total(m14)=sol.x(other_r(m14));
% end

end