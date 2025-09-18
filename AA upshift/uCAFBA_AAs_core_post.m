%date 2025-04-03:this code works for AA upshift:post-shift.
function [v_R,v_C,lambda_FBA,v,v_E,v_met,v_arg,v_leu,v_ilv,v_glt,v_trp,v_his,v_ser,v_lys,v_thr,v_phetyr,v_aro,v_cys,v_other,phi_E] = uCAFBA_AAs_core_post(model,Dir,phi_AAs,phi_Rb,phi_met,phi_arg,phi_leu,phi_ilv,phi_glt,phi_trp,phi_his,phi_ser,phi_lys,phi_thr,phi_phetyr,phi_aro,phi_cys,phi_other,phi_C2,w_C,w_E2,glyc_r,phi_Q,phiE_r,arg_r,aro_r,cys_r,glt_r,his_r,ilv_r,leu_r,lys_r,met_r,other_r,phetyr_r,ser_r,thr_r,trp_r,AAs_r,v_met_i,v_arg_i,v_leu_i,v_ilv_i,v_glt_i,v_trp_i,v_his_i,v_ser_i,v_lys_i,v_thr_i,v_phetyr_i,v_aro_i,v_cys_i,v_other_i,ratios_met,ratios_arg,ratios_leu,ratios_ilv,ratios_glt,ratios_trp,ratios_his,ratios_ser,ratios_lys,ratios_thr,ratios_phetyr,ratios_aro,ratios_cys,ratios_other)

model.protGroup(4).phi0 = phi_Q+phi_C2+phi_Rb;
phi_E=1-model.protGroup(4).phi0;

%met_r
v_met_max = phi_met/w_E2;
for j1=1:length(met_r)
    if Dir(met_r(j1))==1
        model=changeRxnBounds(model,model.rxns{met_r(j1)},0,'l');
        model=changeRxnBounds(model,model.rxns{met_r(j1)},v_met_max,'u');
    else
        model=changeRxnBounds(model,model.rxns{met_r(j1)},-v_met_max,'l');
        model=changeRxnBounds(model,model.rxns{met_r(j1)},v_met_max,'u');  
    end  
end

%arg_r
v_arg_max = phi_arg/w_E2;
for j2=1:length(arg_r)
    if Dir(arg_r(j2))==1
        model=changeRxnBounds(model,model.rxns{arg_r(j2)},0,'l');
        model=changeRxnBounds(model,model.rxns{arg_r(j2)},v_arg_max,'u');
    else
        model=changeRxnBounds(model,model.rxns{arg_r(j2)},-v_arg_max,'l');
        model=changeRxnBounds(model,model.rxns{arg_r(j2)},v_arg_max,'u');  
    end  
end
         
%leu_r
v_leu_max = phi_leu/w_E2;
for j3=1:length(leu_r)
    if Dir(leu_r(j3))==1
        model=changeRxnBounds(model,model.rxns{leu_r(j3)},0,'l');
        model=changeRxnBounds(model,model.rxns{leu_r(j3)},v_leu_max,'u');            
    else
        model=changeRxnBounds(model,model.rxns{leu_r(j3)},-v_leu_max,'l');
        model=changeRxnBounds(model,model.rxns{leu_r(j3)},v_leu_max,'u');  
    end
end
% model=changeRxnBounds(model,'IPPS',1,'u');
% model=changeRxnBounds(model,'IPMD',1,'u');
% model=changeRxnBounds(model,'IPPMIa',-1,'l');
% model=changeRxnBounds(model,'IPPMIb',-1,'l');
% model=changeRxnBounds(model,'IPPS',1,'u');
%ilv_r
v_ilv_max = phi_ilv/w_E2;
for j4=1:length(ilv_r)
    if Dir(ilv_r(j4))==1
        model=changeRxnBounds(model,model.rxns{ilv_r(j4)},0,'l');
        model=changeRxnBounds(model,model.rxns{ilv_r(j4)},v_ilv_max,'u');
    else
        model=changeRxnBounds(model,model.rxns{ilv_r(j4)},-v_ilv_max,'l');
        model=changeRxnBounds(model,model.rxns{ilv_r(j4)},v_ilv_max,'u');  
    end  
end 
 
%glt_r
v_glt_max = phi_glt/w_E2;
for j5=1:length(glt_r)
    if Dir(glt_r(j5))==1
        model=changeRxnBounds(model,model.rxns{glt_r(j5)},0,'l');
        model=changeRxnBounds(model,model.rxns{glt_r(j5)},v_glt_max,'u');
    else
        model=changeRxnBounds(model,model.rxns{glt_r(j5)},-v_glt_max,'l');
        model=changeRxnBounds(model,model.rxns{glt_r(j5)},v_glt_max,'u');  
    end  
end

%trp_r
v_trp_max = phi_trp/w_E2;
for j6=1:length(trp_r)
    if Dir(trp_r(j6))==1
        model=changeRxnBounds(model,model.rxns{trp_r(j6)},0,'l');
        model=changeRxnBounds(model,model.rxns{trp_r(j6)},v_trp_max,'u');
    else
        model=changeRxnBounds(model,model.rxns{trp_r(j6)},-v_trp_max ,'l');
        model=changeRxnBounds(model,model.rxns{trp_r(j6)},v_trp_max,'u');  
    end
end
           
%his_r
v_his_max = phi_his/w_E2;
for j7=1:length(his_r)
    if Dir(his_r(j7))==1
         model=changeRxnBounds(model,model.rxns{his_r(j7)},0,'l');
         model=changeRxnBounds(model,model.rxns{his_r(j7)},v_his_max,'u');
    else
         model=changeRxnBounds(model,model.rxns{his_r(j7)},-v_his_max,'l');
         model=changeRxnBounds(model,model.rxns{his_r(j7)},v_his_max,'u');  
    end  
end

%ser_r
v_ser_max = phi_ser/w_E2;
for j8=1:length(ser_r)
    if Dir(ser_r(j8))==1
        model=changeRxnBounds(model,model.rxns{ser_r(j8)},0,'l');
        model=changeRxnBounds(model,model.rxns{ser_r(j8)},v_ser_max,'u');
    else
        model=changeRxnBounds(model,model.rxns{ser_r(j8)},-v_ser_max,'l');
        model=changeRxnBounds(model,model.rxns{ser_r(j8)},v_ser_max,'u');  
    end  
end

%lys_r
v_lys_max = phi_lys/w_E2;
for j9=1:length(lys_r)
    if Dir(lys_r(j9))==1
        model=changeRxnBounds(model,model.rxns{lys_r(j9)},0,'l');
        model=changeRxnBounds(model,model.rxns{lys_r(j9)},v_lys_max,'u');
    else
        model=changeRxnBounds(model,model.rxns{lys_r(j9)},-v_lys_max,'l');
        model=changeRxnBounds(model,model.rxns{lys_r(j9)},v_lys_max,'u');  
    end  
end

%thr_r
v_thr_max = phi_thr/w_E2;
for j10=1:length(thr_r)
    if Dir(thr_r(j10))==1
        model=changeRxnBounds(model,model.rxns{thr_r(j10)},0,'l');
        model=changeRxnBounds(model,model.rxns{thr_r(j10)},v_thr_max,'u');
    else
        model=changeRxnBounds(model,model.rxns{thr_r(j10)},-v_thr_max,'l');
        model=changeRxnBounds(model,model.rxns{thr_r(j10)},v_thr_max,'u');  
    end  
end
 
%phetyr_r
v_phetyr_max = phi_phetyr/w_E2;
for j11=1:length(phetyr_r)
    if Dir(phetyr_r(j11))==1
        model=changeRxnBounds(model,model.rxns{phetyr_r(j11)},0,'l');
        model=changeRxnBounds(model,model.rxns{phetyr_r(j11)},v_phetyr_max,'u');
    else
        model=changeRxnBounds(model,model.rxns{phetyr_r(j11)},-v_phetyr_max,'l');
        model=changeRxnBounds(model,model.rxns{phetyr_r(j11)},v_phetyr_max,'u');  
    end  
end

%aro_r
v_aro_max = phi_aro/w_E2;
for j12=1:length(aro_r)
     if Dir(aro_r(j12))==1
         model=changeRxnBounds(model,model.rxns{aro_r(j12)},0,'l');
         model=changeRxnBounds(model,model.rxns{aro_r(j12)},v_aro_max,'u');
     else
         model=changeRxnBounds(model,model.rxns{aro_r(j12)},-v_aro_max,'l');
         model=changeRxnBounds(model,model.rxns{aro_r(j12)},v_aro_max,'u');  
     end  
end
        
%cys_r
v_cys_max = phi_cys/w_E2;
for j13=1:length(cys_r)
    if Dir(cys_r(j13))==1
        model=changeRxnBounds(model,model.rxns{cys_r(j13)},0,'l');
        model=changeRxnBounds(model,model.rxns{cys_r(j13)},v_cys_max,'u');
    else
        model=changeRxnBounds(model,model.rxns{cys_r(j13)},-v_cys_max,'l');
        model=changeRxnBounds(model,model.rxns{cys_r(j13)},v_cys_max,'u');  
    end  
end 
        
%other_r
v_other_max = phi_other/w_E2;
for j14=1:length(other_r)
    if Dir(other_r(j14))==1
        model=changeRxnBounds(model,model.rxns{other_r(j14)},0,'l');
        model=changeRxnBounds(model,model.rxns{other_r(j14)},v_other_max,'u');
    else
        model=changeRxnBounds(model,model.rxns{other_r(j14)},-v_other_max,'l');
        model=changeRxnBounds(model,model.rxns{other_r(j14)},v_other_max,'u');  
    end  
end 

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

end
