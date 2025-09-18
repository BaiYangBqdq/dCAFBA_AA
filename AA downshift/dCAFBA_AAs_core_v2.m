%% Date 2024-05-10 this code works for AA downshift using the metabolic model iML515

function [v_R,v_AAs_total,v_C,lambda,v,v_E,v_met,v_arg,v_leu,v_ilv,v_glt,v_trp,v_his,v_ser,v_lys,v_thr,v_phetyr,v_aro,v_cys,v_other,phi_E] = dCAFBA_AAs_core_v2(model,Dir,phi_AAs,phi_Rb,phi_met,phi_arg,phi_leu,phi_ilv,phi_glt,phi_trp,phi_his,phi_ser,phi_lys,phi_thr,phi_phetyr,phi_aro,phi_cys,phi_other,phi_C,w_E,glyc_r,phi_Q,phiE_r,arg_r,aro_r,cys_r,glt_r,his_r,ilv_r,leu_r,lys_r,met_r,other_r,phetyr_r,ser_r,thr_r,trp_r,AAs_r,v_met_f,v_arg_f,v_leu_f,v_ilv_f,v_glt_f,v_trp_f,v_his_f,v_ser_f,v_lys_f,v_thr_f,v_phetyr_f,v_aro_f,v_cys_f,v_other_f,ratios_met,ratios_arg,ratios_leu,ratios_ilv,ratios_glt,ratios_trp,ratios_his,ratios_ser,ratios_lys,ratios_thr,ratios_phetyr,ratios_aro,ratios_cys,ratios_other)

model.protGroup(4).phi0 = phi_Q+phi_C+phi_Rb;
phi_E=1-model.protGroup(4).phi0;

% %met_r
% % v_met_max = ratios_met.*v_AA_f(1);
% % v_met_max = phi_met/w_E;
% %Dir of AAs_r is 1 or 0: Dir==0, reversible; Dir==1,irreversible

for j1=1:length(met_r)
    if Dir(met_r(j1))==1
        model=changeRxnBounds(model,model.rxns{met_r(j1)},0,'l');
        model=changeRxnBounds(model,model.rxns{met_r(j1)},ratios_met.*abs(v_met_f(j1)),'u');%
    else
        model=changeRxnBounds(model,model.rxns{met_r(j1)},-ratios_met.*abs(v_met_f(j1)),'l');
        model=changeRxnBounds(model,model.rxns{met_r(j1)},ratios_met.*abs(v_met_f(j1)),'u');  
    end  
end

% %arg_r
% % v_arg_max = ratios_arg.*v_AA_f(2);
% % v_arg_max = phi_arg/w_E;
for j2=1:length(arg_r)
    if Dir(arg_r(j2))==1
        model=changeRxnBounds(model,model.rxns{arg_r(j2)},0,'l');
        model=changeRxnBounds(model,model.rxns{arg_r(j2)},ratios_arg.*abs(v_arg_f(j2)),'u');
    else
        model=changeRxnBounds(model,model.rxns{arg_r(j2)},-ratios_arg.*abs(v_arg_f(j2)),'l');
        model=changeRxnBounds(model,model.rxns{arg_r(j2)},ratios_arg.*abs(v_arg_f(j2)),'u');  
    end  
end
         
% %leu_r
% % v_leu_max = ratios_leu.*v_AA_f(3);
% % v_leu_max = phi_leu/w_E;
for j3=1:length(leu_r)
    if Dir(leu_r(j3))==1
        model=changeRxnBounds(model,model.rxns{leu_r(j3)},0,'l');
        model=changeRxnBounds(model,model.rxns{leu_r(j3)},ratios_leu.*abs(v_leu_f(j3)),'u');            
    else
        model=changeRxnBounds(model,model.rxns{leu_r(j3)},-ratios_leu.*abs(v_leu_f(j3)),'l');
        model=changeRxnBounds(model,model.rxns{leu_r(j3)},ratios_leu.*abs(v_leu_f(j3)),'u');  
    end
end

% %trp_r
% % v_trp_max = ratios_trp.*v_AA_f(4);
% % v_trp_max = phi_trp/w_E;
for j4=1:length(trp_r)
    if Dir(trp_r(j4))==1
        model=changeRxnBounds(model,model.rxns{trp_r(j4)},0,'l');
        model=changeRxnBounds(model,model.rxns{trp_r(j4)},ratios_trp.*abs(v_trp_f(j4)),'u');
    else
        model=changeRxnBounds(model,model.rxns{trp_r(j4)},-ratios_trp.*abs(v_trp_f(j4)),'l');
        model=changeRxnBounds(model,model.rxns{trp_r(j4)},ratios_trp.*abs(v_trp_f(j4)),'u');  
    end
end
%  
% %glt_r
% % v_glt_max = ratios_glt.*v_AA_f(5);
% % v_glt_max = phi_glt/w_E;
for j5=1:length(glt_r)
    if Dir(glt_r(j5))==1
        model=changeRxnBounds(model,model.rxns{glt_r(j5)},0,'l');
        model=changeRxnBounds(model,model.rxns{glt_r(j5)},ratios_glt.*abs(v_glt_f(j5)),'u');
    else
        model=changeRxnBounds(model,model.rxns{glt_r(j5)},-ratios_glt.*abs(v_glt_f(j5)),'l');
        model=changeRxnBounds(model,model.rxns{glt_r(j5)},ratios_glt.*abs(v_glt_f(j5)),'u');  
    end  
end

% 
% %ilv_r
% % v_ilv_max = ratios_ilv.*v_AA_f(6);
% % v_ilv_max = phi_ilv/w_E;
for j6=1:length(ilv_r)
    if Dir(ilv_r(j6))==1
        model=changeRxnBounds(model,model.rxns{ilv_r(j6)},0,'l');
        model=changeRxnBounds(model,model.rxns{ilv_r(j6)},ratios_ilv.*abs(v_ilv_f(j6)),'u');
    else
        model=changeRxnBounds(model,model.rxns{ilv_r(j6)},-ratios_ilv.*abs(v_ilv_f(j6)),'l');
        model=changeRxnBounds(model,model.rxns{ilv_r(j6)},ratios_ilv.*abs(v_ilv_f(j6)),'u');  
    end  
end 
% 
% %his_r
% % v_his_max = ratios_his.*v_AA_f(7);
% % v_his_max = phi_his/w_E;
for j7=1:length(his_r)
    if Dir(his_r(j7))==1
         model=changeRxnBounds(model,model.rxns{his_r(j7)},0,'l');
         model=changeRxnBounds(model,model.rxns{his_r(j7)},ratios_his.*abs(v_his_f(j7)),'u');
    else
         model=changeRxnBounds(model,model.rxns{his_r(j7)},-ratios_his.*abs(v_his_f(j7)),'l');
         model=changeRxnBounds(model,model.rxns{his_r(j7)},ratios_his.*abs(v_his_f(j7)),'u');  
    end  
end
% 
% %ser_r
% % v_ser_max = ratios_ser.*v_AA_f(8);
% % v_ser_max = phi_ser/w_E;
for j8=1:length(ser_r)
    if Dir(ser_r(j8))==1
        model=changeRxnBounds(model,model.rxns{ser_r(j8)},0,'l');
        model=changeRxnBounds(model,model.rxns{ser_r(j8)},ratios_ser.*abs(v_ser_f(j8)),'u');
    else
        model=changeRxnBounds(model,model.rxns{ser_r(j8)},-ratios_ser.*abs(v_ser_f(j8)),'l');
        model=changeRxnBounds(model,model.rxns{ser_r(j8)},ratios_ser.*abs(v_ser_f(j8)),'u');  
    end  
end
% 
% %lys_r
% % v_lys_max = ratios_lys.*v_AA_f(9);
% % v_lys_max = phi_lys/w_E;
for j9=1:length(lys_r)
    if Dir(lys_r(j9))==1
        model=changeRxnBounds(model,model.rxns{lys_r(j9)},0,'l');
        model=changeRxnBounds(model,model.rxns{lys_r(j9)},ratios_lys.*abs(v_lys_f(j9)),'u');
    else
        model=changeRxnBounds(model,model.rxns{lys_r(j9)},-ratios_lys.*abs(v_lys_f(j9)),'l');
        model=changeRxnBounds(model,model.rxns{lys_r(j9)},ratios_lys.*abs(v_lys_f(j9)),'u');  
    end  
end
% 
% %thr_r
% % v_thr_max = ratios_thr.*v_AA_f(10);
% % v_thr_max = phi_thr/w_E;
for j10=1:length(thr_r)
    if Dir(thr_r(j10))==1
        model=changeRxnBounds(model,model.rxns{thr_r(j10)},0,'l');
        model=changeRxnBounds(model,model.rxns{thr_r(j10)},ratios_thr.*abs(v_thr_f(j10)),'u');
    else
        model=changeRxnBounds(model,model.rxns{thr_r(j10)},-ratios_thr.*abs(v_thr_f(j10)),'l');
        model=changeRxnBounds(model,model.rxns{thr_r(j10)},ratios_thr.*abs(v_thr_f(j10)),'u');  
    end  
end
%  
% %phetyr_r
% % v_phetyr_max = ratios_phetyr.*v_AA_f(11);
% % v_phetyr_max = phi_phetyr/w_E;
for j11=1:length(phetyr_r)
    if Dir(phetyr_r(j11))==1
        model=changeRxnBounds(model,model.rxns{phetyr_r(j11)},0,'l');
        model=changeRxnBounds(model,model.rxns{phetyr_r(j11)},ratios_phetyr.*abs(v_phetyr_f(j11)),'u');
    else
        model=changeRxnBounds(model,model.rxns{phetyr_r(j11)},-ratios_phetyr.*abs(v_phetyr_f(j11)),'l');
        model=changeRxnBounds(model,model.rxns{phetyr_r(j11)},ratios_phetyr.*abs(v_phetyr_f(j11)),'u');  
    end  
end
% 
% %aro_r
% % v_aro_max = ratios_aro.*v_AA_f(12);
% % v_aro_max = phi_aro/w_E;
for j12=1:length(aro_r)
     if Dir(aro_r(j12))==1
         model=changeRxnBounds(model,model.rxns{aro_r(j12)},0,'l');
         model=changeRxnBounds(model,model.rxns{aro_r(j12)},ratios_aro.*abs(v_aro_f(j12)),'u');
     else
         model=changeRxnBounds(model,model.rxns{aro_r(j12)},-ratios_aro.*abs(v_aro_f(j12)),'l');
         model=changeRxnBounds(model,model.rxns{aro_r(j12)},ratios_aro.*abs(v_aro_f(j12)),'u');  
     end  
end
%         
% %cys_r
% % v_cys_max = ratios_cys.*v_AA_f(13);
% % v_cys_max = phi_cys/w_E;
for j13=1:length(cys_r)
    if Dir(cys_r(j13))==1
        model=changeRxnBounds(model,model.rxns{cys_r(j13)},0,'l');
        model=changeRxnBounds(model,model.rxns{cys_r(j13)},ratios_cys.*abs(v_cys_f(j13)),'u');
    else
        model=changeRxnBounds(model,model.rxns{cys_r(j13)},-ratios_cys.*abs(v_cys_f(j13)),'l');
        model=changeRxnBounds(model,model.rxns{cys_r(j13)},ratios_cys.*abs(v_cys_f(j13)),'u');  
    end  
end 
%         
% %other_r
% % v_other_max = ratios_other.*v_AA_f(14);
% % v_other_max = phi_other/w_E;
for j14=1:length(other_r)
    if Dir(other_r(j14))==1
        model=changeRxnBounds(model,model.rxns{other_r(j14)},0,'l');
        model=changeRxnBounds(model,model.rxns{other_r(j14)},ratios_other.*abs(v_other_f(j14)),'u');
    else
        model=changeRxnBounds(model,model.rxns{other_r(j14)},-ratios_other.*abs(v_other_f(j14)),'l');
        model=changeRxnBounds(model,model.rxns{other_r(j14)},ratios_other.*abs(v_other_f(j14)),'u');  
    end  
end 

sol=CAFBA_OptimizeCbModel_glpk(model);
v=sol.x;
v_BOF=sol.f;
lambda=v_BOF;
v_C=abs(sol.x(glyc_r));
% v_R=(0.488+0.281+0.229+0.229+0.087+0.25+0.25+0.582+0.09+0.276+0.428+0.326+0.146+0.176+0.21+0.205+0.241+0.054+0.131+0.402)*v_BOF;
v_R=(0.513689+0.295792+0.241055+0.241055+0.09158+0.26316+0.26316+0.612638+0.094738+0.290529+0.450531+0.343161+0.153686+0.185265+0.221055+0.215792+0.253687+0.056843+0.137896+0.423162)*v_BOF;
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
% v_AAs=v_met+v_arg+v_leu+v_trp+v_glt+v_ilv+v_his+v_ser+v_lys+v_thr+v_aro+v_phetyr+v_cys+v_other;

for m=1:length(AAs_r)
    v_AAs_total(m)=sol.x(AAs_r(m));
end

    
end