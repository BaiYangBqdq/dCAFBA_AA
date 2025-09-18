%% Date 2024-05-10 this code works for AA downshift using the metabolic model iML515;

%% Glycerol+18AAs+/-ser -> Glycerol

clear all
% close all
%% load model and data
load ('CAFBA_iML1515.mat');
load('Experiments.mat');
load('phi_data.mat');
load('v_AAs_f_iML1515.mat');

model=addProteinGroupsToModel_iML1515(model,'glyc');%C is glycerol
model=changeRxnBounds(model,{'HPYRRx','TRSARr','HPYRI'},0,'b');%futile cycle 
model1=changeRxnBounds(model,{'EX_ala__L_e','EX_arg__L_e','EX_asn__L_e','EX_asp__L_e','EX_cys__L_e','EX_glu__L_e','EX_gln__L_e','EX_gly_e','EX_his__L_e','EX_ile__L_e','EX_leu__L_e','EX_lys__L_e','EX_met__L_e','EX_phe__L_e','EX_pro__L_e','EX_ser__L_e','EX_thr__L_e','EX_tyr__L_e','EX_val__L_e','EX_trp__L_e'},-1000,'l');
model1=changeRxnBounds(model1,{'EX_ala__L_e','EX_arg__L_e','EX_asn__L_e','EX_asp__L_e','EX_cys__L_e','EX_glu__L_e','EX_gln__L_e','EX_gly_e','EX_his__L_e','EX_ile__L_e','EX_leu__L_e','EX_lys__L_e','EX_met__L_e','EX_phe__L_e','EX_pro__L_e','EX_ser__L_e','EX_thr__L_e','EX_tyr__L_e','EX_val__L_e','EX_trp__L_e'},0,'u');
model2=changeRxnBounds(model,{'EX_ala__L_e','EX_arg__L_e','EX_asn__L_e','EX_asp__L_e','EX_cys__L_e','EX_glu__L_e','EX_gln__L_e','EX_gly_e','EX_his__L_e','EX_ile__L_e','EX_leu__L_e','EX_lys__L_e','EX_met__L_e','EX_phe__L_e','EX_pro__L_e','EX_ser__L_e','EX_thr__L_e','EX_tyr__L_e','EX_val__L_e','EX_trp__L_e'},0,'b');
glyc_r=find(strcmp(model.rxns,'EX_glyc_e'));
w_E1 = 4.45*10^-4;%18AAs-ser
w_E = 8.3*10^-4;
model1=setWeights(model1,2,w_E1);
model2=setWeights(model2,2,w_E);
model=assignQualDir(model);%Assigns a qualitative direction to each reaction based on the upper and lower bounds; 0, reversible, 1, irreversible-forward, 
Dir=model.qualDir;
phiE_r=model.protGroup(2).rxns;
AAs_group={'met','arg','leu','ilv','glt','trp','his','ser','lys','thr','phetyr','aro','cys','other'};
arg={'ACGS','ACGK','AGPR','ACOTA','ACODA','ARGSS','ARGSL','OCBT'};
aro={'PSCVT','DHQS','CHORS','DHQTi','SHK3Dr','DDPA','SHKK'};
cys={'ADSK','SADT2','SERAT','PAPSR','SULR','CYSS'};
glt={'GLUDy','GLUSy'};
his={'PRMICI','IGPDH','HISTP','HSTPT','HISTD','IG3PS','ATPPRT','PRATPP','PRAMPC'};
ilv={'THRD_L','ACHBS','ACLS','KARA1','KARA2','DHAD1','DHAD2','ILETA','VALTA','ACHBS','ACLS'};
leu={'IPPS','IPMD','OMCDC','IPPMIa','IPPMIb'};
lys={'ASAD','DHDPS','DHDPRy','THDPS','SDPDS','DAPE','DAPDC'};
met={'HSST','SHSL1','CYSTL','METS','ASPK','HSDy'};
other={'ALATA_L','ASNS2','ASNS1','ASPTA','GLNS','THRA2','GHMT2r','G5SD','GLU5K','P5CR'};
phetyr={'CHORM','PPNDH','PHETA1','LEUTAi','TYRTA'};
ser={'PGCD','PSP_L','PSERT'};
thr={'HSK','THRS'};
trp={'TRPS1','TRPS2','TRPS3','PRAIi','IGPS','ANS','ANPRT'};
AAs=[met,arg,leu,ilv,glt,trp,his,ser,lys,thr,phetyr,aro,cys,other];
%extract IDs of each AA
for i=1:length(AAs)
    AAs_r(i)=find(strcmp(model.rxns,AAs(i)));
end
for i=1:length(arg)
    arg_r(i)=find(strcmp(model.rxns,arg(i)));
end
for i=1:length(aro)
    aro_r(i)=find(strcmp(model.rxns,aro(i)));
end
for i=1:length(cys)
    cys_r(i)=find(strcmp(model.rxns,cys(i)));
end
for i=1:length(glt)
    glt_r(i)=find(strcmp(model.rxns,glt(i)));
end 
for i=1:length(his)
    his_r(i)=find(strcmp(model.rxns,his(i)));
end
for i=1:length(ilv)
    ilv_r(i)=find(strcmp(model.rxns,ilv(i)));
end
for i=1:length(leu)
    leu_r(i)=find(strcmp(model.rxns,leu(i)));
end
for i=1:length(lys)
    lys_r(i)=find(strcmp(model.rxns,lys(i)));
end
for i=1:length(met)
    met_r(i)=find(strcmp(model.rxns,met(i)));
end
for i=1:length(other)
    other_r(i)=find(strcmp(model.rxns,other(i)));
end
for i=1:length(phetyr)
    phetyr_r(i)=find(strcmp(model.rxns,phetyr(i)));
end
for i=1:length(ser)
    ser_r(i)=find(strcmp(model.rxns,ser(i)));
end
for i=1:length(thr)
    thr_r(i)=find(strcmp(model.rxns,thr(i)));
end
for i=1:length(trp)
    trp_r(i)=find(strcmp(model.rxns,trp(i)));
end
%% parameters
dt = 0.005;% unit in hour
t = -2:dt:6;%total time
gamma = 8.35;
lambda_C = 1.17;
phi_Rb0 = 0.0445;
% lambda_i = 1.40;%pre-shift 18AA;per h
lambda_i = 1.17;%pre-shift 18AA-minus-ser;
lambda_f = 0.68;%experimental data;
phi_Rb_i=phi_Rb0+lambda_i./gamma;%growth laws
sigma_i = lambda_i./phi_Rb_i; 

phi_Q = 0.35;
phi_C = 0.0;
phi_AAs_max = 0.218;%from Wu et al., 2023, AA without serine
alpha_A = 0.16;%h;from Wu et al., 2023,AA without serine

%Experimental data:pre-shift 18AA-minus-ser
phi_met_i = 1.05*10^-3;
phi_met_f = 32.2*10^-3;
phi_arg_i = 0.20*10^-3;
phi_arg_f = 4.71*10^-3;
phi_leu_i = 0.51*10^-3;
phi_leu_f = 5.85*10^-3;
phi_trp_i = 0.16*10^-3;
phi_trp_f = 2.31*10^-3;
phi_glt_i = 1.12*10^-3;
phi_glt_f = 7.23*10^-3;
phi_ilv_i = 1.8*10^-3;
phi_ilv_f = 13.1*10^-3;
phi_his_i = 1.38*10^-3;
phi_his_f = 4.15*10^-3;
phi_ser_i = 1.25*10^-3;
phi_ser_f = 3.79*10^-3;
phi_lys_i = 2.75*10^-3;
phi_lys_f = 6.32*10^-3;
phi_thr_i = 1.21*10^-3;
phi_thr_f = 4.13*10^-3;
phi_aro_i = 1.67*10^-3;
phi_aro_f = 2.97*10^-3;
phi_phetyr_i = 0.63*10^-3;
phi_phetyr_f = 1.09*10^-3;
phi_cys_i = 6.51*10^-3;
phi_cys_f = 9.88*10^-3;
phi_other_i = 10.5*10^-3;
phi_other_f = 11.5*10^-3;
phi_AAs_total_i = phi_met_i+phi_arg_i+phi_leu_i+phi_ilv_i+phi_glt_i+phi_trp_i+phi_his_i+phi_ser_i+phi_lys_i+phi_thr_i+phi_phetyr_i+phi_aro_i+phi_cys_i+phi_other_i;
phi_AAs_i = [phi_met_i,phi_arg_i,phi_leu_i,phi_ilv_i,phi_glt_i,phi_trp_i,phi_his_i,phi_ser_i,phi_lys_i,phi_thr_i,phi_phetyr_i,phi_aro_i,phi_cys_i,phi_other_i];
phi_AAs_f = [phi_met_f,phi_arg_f,phi_leu_f,phi_ilv_f,phi_glt_f,phi_trp_f,phi_his_f,phi_ser_f,phi_lys_f,phi_thr_f,phi_phetyr_f,phi_aro_f,phi_cys_f,phi_other_f];
%% save data
sigma=nan(size(t));
chi_Rb=nan(size(t));
phi_Rb=nan(size(t));
chi_AAs=nan(size(t));
phi_AAs=nan(size(t));
phi_met=nan(size(t));
phi_arg=nan(size(t));
phi_leu=nan(size(t));
phi_trp=nan(size(t));
phi_glt=nan(size(t));
phi_ilv=nan(size(t));
phi_his=nan(size(t));
phi_ser=nan(size(t));
phi_lys=nan(size(t));
phi_thr=nan(size(t));
phi_aro=nan(size(t));
phi_phetyr=nan(size(t));
phi_cys=nan(size(t));
phi_other=nan(size(t));
chi_met=nan(size(t));
chi_arg=nan(size(t));
chi_leu=nan(size(t));
chi_trp=nan(size(t));
chi_glt=nan(size(t));
chi_ilv=nan(size(t));
chi_his=nan(size(t));
chi_ser=nan(size(t));
chi_lys=nan(size(t));
chi_thr=nan(size(t));
chi_aro=nan(size(t));
chi_phetyr=nan(size(t));
chi_cys=nan(size(t));
chi_other=nan(size(t));
phi_E=nan(size(t));
lambda=nan(size(t));
v_R=nan(size(t));%protein synthesis flux
v_C=nan(size(t));%C uptake flux
v_E=nan(size(t));%C uptake flux
v=nan(length(model.rxns),length(t));% all rxn fluxes
chi_AA_total=nan(length(AAs_group),length(t));
phi_AA_total=nan(length(AAs_group),length(t));
v_AAGroup_total=nan(length(AAs_group),length(t));
v_AA_max_total=nan(length(AAs_group),length(t)-1);
ratios_total=nan(length(AAs_group),length(t));
v_AAs_total=nan(length(AAs_r),length(t));%flux of all AA reactions

%% initialization
sigma(1) = sigma_i;
lambda(1) = lambda_i;
phi_AAs(1) = phi_AAs_total_i;
phi_Rb(1) = phi_Rb_i;
phi_E(1) = 1-phi_C-phi_Rb(1)-phi_Q;
chi_Rb(1) = phi_Rb_i;
chi_AAs(1) = phi_AAs_total_i;
chi_AA_total(:,1)=[phi_met_i;0;0;0;0;0;0;0;0;0;0;0;0;0];
phi_AA_total(:,1)=[phi_met_i;phi_arg_i;phi_leu_i;phi_ilv_i;phi_glt_i;phi_trp_i;phi_his_i;phi_ser_i;phi_lys_i;phi_thr_i;phi_phetyr_i;phi_aro_i;phi_cys_i;phi_other_i];

v_R(1) = 4.6;
v_C(1) = 10;
v_met(1)=2.6563;
v_arg(1)=2.8845;
v_leu(1)=2.5113;
v_trp(1)=0.3168;
v_glt(1)=9.7671;
v_ilv(1)=6.3114;
v_his(1)=0.9505;
v_ser(1)=5.6825;
v_lys(1)=3.6503;
v_thr(1)=1.2134;
v_aro(1)=3.3762;
v_phetyr(1)=1.4293;
v_cys(1)=1.641;
v_other(1)=11.9547;
v_AAGroup_total(:,1)=[v_met(1);v_arg(1);v_leu(1);v_ilv(1);v_glt(1);v_trp(1);v_his(1);v_ser(1);v_lys(1);v_thr(1);v_phetyr(1);v_aro(1);v_cys(1);v_other(1)];
ratios_total(:,1)=[phi_met_i/phi_met_f;phi_arg_i/phi_arg_f;phi_leu_i/phi_leu_f;phi_ilv_i/phi_ilv_f;phi_glt_i/phi_glt_f;phi_trp_i/phi_trp_f;phi_his_i/phi_his_f;phi_ser_i/phi_ser_f;phi_lys_i/phi_lys_f;phi_thr_i/phi_thr_f;phi_phetyr_i/phi_phetyr_f;phi_aro_i/phi_aro_f;phi_cys_i/phi_cys_f;phi_other_i/phi_other_f];

v_met_total(:,1)=[0.1;0.1;0.1;0.1;0.1;0.1];
v_arg_total(:,1)=[0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1];
% v_leu_total(:,1)=[0.1;0.1;0.1;0.1];
v_leu_total(:,1)=[0.1;0.1;0.1;0.1;0.1];
v_ilv_total(:,1)=[0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1];
v_glt_total(:,1)=[0.1;0.1];
v_trp_total(:,1)=[0.1;0.1;0.1;0.1;0.1;0.1;0.1];
v_his_total(:,1)=[0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1];
v_ser_total(:,1)=[0.1;0.1;0.1];
v_lys_total(:,1)=[0.1;0.1;0.1;0.1;0.1;0.1;0.1];
v_thr_total(:,1)=[0.1;0.1];
v_phetyr_total(:,1)=[0.1;0.1;0.1;0.1;0.1];
v_aro_total(:,1)=[0.1;0.1;0.1;0.1;0.1;0.1;0.1];
v_cys_total(:,1)=[0.1;0.1;0.1;0.1;0.1;0.1];
v_other_total(:,1)=[0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1];

%% main loop
tic;  
t0=find(t==0);
eta_total = zeros( size(phi_AA_total) );
eta_total(1, 1) = 1;
onflag = logical( eta_total );
onflag(1, :) = true;
flagact = zeros(1, numel(t));

for i=2:length(t)
    if t(i)<0
        [v_R(i),v_AAs_total(:,i),v_C(i),lambda(i),v(:,i),v_E(i),v_AAGroup_total(1,i),v_AAGroup_total(2,i),v_AAGroup_total(3,i),v_AAGroup_total(4,i),v_AAGroup_total(5,i),v_AAGroup_total(6,i),v_AAGroup_total(7,i),v_AAGroup_total(8,i),v_AAGroup_total(9,i),v_AAGroup_total(10,i),v_AAGroup_total(11,i),v_AAGroup_total(12,i),v_AAGroup_total(13,i),v_AAGroup_total(14,i),phi_E(i)] = dCAFBA_AAs_core_v1(model1,Dir,phi_AAs(i-1),phi_Rb(i-1),phi_AA_total(1,i-1),phi_AA_total(2,i-1),phi_AA_total(3,i-1),phi_AA_total(4,i-1),phi_AA_total(5,i-1),phi_AA_total(6,i-1),phi_AA_total(7,i-1),phi_AA_total(8,i-1),phi_AA_total(9,i-1),phi_AA_total(10,i-1),phi_AA_total(11,i-1),phi_AA_total(12,i-1),phi_AA_total(13,i-1),phi_AA_total(14,i-1),phi_C,w_E1,glyc_r,phi_Q,phiE_r,arg_r,aro_r,cys_r,glt_r,his_r,ilv_r,leu_r,lys_r,met_r,other_r,phetyr_r,ser_r,thr_r,trp_r,AAs_r);
        phi_Rb(i)=sigma(i-1)*phi_Rb(i-1)*(chi_Rb(i-1)-phi_Rb(i-1))*dt+phi_Rb(i-1);
        phi_AAs(i)=sigma(i-1)*phi_Rb(i-1)*(chi_AAs(i-1)-phi_AAs(i-1))*dt+phi_AAs(i-1);
        sigma(i)=lambda(i-1)/phi_Rb(i-1);
        chi_Rb(i)=phi_Rb0/(1-sigma(i)/gamma);
        chi_AAs(i)=phi_AAs_max-alpha_A*sigma(i)*chi_Rb(i);
        chi_AA_total(:,i)=chi_AAs(i)*phi_AAs_i(:)./phi_AAs_total_i;
        phi_AA_total(:,i)=sigma(i-1)*phi_Rb(i-1)*(chi_AA_total(:,i-1)-phi_AA_total(:,i-1))*dt+phi_AA_total(:,i-1);   
        ratios_total(:,i) = phi_AA_total(:,i-1) ./ phi_AAs_f(:);
    else
        [v_R(i),v_AAs_total(:,i),v_C(i),lambda(i),v(:,i),v_E(i),v_AAGroup_total(1,i),v_AAGroup_total(2,i),v_AAGroup_total(3,i),v_AAGroup_total(4,i),v_AAGroup_total(5,i),v_AAGroup_total(6,i),v_AAGroup_total(7,i),v_AAGroup_total(8,i),v_AAGroup_total(9,i),v_AAGroup_total(10,i),v_AAGroup_total(11,i),v_AAGroup_total(12,i),v_AAGroup_total(13,i),v_AAGroup_total(14,i),phi_E(i)] = dCAFBA_AAs_core_v2(model2,Dir,phi_AAs(i-1),phi_Rb(i-1),phi_AA_total(1,i-1),phi_AA_total(2,i-1),phi_AA_total(3,i-1),phi_AA_total(4,i-1),phi_AA_total(5,i-1),phi_AA_total(6,i-1),phi_AA_total(7,i-1),phi_AA_total(8,i-1),phi_AA_total(9,i-1),phi_AA_total(10,i-1),phi_AA_total(11,i-1),phi_AA_total(12,i-1),phi_AA_total(13,i-1),phi_AA_total(14,i-1),phi_C,w_E,glyc_r,phi_Q,phiE_r,arg_r,aro_r,cys_r,glt_r,his_r,ilv_r,leu_r,lys_r,met_r,other_r,phetyr_r,ser_r,thr_r,trp_r,AAs_r,v_met_f,v_arg_f,v_leu_f,v_ilv_f,v_glt_f,v_trp_f,v_his_f,v_ser_f,v_lys_f,v_thr_f,v_phetyr_f,v_aro_f,v_cys_f,v_other_f,ratios_total(1,i-1),ratios_total(2,i-1),ratios_total(3,i-1),ratios_total(4,i-1),ratios_total(5,i-1),ratios_total(6,i-1),ratios_total(7,i-1),ratios_total(8,i-1),ratios_total(9,i-1),ratios_total(10,i-1),ratios_total(11,i-1),ratios_total(12,i-1),ratios_total(13,i-1),ratios_total(14,i-1));
        phi_Rb(i)=sigma(i-1)*phi_Rb(i-1)*(chi_Rb(i-1)-phi_Rb(i-1))*dt+phi_Rb(i-1);
        phi_AAs(i)=sigma(i-1)*phi_Rb(i-1)*(chi_AAs(i-1)-phi_AAs(i-1))*dt+phi_AAs(i-1);
        sigma(i)=lambda(i-1)/phi_Rb(i-1);
        chi_Rb(i)=phi_Rb0/(1-sigma(i)/gamma);
        chi_AAs(i)=phi_AAs_max-alpha_A*sigma(i)*chi_Rb(i);

    %determine which AA is limiting   

        ratios = phi_AA_total(:,i-1) ./ phi_AAs_f(:);

        chi_AA_total(:,i) = zeros(length(AAs_group),1);
        eta_total(:,i) = zeros(length(AAs_group),1);
%         eta_total(ind_min,i)=1;

    % Find the Turn On amino acid at the last moment
       onflagt = onflag(:, i-1);%onflag at the last moment
       indon = find( onflagt );%index of the current active AA
       threshon = mean( ratios(onflagt) );%mean of the all active AAs at the last moment

    % Find next activate candidate AA
       tol = 0.005;
       candidateon = abs(ratios - threshon) < tol;     % Logical flag arr.
       candidateon = candidateon & ~onflagt;
       minr = min( ratios(candidateon) );%ratio of the newly activate AA
       if ~isempty( minr )
           flagact(i) = 1;
           indminr = find( ratios==minr );
  
        % Update onflag
           onflagt( indminr ) = true;
       end
    
    % Update eata_totall
       eta_total( onflagt, i ) = phi_AAs_f(onflagt) ./ sum( phi_AAs_f( onflagt ) );
       chi_AA_total( onflagt, i )=chi_AAs(i)*eta_total( onflagt, i );  

       onflag(:, i) = onflagt;
    
       phi_AA_total(:,i)=sigma(i-1)*phi_Rb(i-1)*(chi_AA_total(:,i-1)-phi_AA_total(:,i-1))*dt+phi_AA_total(:,i-1); 
       ratios_total(:,i) = ratios;
    end
    
end 
T=toc;
% stop;
