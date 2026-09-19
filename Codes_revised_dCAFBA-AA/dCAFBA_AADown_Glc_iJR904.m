%% This script simulates AA downshift using iJR904: Glucose+cAA -> Glucose
clear all
% close all
%% load model and data
load('AAShift_iJR904_glucose.mat');
load('v_AAs_f_iJR904_904.mat');
model = changeRxnBounds(model,'THRAr',0,'l');

model1=changeRxnBounds(model,{'EX_ala_L_e_','EX_arg_L_e_','EX_asn_L_e_','EX_asp_L_e_','EX_cys_L_e_','EX_glu_L_e_','EX_gln_L_e_','EX_gly_e_','EX_his_L_e_','EX_ile_L_e_','EX_leu_L_e_','EX_lys_L_e_','EX_met_L_e_','EX_phe_L_e_','EX_pro_L_e_','EX_ser_L_e_','EX_thr_L_e_','EX_tyr_L_e_','EX_val_L_e_','EX_trp_L_e_'},-1000,'l');
model1=changeRxnBounds(model1,{'EX_ala_L_e_','EX_arg_L_e_','EX_asn_L_e_','EX_asp_L_e_','EX_cys_L_e_','EX_glu_L_e_','EX_gln_L_e_','EX_gly_e_','EX_his_L_e_','EX_ile_L_e_','EX_leu_L_e_','EX_lys_L_e_','EX_met_L_e_','EX_phe_L_e_','EX_pro_L_e_','EX_ser_L_e_','EX_thr_L_e_','EX_tyr_L_e_','EX_val_L_e_','EX_trp_L_e_'},0,'u');

model2=changeRxnBounds(model,{'EX_ala_L_e_','EX_arg_L_e_','EX_asn_L_e_','EX_asp_L_e_','EX_cys_L_e_','EX_glu_L_e_','EX_gln_L_e_','EX_gly_e_','EX_his_L_e_','EX_ile_L_e_','EX_leu_L_e_','EX_lys_L_e_','EX_met_L_e_','EX_phe_L_e_','EX_pro_L_e_','EX_ser_L_e_','EX_thr_L_e_','EX_tyr_L_e_','EX_val_L_e_','EX_trp_L_e_'},0,'b');

glc_r=find(strcmp(model.rxns,'EX_glc_e_'));

w_E1 = 8.69*10^-4;%MOPS+Glucose+cAA;
w_E2 = 9.8*10^-4;%MOPS+Glucose
model1=setWeights(model1,2,w_E1);
model2=setWeights(model2,2,w_E2);
model=assignQualDir(model);% 0, reversible; 1, irreversible-forward, 
Dir=model.qualDir;
phiE_r=model.protGroup(2).rxns;
AAs_group={'met','arg','leu','ilv','glt','trp','his','ser','lys','thr','phetyr','aro','cys','other'};
arg={'ACGS','ACGK','AGPR','ACOTA','ACODA','NACODA','ARGSS','ARGSL','OCBT'};
aro={'PSCVT','DHQS','CHORS','DHQD','SHK3Dr','DDPA','SHKK'};
cys={'ADSK','SADT2','SERAT','PAPSR','SULR','CYSS'};
glt={'GLUDy','GLUSy'};
his={'PRMICIi','IGPDH','HISTP','HSTPT','HISTD','IG3PS','ATPPRT','PRATPP','PRAMPC'};
ilv={'THRD_L','ACHBS','ACLS','KARA1i','KARA2i','DHAD1','DHAD2','ILETA','VALTA','ACHBS','ACLS'};
leu={'IPPS','IPMD','OMCDC','IPPMIa','IPPMIb'};
lys={'ASAD','DHDPS','DHDPRy','THDPS','SDPDS','DAPE','DAPDC'};
met={'HSST','SHSL1','CYSTL','METS','ASPK', 'HSDy'};
other={'ALATA_L','ASNS2','ASNS1','ASPTA', 'GLNS','THRAr','GHMT2','G5SD','GLU5K','P5CR'};
phetyr={'CHORM','PPNDH','PHETA1','LEUTAi','TYRTA'};%(old)
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
gamma = 8.35;% per h; obtain from Supplemental Table 6 of Wu'paper in 2023.
lambda_C = 1.17;%strain-specific constant of growth laws;
phi_Rb0 = 0.0445;%obtain from Supplemental Table 6 of Wu'paper in 2023.
lambda_i = 1.58;%pre-shift glucose+cAA;per h
lambda_f = 0.95;%experimental data;
phi_Rb_i=phi_Rb0+lambda_i./gamma;
sigma_i = lambda_i./phi_Rb_i; 
phi_Rb_f=phi_Rb0+lambda_f./gamma;
sigma_f = lambda_f./phi_Rb_f;
phi_Q = 0.45;% obtain from You et al., 2013
phi_C = 0.03;
phi_AAs_max = 0.322;%from Wu et al., 2023
alpha_A = 0.182;%h;from Wu et al., 2023 

%Experimental data:pre-shift 
phi_met_i = 1.22*10^-3;
phi_met_f = 57*10^-3;
phi_arg_i = 0.559*10^-3;
phi_arg_f = 8.32*10^-3;
phi_leu_i = 0.47*10^-3;
phi_leu_f = 5.74*10^-3;
phi_trp_i = 2.65*10^-3;
phi_trp_f = 6.58*10^-3;
phi_glt_i = 0.94*10^-3;
phi_glt_f = 8.88*10^-3;
phi_ilv_i = 3.22*10^-3;
phi_ilv_f = 15.51*10^-3;
phi_his_i = 0.94*10^-3;
phi_his_f = 3.63*10^-3;
phi_ser_i = 1.0*10^-3;
phi_ser_f = 7.71*10^-3;
phi_lys_i = 2.83*10^-3;
phi_lys_f = 8.69*10^-3;
phi_thr_i = 1.35*10^-3;
phi_thr_f = 5.67*10^-3;
phi_aro_i = 1.98*10^-3;
phi_aro_f = 4.33*10^-3;
phi_phetyr_i = 0.41*10^-3;
phi_phetyr_f = 1.05*10^-3;
phi_cys_i = 1.1*10^-3;
phi_cys_f = 11.48*10^-3;
phi_other_i = 14.9*10^-3;
phi_other_f = 15.37*10^-3;
phi_AAs_total_i = phi_met_i+phi_arg_i+phi_leu_i+phi_ilv_i+phi_glt_i+phi_trp_i+phi_his_i+phi_ser_i+phi_lys_i+phi_thr_i+phi_phetyr_i+phi_aro_i+phi_cys_i+phi_other_i;
phi_AAs_i = [phi_met_i,phi_arg_i,phi_leu_i,phi_ilv_i,phi_glt_i,phi_trp_i,phi_his_i,phi_ser_i,phi_lys_i,phi_thr_i,phi_phetyr_i,phi_aro_i,phi_cys_i,phi_other_i];
phi_AAs_f = [phi_met_f,phi_arg_f,phi_leu_f,phi_ilv_f,phi_glt_f,phi_trp_f,phi_his_f,phi_ser_f,phi_lys_f,phi_thr_f,phi_phetyr_f,phi_aro_f,phi_cys_f,phi_other_f];
%% save data
sigma=nan(size(t));
chi_Rb=nan(size(t));
phi_Rb=nan(size(t));
chi_AAs=nan(size(t));
phi_AAs=nan(size(t));
phi_E=nan(size(t));

chi_AA_total=nan(length(AAs_group),length(t));
phi_AA_total=nan(length(AAs_group),length(t));

lambda=nan(size(t));
v=nan(length(model.rxns),length(t));% all rxn fluxes

ratios_total=nan(length(AAs_group),length(t));
%% initialization
sigma(1) = sigma_i;
lambda(1) = lambda_i;
phi_AAs(1) = phi_AAs_total_i;
phi_Rb(1) = phi_Rb_i;
phi_E(1) = 1-phi_C-phi_Rb(1)-phi_Q;
chi_Rb(1) = phi_Rb_i;
chi_AAs(1) = phi_AAs_total_i;
chi_AA_total(:,1)=phi_AAs_i(:);
phi_AA_total(:,1)=[phi_met_i;phi_arg_i;phi_leu_i;phi_ilv_i;phi_glt_i;phi_trp_i;phi_his_i;phi_ser_i;phi_lys_i;phi_thr_i;phi_phetyr_i;phi_aro_i;phi_cys_i;phi_other_i];
ratios_total(:,1)=[phi_met_i/phi_met_f;phi_arg_i/phi_arg_f;phi_leu_i/phi_leu_f;phi_ilv_i/phi_ilv_f;phi_glt_i/phi_glt_f;phi_trp_i/phi_trp_f;phi_his_i/phi_his_f;phi_ser_i/phi_ser_f;phi_lys_i/phi_lys_f;phi_thr_i/phi_thr_f;phi_phetyr_i/phi_phetyr_f;phi_aro_i/phi_aro_f;phi_cys_i/phi_cys_f;phi_other_i/phi_other_f];
%% Time points selected for CFBA-FVA
FVA_times = [-1, 0, 1, 2, 4];
FVA_indices = zeros(size(FVA_times));

for kFVA = 1:length(FVA_times)
    [~,FVA_indices(kFVA)] = min(abs(t-FVA_times(kFVA)));
end
%actual simulation-grid times
FVA_times_actual = t(FVA_indices);
nFVATimes = length(FVA_indices);

FVAData = repmat(struct( ...
    'requestedTime', [], ...
    'time', [], ...
    'index', [], ...
    'model', [], ...
    'growth', [], ...
    'flux', [], ...
    'phi_Rb', [], ...
    'phi_AAs', [], ...
    'phi_AA_total', [], ...
    'phi_E', [], ...
    'saved', false), nFVATimes, 1);
% Fill information that is already known before simulation
for kFVA = 1:nFVATimes
    FVAData(kFVA).requestedTime = FVA_times(kFVA);
    FVAData(kFVA).time = FVA_times_actual(kFVA);
    FVAData(kFVA).index = FVA_indices(kFVA);
end
%% Compact configuration for CAFBA core functions
% Group order must be identical to AAs_group and ratios_total:
% met, arg, leu, ilv, glt, trp, his, ser, lys, thr, phetyr, aro, cys, other
coreInput.Dir = Dir;
coreInput.phi_C = phi_C;
coreInput.phi_Q = phi_Q;
coreInput.groupNames = AAs_group;
coreInput.groupRxns = { ...
    met_r, arg_r, leu_r, ilv_r, glt_r, trp_r, his_r, ...
    ser_r, lys_r, thr_r, phetyr_r, aro_r, cys_r, other_r};
coreInput.steadyGroupFlux = { ...
    v_met_f, v_arg_f, v_leu_f, v_ilv_f, v_glt_f, v_trp_f, ...
    v_his_f, v_ser_f, v_lys_f, v_thr_f, v_phetyr_f, ...
    v_aro_f, v_cys_f, v_other_f};

%% main loop
tic;  
t0=find(t==0);
q0 = phi_AAs_i(:) ./ phi_AAs_f(:);
[~,idxFirst] = min(q0);
eta_total = zeros( size(phi_AA_total) );
eta_total(idxFirst,1) = 1;
onflag = logical( eta_total );
onflag(idxFirst,:) = true;
flagact = zeros(1, numel(t));

for i=2:length(t)
  
    if t(i)<0
        
        coreState.phi_Rb = phi_Rb(i-1);

        [lambda(i),v(:,i),phi_E(i),model_t_current] = dCAFBA_AADown_Glc_iJR904_pre(model1,coreState,coreInput);

        phi_Rb(i)=sigma(i-1)*phi_Rb(i-1)*(chi_Rb(i-1)-phi_Rb(i-1))*dt+phi_Rb(i-1);
        phi_AAs(i)=sigma(i-1)*phi_Rb(i-1)*(chi_AAs(i-1)-phi_AAs(i-1))*dt+phi_AAs(i-1);
        sigma(i)=lambda(i-1)/phi_Rb(i-1);
        chi_Rb(i)=phi_Rb0/(1-sigma(i)/gamma);
        chi_AAs(i)=phi_AAs_max-alpha_A*sigma(i)*chi_Rb(i);

        chi_AA_total(:,i)=chi_AAs(i)*phi_AAs_i(:)./phi_AAs_total_i;
        phi_AA_total(:,i)=sigma(i-1)*phi_Rb(i-1)*(chi_AA_total(:,i-1)-phi_AA_total(:,i-1))*dt+phi_AA_total(:,i-1);   
        ratios_total(:,i) = phi_AA_total(:,i-1) ./ phi_AAs_f(:);
            
    else
        coreState.phi_Rb = phi_Rb(i-1);
        coreState.ratios = ratios_total(:,i-1);

        [lambda(i),v(:,i),phi_E(i),model_t_current] = dCAFBA_AADown_Glc_iJR904_post(model2,coreState,coreInput);
        phi_Rb(i)=sigma(i-1)*phi_Rb(i-1)*(chi_Rb(i-1)-phi_Rb(i-1))*dt+phi_Rb(i-1);
        phi_AAs(i)=sigma(i-1)*phi_Rb(i-1)*(chi_AAs(i-1)-phi_AAs(i-1))*dt+phi_AAs(i-1);
        sigma(i)=lambda(i-1)/phi_Rb(i-1);
        chi_Rb(i)=phi_Rb0/(1-sigma(i)/gamma);
        chi_AAs(i)=phi_AAs_max-alpha_A*sigma(i)*chi_Rb(i);

    %determine which AA is limiting   

        ratios = phi_AA_total(:,i-1) ./ phi_AAs_f(:);

        chi_AA_total(:,i) = zeros(length(AAs_group),1);
        eta_total(:,i) = zeros(length(AAs_group),1);

    % Find the Turn On amino acid at the last moment
        onflagt = onflag(:, i-1);%onflag at the last moment
        indon = find( onflagt );%index of the current active AA
        threshon = mean( ratios(onflagt) );%mean of the all active AAs at the last moment

    % Find next activate candidate AA
        tol = 0.005;
        candidateon = abs(ratios - threshon) < tol;     
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

%-------------------------------------------------------------
       %Save the exact model and state used at selected FVA time points
%-------------------------------------------------------------
       
    savePositions = find(FVA_indices == i);
    for kSave = 1:length(savePositions)
        dataIndex = savePositions(kSave);
        FVAData(dataIndex).model = model_t_current;
        FVAData(dataIndex).growth = lambda(i);
        FVAData(dataIndex).flux = v(:, i);
        FVAData(dataIndex).phi_Rb = phi_Rb(i);
        FVAData(dataIndex).phi_AAs = phi_AAs(i);
        FVAData(dataIndex).phi_AA_total = phi_AA_total(:, i);
        FVAData(dataIndex).phi_E = phi_E(i);
        FVAData(dataIndex).saved = true;
    end
    
 end
    
T=toc;
% stop;
save('dCAFBA_AA_FVAData.mat','FVAData','t','-v7.3');

