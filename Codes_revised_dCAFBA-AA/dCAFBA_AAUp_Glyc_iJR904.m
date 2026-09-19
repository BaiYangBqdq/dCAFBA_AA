%% This script simulates AA upshift using iJR904: Glycerol -> Glycerol+18AAs-serine
clear all
% close all
%% load model and data
load('AAShift_iJR904_glycerol.mat');
load('ExperimentalData.mat');


model = changeRxnBounds(model,'THRAr',0,'l');

%without serine: EX_ser_L_e_
model1=changeRxnBounds(model,{'EX_ala_L_e_','EX_arg_L_e_','EX_asn_L_e_','EX_asp_L_e_','EX_cys_L_e_','EX_glu_L_e_','EX_gln_L_e_','EX_gly_e_','EX_his_L_e_','EX_ile_L_e_','EX_leu_L_e_','EX_lys_L_e_','EX_met_L_e_','EX_phe_L_e_','EX_pro_L_e_','EX_ser_L_e_','EX_thr_L_e_','EX_tyr_L_e_','EX_val_L_e_','EX_trp_L_e_'},0,'b');

model2=changeRxnBounds(model,{'EX_ala_L_e_','EX_arg_L_e_','EX_asn_L_e_','EX_asp_L_e_','EX_glu_L_e_','EX_gln_L_e_','EX_gly_e_','EX_his_L_e_','EX_ile_L_e_','EX_leu_L_e_','EX_lys_L_e_','EX_met_L_e_','EX_phe_L_e_','EX_pro_L_e_','EX_thr_L_e_','EX_val_L_e_','EX_trp_L_e_'},-1000,'l');
model2=changeRxnBounds(model2,{'EX_ala_L_e_','EX_arg_L_e_','EX_asn_L_e_','EX_asp_L_e_','EX_glu_L_e_','EX_gln_L_e_','EX_gly_e_','EX_his_L_e_','EX_ile_L_e_','EX_leu_L_e_','EX_lys_L_e_','EX_met_L_e_','EX_phe_L_e_','EX_pro_L_e_','EX_thr_L_e_','EX_val_L_e_','EX_trp_L_e_'},0,'u');
model2=changeRxnBounds(model2,{'EX_cys_L_e_','EX_tyr_L_e_','EX_ser_L_e_'},0,'b');

glyc_r=find(strcmp(model.rxns,'EX_glyc_e_'));
w_E1 = 1.70*10^-3;
w_E2 = 1.37*10^-3;
model1=setWeights(model1,2,w_E1);%MOPS+Glycerol
model2=setWeights(model2,2,w_E2);%MOPS+Glycerol+18AA-serine
model=assignQualDir(model);% 0, reversible; 1, irreversible-forward
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
other={'ALATA_L','ASNS2','ASNS1','ASPTA', 'GLNS','THRAr','GHMT2','G5SD','GLU5K','P5CR'};%including Alanine &Asparate, Glutmate,Threonine&Lysine, Glycerine&Serine, Argine;
phetyr={'CHORM','PPNDH','PHETA1','LEUTAi','TYRTA'};
ser={'PGCD','PSP_L','PSERT'};
thr={'HSK','THRS'};
trp={'TRPS1','TRPS2','TRPS3','PRAIi','IGPS','ANS','ANPRT'};
%% parameters
dt = 0.005;% unit in hour
t = -3:dt:8;%total time
gamma = 8.35;% per h; obtain from Supplemental Table 6 of Wu'paper in 2023.
lambda_C = 1.17;%strain-specific constant of growth laws;
phi_Rb0 = 0.0445;%obtain from Supplemental Table 6 of Wu'paper in 2023.
lambda_i = 0.68;%experimental data;
lambda_f = 1.17;%pre-shift;1.37
phi_Rb_i=phi_Rb0+lambda_i./gamma;
sigma_i = lambda_i./phi_Rb_i; 
phi_Rb_f=phi_Rb0+lambda_f./gamma;
sigma_f = lambda_f./phi_Rb_f; 
phi_Q = 0.45;% obtain from You et al., 2013
phi_C = 0.03;
phi_AAs_max = 0.218;%from Wu et al., 2023, AA with serine;0.181
alpha_A = 0.16;%h;from Wu et al., 2023 AA with serine

tau_sigma = 10/60;   %10/60,[1,10,30]/60 hours

% k_sigma=3.4;%0.36,3.4,4.8

%Experimental data:pre-shift 18AA-serine
phi_met_f = 1.05*10^-3;
phi_met_i = 32.2*10^-3;
phi_arg_f = 0.20*10^-3;
phi_arg_i = 4.71*10^-3;
phi_leu_f = 0.51*10^-3;
phi_leu_i = 5.85*10^-3;
phi_trp_f = 0.16*10^-3;
phi_trp_i = 2.31*10^-3;
phi_glt_f = 1.12*10^-3;
phi_glt_i = 7.23*10^-3;
phi_ilv_f = 1.8*10^-3;
phi_ilv_i = 13.1*10^-3;
phi_his_f = 1.38*10^-3;
phi_his_i = 4.15*10^-3;
phi_ser_f = 1.25*10^-3;
phi_ser_i = 3.79*10^-3;
phi_lys_f = 2.75*10^-3;
phi_lys_i = 6.32*10^-3;
phi_thr_f = 1.21*10^-3;
phi_thr_i = 4.13*10^-3;
phi_aro_f = 1.67*10^-3;
phi_aro_i = 2.97*10^-3;
phi_phetyr_f = 0.63*10^-3;
phi_phetyr_i = 1.09*10^-3;
phi_cys_f = 6.51*10^-3;
phi_cys_i = 9.88*10^-3;
phi_other_f = 10.5*10^-3;
phi_other_i = 11.5*10^-3;


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
v_bio=nan(size(t));
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
%% Time points selected for CAFVA
FVA_times = [-1, 0, 1, 2 , 4];
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
coreInput.phi_C = phi_C;
coreInput.phi_Q = phi_Q;
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
        phi_Rb(i)=phi_Rb_i;
        phi_AAs(i)=phi_AAs_total_i;
        sigma(i)=sigma_i;
        chi_Rb(i)=phi_Rb_i;
        chi_AAs(i)=phi_AAs_total_i;
        chi_AA_total(:,i)=chi_AA_total(:,1);
        phi_AA_total(:,i)=phi_AA_total(:,1); 
        ratios_total(:,i) = phi_AA_total(:,i) ./ phi_AAs_f(:);
        lambda(i)=lambda_i;
        
        coreState.phi_Rb = phi_Rb(i-1);

        [v_bio(i),v(:,i),phi_E(i),model_t_current] = dCAFBA_AAUp_Glyc_iJR904_pre(model1,coreState,coreInput);    
    else
  
%         if t(i)>= 0 && t(i)< 1/3
%             sigma(i)=sigma_i+k_sigma*t(i);%k=3.6, t=1/3
%         else
%             sigma(i)=sigma_f;
%         end  

       phi_Rb(i)=sigma(i-1)*phi_Rb(i-1)*(chi_Rb(i-1)-phi_Rb(i-1))*dt+phi_Rb(i-1);
       phi_AAs(i)=sigma(i-1)*phi_Rb(i-1)*(chi_AAs(i-1)-phi_AAs(i-1))*dt+phi_AAs(i-1);
       sigma(i) = sigma_f-(sigma_f-sigma_i)*exp(-t(i)/tau_sigma);% Michaelis-Menten-type translational activity
       chi_Rb(i)=phi_Rb0/(1-sigma(i)/gamma);
       chi_AAs(i)=phi_AAs_max-alpha_A*sigma(i)*chi_Rb(i);
       lambda(i)=sigma(i)*phi_Rb(i);  
       coreState.phi_Rb = phi_Rb(i-1);

        [v_bio(i),v(:,i),phi_E(i),model_t_current] = dCAFBA_AAUp_Glyc_iJR904_post(model2,coreState,coreInput);

        %determine which AA is limiting

       ratios = phi_AA_total(:,i-1) ./ phi_AAs_f(:);

       chi_AA_total(:,i) = zeros(length(AAs_group),1);
       eta_total(:,i) = zeros(length(AAs_group),1);

    % Find the Turn On amino acid at the last moment
       onflagt = onflag(:, i-1);%onflag at the last moment     

    % Find next activate candidate AA
       tol = 0.005;
       candidateon = ratios < 1 + tol ;     % Logical flag arr.
       candidateon = candidateon & ~onflagt;
       if any(candidateon)
           flagact(i) = sum(candidateon); 
           onflagt(candidateon) = true;  
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
save('dCAFBA_AA_FVAData.mat','FVAData', 't', '-v7.3');
