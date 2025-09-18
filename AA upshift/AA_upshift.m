%% Date 2025-09: Dynamic simulation of amino acid upshift (Glycerol -> Glycerol + 18 AAs)
%% assume sigma dynamics is a step function; and add adaptation time to sigma
% clear all
% close all
%% load model and data
% fprintf('Loading metabolic model and experimental data...\n')
load('model.mat.mat');
load('Experiments.mat');
load('v_AAs_i.mat');

% Define the initial and final metabolic models
% Model for pre-shift (glyc is the sole carbon source)
model1=changeRxnBounds(model,{'EX_ala_L_e_','EX_arg_L_e_','EX_asn_L_e_','EX_asp_L_e_','EX_cys_L_e_','EX_glu_L_e_','EX_gln_L_e_','EX_gly_e_','EX_his_L_e_','EX_ile_L_e_','EX_leu_L_e_','EX_lys_L_e_','EX_met_L_e_','EX_phe_L_e_','EX_pro_L_e_','EX_ser_L_e_','EX_thr_L_e_','EX_tyr_L_e_','EX_val_L_e_','EX_trp_L_e_'},0,'b');
% Model for post-shift (glyc and 18 AAs are available)
model2=changeRxnBounds(model,{'EX_ala_L_e_','EX_arg_L_e_','EX_asn_L_e_','EX_asp_L_e_','EX_cys_L_e_','EX_glu_L_e_','EX_gln_L_e_','EX_gly_e_','EX_his_L_e_','EX_ile_L_e_','EX_leu_L_e_','EX_lys_L_e_','EX_met_L_e_','EX_phe_L_e_','EX_pro_L_e_','EX_ser_L_e_','EX_thr_L_e_','EX_tyr_L_e_','EX_val_L_e_','EX_trp_L_e_'},-1000,'l');
model2=changeRxnBounds(model2,{'EX_ala_L_e_','EX_arg_L_e_','EX_asn_L_e_','EX_asp_L_e_','EX_cys_L_e_','EX_glu_L_e_','EX_gln_L_e_','EX_gly_e_','EX_his_L_e_','EX_ile_L_e_','EX_leu_L_e_','EX_lys_L_e_','EX_met_L_e_','EX_phe_L_e_','EX_pro_L_e_','EX_ser_L_e_','EX_thr_L_e_','EX_tyr_L_e_','EX_val_L_e_','EX_trp_L_e_'},0,'u');

% Get reaction indices for key reactions
glyc_r=find(strcmp(model.rxns,'EX_glyc_e_'));
phiE_r=model.protGroup(2).rxns;

%Assigns a qualitative direction to each reaction based on the upper and lower bounds; 0, reversible, 1, irreversible-forward,
model=assignQualDir(model); 
Dir=model.qualDir;
% Set weight for pre and post shift
w_C = 3.07*10^-3;
w_E1 = 1.7*10^-3;
w_E2 = 2.5*10^-4;
model1=setWeights(model1,2,w_E1);
model2=setWeights(model2,2,w_E2);

% Define AA groups and their corresponding reaction indices
AAs_group={'met','arg','leu','ilv','glt','trp','his','ser','lys','thr','phetyr','aro','cys','other'};
arg={'ACGS','ACGK','AGPR','ACOTA','ACODA','NACODA','ARGSS','ARGSL','OCBT'};
aro={'PSCVT','DHQS','CHORS','DHQD','SHK3Dr','DDPA','SHKK'};
cys={'ADSK','SADT2','SERAT','PAPSR','SULR','CYSS'};
glt={'GLUDy','GLUSy'};
his={'PRMICIi','IGPDH','HISTP','HSTPT','HISTD','IG3PS','ATPPRT','PRATPP','PRAMPC'};
% ilv={'THRD_L','ACHBS','ACLS','KARA1i','KARA2i','DHAD1','DHAD2','PHETA1','LEUTAi','ILETA','VALTA','ACHBS','ACLS'};
ilv={'THRD_L','ACHBS','ACLS','KARA1i','KARA2i','DHAD1','DHAD2','ILETA','VALTA','ACHBS','ACLS'};
% leu={'IPPS','IPMD','IPPMIa','IPPMIb'};
leu={'IPPS','IPMD','OMCDC','IPPMIa','IPPMIb'};
lys={'ASAD','DHDPS','DHDPRy','THDPS','SDPDS','DAPE','DAPDC'};
% lys={'ASAD','DHDPS','DHDPRy','THDPS','SDPDS','DAPE','DAPDC','ASPK'};
met={'HSST','SHSL1','CYSTL','METS','ASPK', 'HSDy'};
other={'ALATA_L','ASNS2','ASNS1','ASPTA', 'GLNS','THRAr','GHMT2','G5SD','GLU5K','P5CR'};%including Alanine &Asparate, Glutmate,Threonine&Lysine, Glycerine&Serine, Argine
% other={'ALATA_L','ASNS2','ASNS1','TYRTA','ASPTA', 'PHETA1','GLNS','THRAr','GHMT2','G5SD','GLU5K','P5CR'};
phetyr={'CHORM','PPNDH','PHETA1','LEUTAi','TYRTA'};
ser={'PGCD','PSP_L','PSERT'};
thr={'HSK','THRS'};
% thr={'ASPK', 'HSDy','HSK','THRS'};
trp={'TRPS1','TRPS2','TRPS3','PRAIi','IGPS','ANS','ANPRT'};
AAs=[met,arg,leu,ilv,glt,trp,his,ser,lys,thr,phetyr,aro,cys,other];

% Get reaction indices for each AA group
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

%% define parameters
fprintf('Defining model parameters...\n');
dt = 0.005;% Time step in hours
t = -2:dt:6;% Total simulation time
gamma = 8.35;% per h;
lambda_C = 1.17;%strain-specific constant of growth laws;
phi_Rb0 = 0.0445;% Minimum ribosome fraction
phi_Q = 0.45;% Housekeeping protein fraction from You et al., 2013
phi_C1 = 0.03;% Carbon catabolism protein fraction
phi_C2 = 0.03;% Carbon catabolism protein fraction
phi_AAs_max = 0.181;%from Wu et al., 2023
alpha_A = 0.11;%from Wu et al., 2023 
k=3.6;% Adaptation rate constant for sigma

% Growth rates and protein fractions from experimental data
lambda_i = 0.68;%Initial growth rate (per h)
lambda_f = 1.40;% Final growth rate (per h)
phi_Rb_i=phi_Rb0+lambda_i./gamma;% Initial ribosome fracti
sigma_i = lambda_i./phi_Rb_i; % Initial sigma
phi_Rb_f=phi_Rb0+lambda_f./gamma;% Final ribosome fraction
sigma_f = lambda_f./phi_Rb_f; % Final sigma

% Experimental AA biosynthesis protein fractions (phi_AA)
phi_met_f = 0.6*10^-3;
phi_met_i = 32.2*10^-3;
phi_arg_f = 0.14*10^-3;
phi_arg_i = 4.71*10^-3;
phi_leu_f = 0.25*10^-3;
phi_leu_i = 5.85*10^-3;
phi_trp_f = 0.12*10^-3;
phi_trp_i = 2.31*10^-3;
phi_glt_f = 0.72*10^-3;
phi_glt_i = 7.23*10^-3;
phi_ilv_f = 1.6*10^-3;
phi_ilv_i = 13.1*10^-3;
phi_his_f = 0.74*10^-3;
phi_his_i = 4.15*10^-3;
phi_ser_f = 1.0*10^-3;
phi_ser_i = 3.79*10^-3;
phi_lys_f = 2.63*10^-3;
phi_lys_i = 6.32*10^-3;
phi_thr_f = 2.4*10^-3;
phi_thr_i = 4.13*10^-3;
phi_aro_f = 1.97*10^-3;
phi_aro_i = 2.97*10^-3;
phi_phetyr_f = 0.65*10^-3;
phi_phetyr_i = 1.09*10^-3;
phi_cys_f = 8.66*10^-3;
phi_cys_i = 9.88*10^-3;
phi_other_f = 11.1*10^-3;
phi_other_i = 11.5*10^-3;
phi_AAs_total_i = phi_met_i+phi_arg_i+phi_leu_i+phi_ilv_i+phi_glt_i+phi_trp_i+phi_his_i+phi_ser_i+phi_lys_i+phi_thr_i+phi_phetyr_i+phi_aro_i+phi_cys_i+phi_other_i;
phi_AAs_total_f = phi_met_f+phi_arg_f+phi_leu_f+phi_ilv_f+phi_glt_f+phi_trp_f+phi_his_f+phi_ser_f+phi_lys_f+phi_thr_f+phi_phetyr_f+phi_aro_f+phi_cys_f+phi_other_f;

% % chi_AAs_i = phi_AAs_total_i;
phi_AAs_i = [phi_met_i,phi_arg_i,phi_leu_i,phi_ilv_i,phi_glt_i,phi_trp_i,phi_his_i,phi_ser_i,phi_lys_i,phi_thr_i,phi_phetyr_i,phi_aro_i,phi_cys_i,phi_other_i];
phi_AAs_f = [phi_met_f,phi_arg_f,phi_leu_f,phi_ilv_f,phi_glt_f,phi_trp_f,phi_his_f,phi_ser_f,phi_lys_f,phi_thr_f,phi_phetyr_f,phi_aro_f,phi_cys_f,phi_other_f];
%% Initialize arrays to store simulation results
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
lambda_FBA=nan(size(t));
v_R=nan(size(t));%protein synthesis flux
v_C=nan(size(t));
v_E=nan(size(t));%C uptake flux
v=nan(length(model.rxns),length(t));% all rxn fluxes
chi_AA_total=nan(length(AAs_group),length(t));
phi_AA_total=nan(length(AAs_group),length(t));
v_AAGroup_total=nan(length(AAs_group),length(t));
ratios_total=nan(length(AAs_group),length(t));
v_AAs_total=nan(length(AAs_r),length(t));%flux of all AA reactions

v_met_total=nan(length(met_r),length(t));
v_arg_total=nan(length(arg_r),length(t));
v_leu_total=nan(length(leu_r),length(t));
v_ilv_total=nan(length(ilv_r),length(t));
v_glt_total=nan(length(glt_r),length(t));
v_trp_total=nan(length(trp_r),length(t));
v_his_total=nan(length(his_r),length(t));
v_ser_total=nan(length(ser_r),length(t));
v_lys_total=nan(length(lys_r),length(t));
v_thr_total=nan(length(thr_r),length(t));
v_phetyr_total=nan(length(phetyr_r),length(t));
v_aro_total=nan(length(aro_r),length(t));
v_cys_total=nan(length(cys_r),length(t));
v_other_total=nan(length(other_r),length(t));

%% initialization
sigma(1) = sigma_i;
lambda(1) = lambda_i;
phi_AAs(1) = phi_AAs_total_i;
phi_Rb(1) = phi_Rb_i;
phi_E(1) = 1-phi_C1-phi_Rb(1)-phi_Q;
chi_Rb(1) = phi_Rb_i;
chi_AAs(1) = phi_AAs_total_i;
chi_AA_total(:,1)=[phi_met_i;phi_arg_i;phi_leu_i;phi_ilv_i;phi_glt_i;phi_trp_i;phi_his_i;phi_ser_i;phi_lys_i;phi_thr_i;phi_phetyr_i;phi_aro_i;phi_cys_i;phi_other_i];
phi_AA_total(:,1)=[phi_met_i;phi_arg_i;phi_leu_i;phi_ilv_i;phi_glt_i;phi_trp_i;phi_his_i;phi_ser_i;phi_lys_i;phi_thr_i;phi_phetyr_i;phi_aro_i;phi_cys_i;phi_other_i];
% Dummy initial fluxes for first dCAFBA call
v_R(1) = 4.6;
v_C(1) = 10;
v_met(1)=1;
v_arg(1)=1;
v_leu(1)=1;
v_trp(1)=1;
v_glt(1)=1;
v_ilv(1)=1;
v_his(1)=1;
v_ser(1)=1;
v_lys(1)=1;
v_thr(1)=1;
v_aro(1)=1;
v_phetyr(1)=1;
v_cys(1)=1;
v_other(1)=1;
v_AAGroup_total(:,1)=[v_met(1);v_arg(1);v_leu(1);v_ilv(1);v_glt(1);v_trp(1);v_his(1);v_ser(1);v_lys(1);v_thr(1);v_phetyr(1);v_aro(1);v_cys(1);v_other(1)];
ratios_total(:,1)=[phi_met_i/phi_met_f;phi_arg_i/phi_arg_f;phi_leu_i/phi_leu_f;phi_ilv_i/phi_ilv_f;phi_glt_i/phi_glt_f;phi_trp_i/phi_trp_f;phi_his_i/phi_his_f;phi_ser_i/phi_ser_f;phi_lys_i/phi_lys_f;phi_thr_i/phi_thr_f;phi_phetyr_i/phi_phetyr_f;phi_aro_i/phi_aro_f;phi_cys_i/phi_cys_f;phi_other_i/phi_other_f];
v_met_total(:,1)=[0.1;0.1;0.1;0.1;0.1;0.1];
v_met_total(:,1) = 0.1 * ones(1, length(met_r));
v_arg_total(:,1)= 0.1 * ones(1, length(arg_r));
v_leu_total(:,1)= 0.1 * ones(1, length(leu_r));
v_ilv_total(:,1)= 0.1 * ones(1, length(ilv_r));
v_glt_total(:,1)= 0.1 * ones(1, length(glt_r));
v_trp_total(:,1)= 0.1 * ones(1, length(trp_r));
v_his_total(:,1)= 0.1 * ones(1, length(his_r));
v_ser_total(:,1)= 0.1 * ones(1, length(ser_r));
v_lys_total(:,1)= 0.1 * ones(1, length(lys_r));
v_thr_total(:,1)= 0.1 * ones(1, length(thr_r));
v_phetyr_total(:,1)= 0.1 * ones(1, length(phetyr_r));
v_aro_total(:,1)= 0.1 * ones(1, length(aro_r));
v_cys_total(:,1)= 0.1 * ones(1, length(cys_r));
v_other_total(:,1)=0.1 * ones(1, length(other_r));
%% main simulation loop
tic;  
t0=find(t==0);
eta_total = zeros( size(phi_AA_total) );
onflag = logical( eta_total );
flagact = zeros(1, numel(t));

for i=2:length(t)
    % Pre-shift steady state (t < 0)
    if t(i)<0
        % System remains in initial steady stat
        phi_Rb(i)=phi_Rb_i;
        phi_AAs(i)=phi_AAs_total_i;
        sigma(i)=sigma_i;
        chi_Rb(i)=phi_Rb_i;
        chi_AAs(i)=phi_AAs_total_i;
        chi_AA_total(:,i)=chi_AA_total(:,1);
        phi_AA_total(:,i)=phi_AA_total(:,1); 
        ratios_total(:,i) = phi_AA_total(:,i) ./ phi_AAs_f(:);
        lambda(i)=lambda_i;
        % Run CAFBA model for pre-shift condition
       [v_R(i),v_C(i),lambda_FBA(i),v(:,i),v_E(i),v_AAGroup_total(1,i),v_AAGroup_total(2,i),v_AAGroup_total(3,i),v_AAGroup_total(4,i),v_AAGroup_total(5,i),v_AAGroup_total(6,i),v_AAGroup_total(7,i),v_AAGroup_total(8,i),v_AAGroup_total(9,i),v_AAGroup_total(10,i),v_AAGroup_total(11,i),v_AAGroup_total(12,i),v_AAGroup_total(13,i),v_AAGroup_total(14,i),phi_E(i)] = uCAFBA_AAs_core_pre(model1,Dir,phi_AAs(i-1),phi_Rb(i-1),phi_AA_total(1,i-1),phi_AA_total(2,i-1),phi_AA_total(3,i-1),phi_AA_total(4,i-1),phi_AA_total(5,i-1),phi_AA_total(6,i-1),phi_AA_total(7,i-1),phi_AA_total(8,i-1),phi_AA_total(9,i-1),phi_AA_total(10,i-1),phi_AA_total(11,i-1),phi_AA_total(12,i-1),phi_AA_total(13,i-1),phi_AA_total(14,i-1),phi_C1,w_C,w_E1,glyc_r,phi_Q,phiE_r,arg_r,aro_r,cys_r,glt_r,his_r,ilv_r,leu_r,lys_r,met_r,other_r,phetyr_r,ser_r,thr_r,trp_r,AAs_r);
               
     % Post-shift dynamic adaptation (t >= 0)   
    else 
        % Update sigma based on a step-like function
        if t(i)>= 0 && t(i)< 1/3
            sigma(i)=sigma_i+k*t(i);
        else
            sigma(i)=sigma_f;
        end  
        
        % Update protein fractions using a dynamic model
       chi_Rb(i)=phi_Rb0/(1-sigma(i)/gamma);
       chi_AAs(i)=phi_AAs_max-alpha_A*sigma(i)*chi_Rb(i);
       
       % Euler integration for protein fractions
       phi_Rb(i)=sigma(i-1)*phi_Rb(i-1)*(chi_Rb(i-1)-phi_Rb(i-1))*dt+phi_Rb(i-1);
       phi_AAs(i)=sigma(i-1)*phi_Rb(i-1)*(chi_AAs(i-1)-phi_AAs(i-1))*dt+phi_AAs(i-1);
       lambda(i)=sigma(i)*phi_Rb(i);  
       
       % Run CAFBA model for post-shift condition
       [v_R(i),v_C(i),lambda_FBA(i),v(:,i),v_E(i),v_AAGroup_total(1,i),v_AAGroup_total(2,i),v_AAGroup_total(3,i),v_AAGroup_total(4,i),v_AAGroup_total(5,i),v_AAGroup_total(6,i),v_AAGroup_total(7,i),v_AAGroup_total(8,i),v_AAGroup_total(9,i),v_AAGroup_total(10,i),v_AAGroup_total(11,i),v_AAGroup_total(12,i),v_AAGroup_total(13,i),v_AAGroup_total(14,i),phi_E(i)] = uCAFBA_AAs_core_post(model2,Dir,phi_AAs(i-1),phi_Rb(i-1),phi_AA_total(1,i-1),phi_AA_total(2,i-1),phi_AA_total(3,i-1),phi_AA_total(4,i-1),phi_AA_total(5,i-1),phi_AA_total(6,i-1),phi_AA_total(7,i-1),phi_AA_total(8,i-1),phi_AA_total(9,i-1),phi_AA_total(10,i-1),phi_AA_total(11,i-1),phi_AA_total(12,i-1),phi_AA_total(13,i-1),phi_AA_total(14,i-1),phi_C2,w_C,w_E2,glyc_r,phi_Q,phiE_r,arg_r,aro_r,cys_r,glt_r,his_r,ilv_r,leu_r,lys_r,met_r,other_r,phetyr_r,ser_r,thr_r,trp_r,AAs_r,v_met_i,v_arg_i,v_leu_i,v_ilv_i,v_glt_i,v_trp_i,v_his_i,v_ser_i,v_lys_i,v_thr_i,v_phetyr_i,v_aro_i,v_cys_i,v_other_i,ratios_total(1,i-1),ratios_total(2,i-1),ratios_total(3,i-1),ratios_total(4,i-1),ratios_total(5,i-1),ratios_total(6,i-1),ratios_total(7,i-1),ratios_total(8,i-1),ratios_total(9,i-1),ratios_total(10,i-1),ratios_total(11,i-1),ratios_total(12,i-1),ratios_total(13,i-1),ratios_total(14,i-1));   
        
        %determine which AA is limiting and update allocation

       ratios = phi_AA_total(:,i-1) ./ phi_AAs_f(:);

       chi_AA_total(:,i) = zeros(length(AAs_group),1);
       eta_total(:,i) = zeros(length(AAs_group),1);

    % Find the Turn On amino acid at the last moment
       onflagt = onflag(:, i-1);%onflag at the last moment     

    % Find next activate candidate AA
       tol = 0.005;
       
       candidateon = ratios < 1 + tol ;    
       candidateon = candidateon & ~onflagt;
       if any(candidateon)
           flagact(i) = sum(candidateon); 
           onflagt(candidateon) = true;  
       end

    % Update eata_totall
       eta_total( onflagt, i ) = phi_AAs_f(onflagt) ./ sum( phi_AAs_f( onflagt ) );
       chi_AA_total( onflagt, i )=chi_AAs(i)*eta_total( onflagt, i );  

       onflag(:, i) = onflagt;
     % Euler integration for individual AA protein fractions
       phi_AA_total(:,i)=sigma(i-1)*phi_Rb(i-1)*(chi_AA_total(:,i-1)-phi_AA_total(:,i-1))*dt+phi_AA_total(:,i-1); 
       ratios_total(:,i) = ratios;
   end
    
end 
T=toc;
%% Visualization (example plots)
% This section can be added to visualize the results, e.g., plotting
% growth rate over time.
figure;
plot(t, lambda, 'b-', 'LineWidth', 2);
xlabel('Time (h)');
ylabel('Growth Rate (h^{-1})');
title('Growth Rate during AA Upshift');