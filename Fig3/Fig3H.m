%%%%%%%%%%%% modeled phase lengths over time of TRAIL addition in S/G2/M %%%
%%%%%%%%%%%% Attention: example file is loaded %%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

%%%%%%%%%%%%%%%%%%% data  unstimulated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
run Import_Cycle_Data.m;               %cell cycle data of untreated cells
cd Fig3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Faktor_Survivors=0.3831;  %here: survivors/All_Apoptotic  0.1832 for HCT
Proportion_G1=0.45;       %0.39 for HCT
paracycle1(1)=1.2506;     %tdeath SG2M, data   2.3214 for HCT
paracycle1(2)=0.5730;     %0.5956 for HCT 
z=3.4;                    %1.95 for HCT
n=1200;

%%%%%%%%% Ensemble model %%%%%%%%%%%%%%%

a=log(2);
y = rand(n,1); 
F_inv_cdf_eval=-log((2-y)./2)/a;
C_0=F_inv_cdf_eval;
C_0_SG2M=C_0(C_0>=Proportion_G1);
exp(log_SG2M(1) + (log_SG2M(2)^2)/2)
g_zw=lognrnd(log_SG2M(1),log_SG2M(2),length(C_0_SG2M),1);
g=(1-Proportion_G1)./g_zw;   
g_zw_NEW=g_zw+z;
g_NEW=(1-Proportion_G1)./g_zw_NEW;   
T_death = lognrnd(paracycle1(1),paracycle1(2), length(C_0_SG2M),1);
C_end= C_0_SG2M + g_NEW.*T_death;
t_0= (C_0_SG2M - Proportion_G1)./g;         %time(0): start of S/G2/M
t_SG2M_new = t_0 + (1-C_0_SG2M)./(g_NEW);
t_0_without_censored=t_0(C_end>1);
t_SG2M_new_without_censored=t_SG2M_new(C_end>1);


%%%%%%%%%%%%%%%Including Survivors%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_Surv=round(Faktor_Survivors*n)   
a=log(2);
y = rand(n_Surv,1);  
F_inv_cdf_eval=-log((2-y)./2)/a;
C_0=F_inv_cdf_eval;
C_0_SG2M_Surv=C_0(C_0>=Proportion_G1);
g_zw=lognrnd(log_SG2M(1),log_SG2M(2),length(C_0_SG2M_Surv),1);
g=(1-Proportion_G1)./g_zw;   %
g_zw_NEW=g_zw+z;
g_NEW=(1-Proportion_G1)./g_zw_NEW;
t_0_Surv= (C_0_SG2M_Surv - Proportion_G1)./g;
t_SG2M_new_Surv = t_0_Surv + (1-C_0_SG2M_Surv)./(g_NEW);
t_0_without_censored=[t_0_without_censored;t_0_Surv];
t_SG2M_new_without_censored=[t_SG2M_new_without_censored;t_SG2M_new_Surv];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%save('model_data_example_SG2M.mat','t_0_without_censored','t_SG2M_new_without_censored')
load('model_data_example_SG2M.mat')

uni=linspace(0,16,(length(t_0_without_censored)))';
x=t_0_without_censored;
y=t_SG2M_new_without_censored;
gprMdl=fitrgp(x,y,'Basis','linear',...
      'FitMethod','exact','PredictMethod','exact');
[ypred_Data,ysd_Data,conf] = predict(gprMdl,uni);

c1=[0:18];
figure()
h = area(uni,conf,'LineStyle',':');
h(1).FaceColor = [1 1 1];
h(1).FaceAlpha =0;
h(2).FaceColor = [0.8 1 0.8];
h(2).FaceAlpha =0.5;
hold on;
plot(uni,ypred_Data,'k','LineWidth',1.5);
hold on;
scatter(t_0_without_censored,t_SG2M_new_without_censored,80,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0.4,0],'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
grid on;
box on;
set(gca,'FontSize',28);
xlabel('TRAIL in S/G_2/M [h]');
ylabel('S/G_2/M [h]');
xlim([0,12.25]);
ylim([0,25]);   
