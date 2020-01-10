%%%%%%%%%%%% modeled phase lengths over time for TRAIL addition in G1 %%%%%%
%%%%%%%%%%%% Attention: example file is loaded %%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

%%%%%%%%%%%%%%%%%%% data  unstimulated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
run Import_Cycle_Data.m;       %cell cycle data of untreated cells
cd Fig3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z=3.6;                     %prolongation of cell cycle, in hours, 2.5 for HCT
n=600;                     %number of cells

Faktor_Survivors=0.3201;   %here: survivors/All_Apoptotic from data   0.2639 for HCT
Proportion_G1=0.45;        %0.39 for HCT
paracycle1(1)=1.1062;      %t_death G1    HCT: 2.2190
paracycle1(2)=0.5649;      %HCT: 0.6724
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Ensemble model %%%%%%%%%%%%%%%
a=log(2);
y = rand(n,1);  
F_inv_cdf_eval=-log((2-y)./2)/a;
C_0=F_inv_cdf_eval;               %intitial position in cell cycle, all cells
C_0_G1=C_0(C_0<=Proportion_G1);   %intitial position in cell cycle, only G1 cells
g_zw=lognrnd(log_G1(1),log_G1(2),length(C_0_G1),1); %
g=Proportion_G1./g_zw;                              %growth rate in G1
g_zw_NEW=g_zw+z;                                    %prolongation of cell cycle length
g_NEW=Proportion_G1./g_zw_NEW;   
T_death = lognrnd(paracycle1(1),paracycle1(2), length(C_0_G1),1);
C_end= C_0_G1 + g_NEW.*T_death;
t_0= (C_0_G1)./g;                                         %time(0): start of G1
t_G1_new = t_0 + (Proportion_G1-C_0_G1)./(g_NEW);
t_0_without_censored=t_0(C_end>Proportion_G1);
t_G1_new_without_censored=t_G1_new(C_end>Proportion_G1);  %length of G1 in uncensored cells


%%%%%%%%%%%%%%%Including Survivors%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=log(2);
y = rand(round(n*Faktor_Survivors),1);  
F_inv_cdf_eval=-log((2-y)./2)/a;
C_0=F_inv_cdf_eval;               %intitial position in cell cycle, all cells
C_0_G1_Surv=C_0(C_0<=Proportion_G1); 
g_zw=lognrnd(log_G1(1),log_G1(2),length(C_0_G1_Surv),1);
g=Proportion_G1./g_zw;           
g_zw_NEW=g_zw+z;
g_NEW=Proportion_G1./g_zw_NEW;
t_0_Surv= (C_0_G1_Surv)./g; 
t_G1_new_Surv = t_0_Surv + (Proportion_G1-C_0_G1_Surv)./(g_NEW);
t_0_without_censored=[t_0_without_censored;t_0_Surv];
t_G1_new_without_censored=[t_G1_new_without_censored;t_G1_new_Surv];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save('model_data_example_G1.mat','t_0_without_censored','t_G1_new_without_censored')
load('model_data_example_G1.mat')

%%%%%%%%%%%% estimate confidence bounds with linear gaussian process %%%%%%
uni=[0:0.1:18]';
x=t_0_without_censored;
y=t_G1_new_without_censored;
gprMdl=fitrgp(x,y,'Basis','linear',...
      'FitMethod','exact','PredictMethod','exact');
[ypred_Data,ysd_Data,conf] = predict(gprMdl,uni);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c1=[0:18];
figure()
h = area(uni,conf,'LineStyle',':');
h(1).FaceColor = [1 1 1];
h(1).FaceAlpha = 0;
h(2).FaceColor = [0.8 0.8 0.8];
h(2).FaceAlpha = 0.5;
hold on;
plot(uni,ypred_Data,'k','LineWidth',1.5);
hold on;
scatter(t_0_without_censored,t_G1_new_without_censored,80,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0.4,0.4,0.4],'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
grid on;
box on;
set(gca,'FontSize',28);
xlabel('TRAIL in G_1 [h]');
ylabel('G_1 [h]');
xlim([0,18]);
ylim([0,25]);



