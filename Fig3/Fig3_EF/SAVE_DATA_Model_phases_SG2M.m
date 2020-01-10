clear all

%%%%%%%%%%%%%%%%%%% data  unstimulated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
cd ..
run Import_Cycle_Data.m;       %% cell cycle data of untreated cells
cd Fig3/Fig3_EF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Faktor_Survivors=0.3831;  %here survivors/ALL_apoptotic
Proportion_G1=0.45;
paracycle1(1)=1.2506;      %tdeath SG2M
paracycle1(2)=0.5730;

%%%%%%%%% Ensemble model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numb=20;
z=linspace(0,10,numb);
n=30000000;

for i=1:length(z)
    a=log(2);
    y = rand(n,1); 
    F_inv_cdf_eval=-log((2-y)./2)/a;
    C_0=F_inv_cdf_eval;
    C_0_SG2M=C_0(C_0>=Proportion_G1);
    g_zw=lognrnd(log_SG2M(1),log_SG2M(2),length(C_0_SG2M),1);
    g=(1-Proportion_G1)./g_zw;   
    g_zw_NEW=g_zw+z(i);
    g_NEW=(1-Proportion_G1)./g_zw_NEW;   
    T_death = lognrnd(paracycle1(1),paracycle1(2), length(C_0_SG2M),1);
    C_end= C_0_SG2M + g_NEW.*T_death;
    t_0= (C_0_SG2M - Proportion_G1)./g;         %time(0): start of S/G2/M
    t_SG2M_new = t_0 + (1-C_0_SG2M)./(g_NEW);
    t_0_without_censored=t_0(C_end>1);
    t_SG2M_new_without_censored=t_SG2M_new(C_end>1);

    %%%%%%%%%%%%%%%Including Survivors%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_Surv=round(Faktor_Survivors*n);   %fraction of survivors
    a=log(2);
    y = rand(n_Surv,1);  
    F_inv_cdf_eval=-log((2-y)./2)/a;
    C_0=F_inv_cdf_eval;
    C_0_SG2M_Surv=C_0(C_0>=Proportion_G1);
    g_zw=lognrnd(log_SG2M(1),log_SG2M(2),length(C_0_SG2M_Surv),1);
    g=(1-Proportion_G1)./g_zw;   %
    g_zw_NEW=g_zw+z(i);
    g_NEW=(1-Proportion_G1)./g_zw_NEW;
    t_0_Surv= (C_0_SG2M_Surv - Proportion_G1)./g;
    t_SG2M_new_Surv = t_0_Surv + (1-C_0_SG2M_Surv)./(g_NEW);
    t_0_without_censored =[t_0_without_censored;t_0_Surv];
    t_SG2M_new_without_censored =[t_SG2M_new_without_censored;t_SG2M_new_Surv];
    my_field = strcat('z',num2str(round(z(i)*100)));
    t_0_model.(my_field)=round(t_0_without_censored*4)/4;
    t_SG2M_model.(my_field)=t_SG2M_new_without_censored;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('Data_Model_SG2M_0.mat','t_0_model','t_SG2M_model');



