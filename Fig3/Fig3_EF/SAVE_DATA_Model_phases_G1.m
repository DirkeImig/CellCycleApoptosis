clear all

%%%%%%%%%%%%%%%%%%% data  unstimulated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
cd ..
run Import_Cycle_Data.m;       %% cell cycle data of untreated cells
cd Fig3/Fig3_EF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Faktor_Survivors=0.3201; %survivors/all apoptotic
Proportion_G1=0.45;
paracycle1(1)=1.1062;   %t_death only G1
paracycle1(2)=0.5649;


%%%%%%%%% Ensemble model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numb=20;
z=linspace(0,10,numb);
n=30000000;   %number of initial cells




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(z)
    a=log(2);
    y = rand(n,1);  
    F_inv_cdf_eval=-log((2-y)./2)/a;
    C_0=F_inv_cdf_eval;                                 %intitial position in cell cycle, all cells
    C_0_G1=C_0(C_0<Proportion_G1);                      %intitial position in cell cycle, only G1 cells
    g_zw=lognrnd(log_G1(1),log_G1(2),length(C_0_G1),1); 
    g=Proportion_G1./g_zw;                              %growth rate in G1
    g_zw_NEW=g_zw+z(i);                                 %prolongation of cell cycle length
    g_NEW=Proportion_G1./g_zw_NEW;   
    T_death = lognrnd(paracycle1(1),paracycle1(2), length(C_0_G1),1);
    C_end= C_0_G1 + g_NEW.*T_death;
    t_0= (C_0_G1)./g;                                   %time(0) means start of G1
    t_G1_new = t_0 + (Proportion_G1-C_0_G1)./(g_NEW);
    t_0_without_censored=t_0(C_end>Proportion_G1);
    t_G1_new_without_censored=t_G1_new(C_end>Proportion_G1);

%%%%%%%%%%%%%%% Including Survivors%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_Surv=round(Faktor_Survivors*n);  %fraction of survivor from data
    a=log(2);
    y = rand(n_Surv,1);  
    F_inv_cdf_eval=-log((2-y)./2)/a;
    C_0=F_inv_cdf_eval;
    C_0_G1_Surv=C_0(C_0<Proportion_G1);
    g_zw=lognrnd(log_G1(1),log_G1(2),length(C_0_G1_Surv),1);
    g=Proportion_G1./g_zw;   
    g_zw_NEW=g_zw+z(i);
    g_NEW=Proportion_G1./g_zw_NEW;
    t_0_Surv= (C_0_G1_Surv)./g; 
    t_G1_new_Surv = t_0_Surv + (Proportion_G1-C_0_G1_Surv)./(g_NEW);
    t_0_without_censored=[t_0_without_censored;t_0_Surv];
    t_G1_new_without_censored=[t_G1_new_without_censored;t_G1_new_Surv];
    my_field = strcat('z',num2str(round(z(i)*100)));
    t_0_model.(my_field)=round(t_0_without_censored*4)/4;
    t_G1_model.(my_field)=t_G1_new_without_censored;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('Data_Model_G1_0.mat','t_0_model','t_G1_model');



