%%%%%%%%%%%%% influence of different t_death distributions on expected
%%%%%%%%%%%%% distributions of phase lengths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

Proportion_G1=0.45;

%%%%%%%%%%%%%%%%%%% data  unstimulated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
run Import_Cycle_Data.m;       %% cell cycle data of untreated cells
cd Fig_Supplement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=log(2);
y = rand(1000000,1);  
F_inv_cdf_eval=-log((2-y)./2)/a;
C_0=F_inv_cdf_eval;
C_0_G1=C_0(C_0<=Proportion_G1);
g_zw=lognrnd(log_G1(1),log_G1(2),length(C_0_G1),1);
g=Proportion_G1./g_zw;   


paracycle_mu=[0.3, 0.8, 1.3, 1.8,  2.1];   %parameters for different t_death distributions
paracycle_sigma=[0.4, 0.5, 0.6, 0.7, 0.8];

Mean_old=exp(log_G1(1) + ((log_G1(2)^2)/2));

Delta_Mean=zeros(5);

for i=1:5
    for j=1:5
        paracycle1(1)=paracycle_mu(i)   ;   %t_death only G1
        paracycle1(2)=paracycle_sigma(j);
%%%%%%%%%%%%%%% ensemble model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T_death = lognrnd(paracycle1(1),paracycle1(2), length(C_0_G1),1);
        C_end= C_0_G1 + g.*T_death;
        Length_Phase_Uncensored=g_zw(C_end>Proportion_G1);
        New_Para_white_MODEL=lognfit([Length_Phase_Uncensored]');
        Mean_new=exp(New_Para_white_MODEL(1) + ((New_Para_white_MODEL(2)^2)/2));
        Delta_Mean(i,j)=Mean_old-Mean_new;
    end
end

Delta_Mean=round(100*Delta_Mean)/100;

map = [0.8,1,1
    0.7,0.9,1
    0.6 0.8 1
    0.5 0.7 1
    0.4,0.6,1
    0.3 0.5 1
    0.2,0.4,1
    0.1 0.3 1
    ];

figure()
xvalues = {'0.3', '0.8', '1.3', '1.8','2.1'};
yvalues = {'0.4', '0.5', '0.6', '0.7', '0.8'};
h=heatmap(xvalues,yvalues,Delta_Mean', 'FontSize',20,'colormap', map);
h.XLabel = 'mu';
h.YLabel = 'sigma'; 






