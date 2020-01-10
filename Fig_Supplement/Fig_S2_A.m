%%%%%%%%%%%%%% comparison of expected distribution and control distribution
%%%%%%%%%%%%%% of length of G1 phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

Proportion_G1=0.45;

%%%%%%%%%%%%%%%%%%% data  unstimulated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
run Import_Cycle_Data.m;       % cell cycle data of untreated cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% data stimulated G1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'Dead_Survivor_white.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

generation = dataArray{:, 1};
cell = dataArray{:, 2};
birthwhitestart= dataArray{:, 3};
whiteendgreenstart = dataArray{:, 4};
greenend = dataArray{:, 5};
roundstartMstart = dataArray{:, 6};
celldivision = dataArray{:, 7};
celldeath = dataArray{:, 8};
celllost = dataArray{:, 9};
survived = dataArray{:, 10};
time = dataArray{:, 11};
videoend = dataArray{:, 12};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

survived((celldeath-time)*0.25>19.25)=1;
celldeath((celldeath-time)*0.25>19.25)=nan;

cd Fig_Supplement
%%%%%%%%%%%%%%%%%%% data  end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% time of cell death (only cells with known history) %%%%%%%%
cd ..
run Mean_Cycle_Start.m
cd Fig_Supplement
t_TRAIL_Death_white_gen0=t_TRAIL_Death(~isnan(t_Birth_TRAIL_white) & isnan(t_TRAIL_Div_F0));
t_TRAIL_Death_white_gen1=t_TRAIL_Death(~isnan(t_Birth_TRAIL_white) & t_TRAIL_Div_F0>=0);
t_TRAIL_Death_white_weighted=[t_TRAIL_Death_white_gen0,t_TRAIL_Death_white_gen0,t_TRAIL_Death_white_gen1];
paracycle1=lognfit(t_TRAIL_Death_white_weighted);
mean_t_TRAIL_Death_white=mean(t_TRAIL_Death_white_weighted);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paracycle1(1)=1.1062;   %parameters lognormal distribution t_death G1
paracycle1(2)=0.5649;




%%%%%%%%%%%%%%%%%% Ensemble model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=log(2);
varax=[0:0.01:Proportion_G1];  
my_pdf=2*a*exp(-varax*a);    
y = rand(1000000,1);   
F_inv_cdf_eval=-log((2-y)./2)/a;
C_0=F_inv_cdf_eval;
C_0_G1=C_0(C_0<Proportion_G1);   %only stimulation in G1
g_zw=lognrnd(log_G1(1),log_G1(2),length(C_0_G1),1); %phase length untreated
g=Proportion_G1./g_zw;                              %growth rate untreated
T_death = lognrnd(paracycle1(1),paracycle1(2), length(C_0_G1),1);
C_end= C_0_G1 + g.*T_death;
Length_Phase_Uncensored=g_zw(C_end>Proportion_G1);
New_Para_white_MODEL=lognfit([Length_Phase_Uncensored]');

%%%%%%%%%%%%%%% expected cell cycle length in our scenario (with survivors and apoptotic), without censored cells
Faktor_Survivors=length(cell(~isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart) & survived==1))/length(cell(~isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart) & survived~=1));
Length_Phase_Survivors=lognrnd(log_G1(1),log_G1(2),(round(Faktor_Survivors*length(Length_Phase_Uncensored))) ,1);
Para_white_MODEL_Surv_AND_Apop=lognfit([Length_Phase_Uncensored;Length_Phase_Survivors]);

Sample_old=lognrnd(log_G1(1),log_G1(2),length(Length_Phase_Uncensored),1); %without TRAIL
mean_new=exp(New_Para_white_MODEL(1) + (New_Para_white_MODEL(2)^2)/2)
mean(Length_Phase_Uncensored)
mean_old=exp(log_G1(1) + (log_G1(2)^2)/2)

figure()
histogram(Sample_old,'FaceColor',[0,0,0]);
hold on;
histogram(Length_Phase_Uncensored,'FaceColor',[1,1,1]);
set(gca,'FontSize',20);
xlim([2,15]);
ylim([0,8000]);
leg{1}=['control (mean=7.4h)'];  
leg{2}=['estimated (mean=6.9h)'];
h=legend(leg, 'Location','northeast');
set(h,'FontSize',16);
xlabel('G_1 [h]');
ylabel('# cells');
grid on;













