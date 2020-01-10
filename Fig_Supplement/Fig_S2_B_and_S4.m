%%%%%%%%%%%%%% comparison of expected distribution and control distribution
%%%%%%%%%%%%%% of length of G1 phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% and validation with 'untreated' cells %%%%%%%%%%%%%%%%%%%%%%

clear all

Proportion_G1=0.45;

%%%%%%%%%%%%%%%%%%% data unstimulated %%%%%%%%%%%%
cd ..
run Import_Cycle_Data.m;       %% cell cycle data of untreated cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% data  stimulated S/G2/M  %%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'Dead_Survivor_green.csv';
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

survived((celldeath-time)*0.25>19.25)=1;
celldeath((celldeath-time)*0.25>19.25)=nan;

clearvars filename delimiter startRow formatSpec fileID dataArray ans;

cd Fig_Supplement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd ..
run Mean_Cycle_Start.m
cd Fig_Supplement


paracycle1(1)=1.2506;  %parameters lognormal distribution for S/G2/M
paracycle1(2)=0.5730;



%%%%%%%%%%% TEST if G1 'untreated' cells are like the other untreated cells %%%%%%%%%%%%
TEST_Uncens_Length_G1_F0 = (whiteendgreenstart(~isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart)) - birthwhitestart(~isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart)) )*0.25;
mean(G1(cens_G1==0));
mean(TEST_Uncens_Length_G1_F0);
[probTest,hypTest]=ranksum(TEST_Uncens_Length_G1_F0',G1(cens_G1==0)')

figure()
h1=histogram(TEST_Uncens_Length_G1_F0,50)
set(h1,'FaceColor',[0.6,0.6,0.6],'FaceAlpha',1) 
hold on;
h2=histogram(G1(cens_G1==0),50);
set(h2,'FaceColor',[0.4,0.6,1],'FaceAlpha',0.5) 
xlabel('G1 [h]');
ylabel('# cells');
set(gca,'FontSize',20);
leg{1}=['control'];  
leg{2}=['treated in S/G_2/M'];
h=legend(leg, 'Location','northeast');
set(h,'FontSize',18);
grid on;
ax = gca;
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [5 4]);  
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 5 4]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%% Ensemble model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=log(2);
y = rand(1000000,1);  
F_inv_cdf_eval=-log((2-y)./2)/a;
C_0=F_inv_cdf_eval;
C_0_SG2M=C_0(C_0>=Proportion_G1);
g_zw=lognrnd(log_SG2M(1),log_SG2M(2),length(C_0_SG2M),1);
g=(1-Proportion_G1)./g_zw;   
T_death = lognrnd(paracycle1(1),paracycle1(2), length(C_0_SG2M),1);
C_end= C_0_SG2M + g.*T_death;
Length_Phase_Uncensored=g_zw(C_end>1);
New_Para_green_MODEL=lognfit([Length_Phase_Uncensored]);

%%%%%%%%%%%%%%% expected cell cycle length in our scenario with survivors and apoptotic (without censored cells)
Faktor_Survivors=length(cell(~isnan(whiteendgreenstart) & generation==0 & ~isnan(celldivision) & survived==1))/length(cell(~isnan(whiteendgreenstart) & generation==0 & ~isnan(celldivision) & survived~=1));
Length_Phase_Survivors=lognrnd(log_SG2M(1),log_SG2M(2),(round(Faktor_Survivors*length(Length_Phase_Uncensored))) ,1);
Para_green_MODEL_Surv_AND_Apop=lognfit([Length_Phase_Uncensored;Length_Phase_Survivors]);


Sample_old=lognrnd(log_SG2M(1),log_SG2M(2),length(Length_Phase_Uncensored),1); %without TRAIL
mean_new=exp(New_Para_green_MODEL(1) + (New_Para_green_MODEL(2)^2)/2);
mean_old=exp(log_SG2M(1) + (log_SG2M(2)^2)/2);


figure()
histogram(Sample_old,'FaceColor',[0,0.4,0]);
hold on;
histogram(Length_Phase_Uncensored,'FaceColor',[1,1,1]);
set(gca,'FontSize',20);
xlim([3,16]);
ylim([0,5500]);
leg{1}=['control (mean=9.0h)'];  
leg{2}=['estimated (mean=8.7h)'];
h=legend(leg, 'Location','northeast');
set(h,'FontSize',16);
xlabel('S/G_2/M [h]');
ylabel('# cells');
grid on;


