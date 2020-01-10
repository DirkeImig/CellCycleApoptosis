%%%%%%%%%%%% experimental data of phase lengths over time for TRAIL addition in S/G2/M %%%%%%

clear all

%%%%%%%%%%%%%%%%%%% data  unstimulated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
run Import_Cycle_Data.m;       %% cell cycle data of untreated cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% data  stimulated in S/G2/M  %%%%%%%%%%%%%%%%%%%%%%%%
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
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
cd Fig3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% 1: cell division: end of phase is known
Uncens_Length_SG2M_F0 = (celldivision(~isnan(celldivision) & generation==0 & ~isnan(whiteendgreenstart)) - whiteendgreenstart(~isnan(celldivision) & generation==0 & ~isnan(whiteendgreenstart)) )*0.25;
Time_TRAIL_SG2M = (time(~isnan(celldivision) & generation==0 & ~isnan(whiteendgreenstart)) - whiteendgreenstart(~isnan(celldivision) & generation==0 & ~isnan(whiteendgreenstart)) )*0.25;

%%%%%% gaussian process for confidence interval
uni=unique(Time_TRAIL_SG2M);
x=Time_TRAIL_SG2M;
y=Uncens_Length_SG2M_F0;
gprMdl=fitrgp(x,y,'Basis','linear',...
      'FitMethod','exact','PredictMethod','exact');
[ypred_Data,ysd_Data,conf] = predict(gprMdl,uni);

c1=[0:17];
figure()
h = area(uni,conf,'LineStyle',':');
h(1).FaceColor = [1 1 1];
h(1).FaceAlpha =0;
h(2).FaceColor = [0.8 1 0.8];
h(2).FaceAlpha =0.5;
hold on;
plot(uni,ypred_Data,'k','LineWidth',1.5);
hold on;
scatter(Time_TRAIL_SG2M,Uncens_Length_SG2M_F0,80,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0.4,0],'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
grid on;
box on;
set(gca,'FontSize',28);
xlabel('TRAIL in S/G_2/M [h]');
ylabel('S/G_2/M [h]');
xlim([0,max(Time_TRAIL_SG2M)]);
ylim([0,25]);




