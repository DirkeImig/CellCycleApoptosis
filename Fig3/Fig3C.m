%%%%%%%%%%%% experimental data of phase lengths over time for TRAIL addition in G1 %%%%%%

clear all

max_Exp_time=19.25; 


%%%%%%%%%%%%%%%%%%% data unstimulated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
run Import_Cycle_Data.m;       %% cell cycle data of untreated cells
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

survived((celldeath-time)*0.25>max_Exp_time)=1;
celldeath((celldeath-time)*0.25>max_Exp_time)=nan;

clearvars filename delimiter startRow formatSpec fileID dataArray ans;
cd Fig3


%%%%%% end of G1 phase is known %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Uncens_Length_G1_F0 = (whiteendgreenstart(~isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart)) - birthwhitestart(~isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart)) )*0.25;
Time_TRAIL_G1 =(time(~isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart)) - birthwhitestart(~isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart)) )*0.25;

%%%%%%%% estimate confidence inteval with linear gaussian process %%%%%%%%%
uni=unique(Time_TRAIL_G1);
x=Time_TRAIL_G1;
y=Uncens_Length_G1_F0;
gprMdl=fitrgp(x,y,'Basis','linear',...
      'FitMethod','exact','PredictMethod','exact');
[ypred_Data,ysd_Data,conf] = predict(gprMdl,uni);


c1=[0:20];
figure()
h = area(uni,conf,'LineStyle',':');
h(1).FaceColor = [1 1 1];
h(1).FaceAlpha = 0;
h(2).FaceColor = [0.8 0.8 0.8];
h(2).FaceAlpha = 0.5;
hold on;
plot(uni,ypred_Data,'k','LineWidth',1.5);
hold on;
scatter(Time_TRAIL_G1,Uncens_Length_G1_F0,80,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0.4,0.4,0.4],'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
set(gca,'FontSize',28);
xlim([0,max(Time_TRAIL_G1)]);
ylim([0,25]);
hold on;
xlabel('TRAIL in G_1 [h]');
ylabel('G_1 [h]');
grid on;




