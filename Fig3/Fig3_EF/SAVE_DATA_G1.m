clear all

%%%%%%%%%%%%%%%%%%% data  unstimulated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
cd ..
run Import_Cycle_Data.m;       %% cell cycle data of untreated cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% data  stimulated in G1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

cd Fig3/Fig3_EF

%%%%%% 1.: end of G1 phase is known
Uncens_Length_G1_F0 = (whiteendgreenstart(~isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart)) - birthwhitestart(~isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart)) )*0.25;
Time_TRAIL_G1 =(time(~isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart)) - birthwhitestart(~isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart)) )*0.25;
save('Data_real_G1.mat','Time_TRAIL_G1','Uncens_Length_G1_F0');






