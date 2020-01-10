clear all

%%%%%%%%%%%%%%%%%%% data  unstimulated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
cd ..
run Import_Cycle_Data.m;       %% cell cycle data of untreated cells

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
cd Fig3/Fig3_EF

%%%%%% 1: cell division: end of phase is known
Uncens_Length_SG2M_F0 = (celldivision(~isnan(celldivision) & generation==0 & ~isnan(whiteendgreenstart)) - whiteendgreenstart(~isnan(celldivision) & generation==0 & ~isnan(whiteendgreenstart)) )*0.25;
Time_TRAIL_SG2M = (time(~isnan(celldivision) & generation==0 & ~isnan(whiteendgreenstart)) - whiteendgreenstart(~isnan(celldivision) & generation==0 & ~isnan(whiteendgreenstart)) )*0.25;
save('Data_real_SG2M.mat','Time_TRAIL_SG2M','Uncens_Length_SG2M_F0');











