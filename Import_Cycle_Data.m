filename = 'Untreated_CellCycle.csv';
%filename = 'Untreated_CellCycle_HCT.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
CC =   dataArray{:, 1};
G1 =   dataArray{:, 2};
SG2 =  dataArray{:, 3};
M =    dataArray{:, 4};
cens = dataArray{:, 5};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%%%%%% cens: information about censoring 
%%%%%% cens: not changed
cens_Cycle= cens;                   
cens_Cycle(isnan(cens_Cycle))=0;   %if nan: no censoring
cens_Cycle(cens_Cycle~=0) = 1;     %if number, censoring of cycle
CC(isnan(CC))=cens(isnan(CC));     %how long was cell cycle at least (censoring)
%%%%% CC contains both, censored and non-censored data %%%%%%%%%%%%%%%%%%%%%

cens_G1= cens;
cens_G1(G1<0) = 1;
cens_G1(G1>0) = 0;
cens_G1(G1==0) = nan;   %no information about G1
G1(G1<0)=cens(G1<0); 
G1(G1==0)= nan;  %not defined, beginning not known.              

SG2M=SG2+M;                           
cens_SG2M= cens;
cens_SG2M(SG2M==0) = nan;             %no information about SG2M (G1 was already censored)
cens_SG2M(SG2M<0) = 1;                %censoring in SG2M
cens_SG2M(SG2M>0) = 0;                %no censoring
SG2M(SG2M<0)=cens(SG2M<0)-G1(SG2M<0); %how long was S/G2/M at least? Censoring minus G1. 
SG2M(SG2M==0)= nan;                   %not defined, beginning not known.                   

%%%%% exclude censored G1 data if G1>35 (non-cycling) --> for cycle and G1
%%%%% exclude censored SG2M data if SG2M>30 --> for SG2M
%%%%% exclude G1 cells where only S/G2/M was measured
log_cycle=lognfit(CC(~(G1>35 & cens_G1==1) ),0.05,cens_Cycle(~(G1>35 & cens_G1==1)));
log_G1=lognfit(G1(~isnan(G1) & ~(G1>35 & cens_G1==1)),0.05,cens_G1(~isnan(G1) & ~(G1>35 & cens_G1==1)));  
log_SG2M=lognfit(SG2M(~isnan(SG2M) & ~(SG2M>30 & cens_SG2M==1)),0.05,cens_SG2M(~isnan(SG2M) & ~(SG2M>30 & cens_SG2M==1)));      %in HCT: outliers if >30

