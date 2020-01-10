%%%%%%%%%% p-values for apoptosis or survival in dependence on C_0 %%%%%%%%

clear all

%%%%%%%%%%%%%%%%%%% data  unstimulated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
run Import_Cycle_Data.m;       %cell cycle data of untreated cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% Data untreated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_t0=[];
t_CycleF0=[];
t_white_Zw = lognrnd(log_G1(1), log_G1(2), 1, 1000);     %length G1 phase of untreated cells
t_green_Zw= lognrnd(log_SG2M(1), log_SG2M(2), 1, 1000);  %length SG2M phase of untreated cells
t_cycle_Zw=lognrnd(log_cycle(1),log_cycle(2), 1, 1000);
C_WE_zw = 0.45;
x_G1=3.6;
x_SG2M=3.4;
max_Exp_time=19.25; 
pdflog_white=makedist('Lognormal','mu',log_G1(1),'sigma',log_G1(2));
pdflog_green=makedist('Lognormal','mu',log_SG2M(1),'sigma',log_SG2M(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% Szenario 1: TRAIL in G1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'Dead_Survivor_white.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
generation = dataArray{:, 1};
birthwhitestart = dataArray{:, 3};
whiteendgreenstart = dataArray{:, 4};
survived = dataArray{:, 10};
time = dataArray{:, 11};
celldeath=dataArray{:, 8};

survived((celldeath-time)*0.25>max_Exp_time)=1;
celldeath((celldeath-time)*0.25>max_Exp_time)=nan;

clearvars filename delimiter startRow formatSpec fileID dataArray ans;
cd Fig4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t_Birth_TRAIL_white=[];   %TRAIL after Birth
Birth_to_WhiteEnd =[];
survival_white=[];

numb1=0;
for i=1:length(generation)
    if generation(i) == 0 && ~isnan(birthwhitestart(i))
        t_Birth_TRAIL_white=[t_Birth_TRAIL_white, (time(i)-birthwhitestart(i))*0.25];
        survival_white=[survival_white,survived(i)];
        Birth_to_WhiteEnd =[Birth_to_WhiteEnd,(whiteendgreenstart(i)-birthwhitestart(i))*0.25];
    elseif generation(i) == 0 && isnan(birthwhitestart(i))
        numb1=numb1+1;
    end
end


Mean_CycleStart_TR_in_white=[];
Std_CycleStart_TR_in_white=[];

for i = 1:length(t_Birth_TRAIL_white)  
    if t_Birth_TRAIL_white(i)==0
        Mean_CycleStart_TR_in_white(i)=0;
        Std_CycleStart_TR_in_white(i)=0;
    elseif Birth_to_WhiteEnd(i) > 0
        p=-(C_WE_zw/t_Birth_TRAIL_white(i) - Birth_to_WhiteEnd(i)*C_WE_zw/(x_G1*t_Birth_TRAIL_white(i)));
        q=-(C_WE_zw^2)/(x_G1*t_Birth_TRAIL_white(i));
        g1=-p/2+sqrt((p/2)^2 -q);
        C_TR_zw=t_Birth_TRAIL_white(i)*g1;   
        if isnan(C_TR_zw)
            i
        end
        Mean_CycleStart_TR_in_white(i)=mean(C_TR_zw);
        Std_CycleStart_TR_in_white(i)=std(C_TR_zw);
    else
        Trunc_gen_t_white=truncate(pdflog_white,t_Birth_TRAIL_white(i),Inf);
        t_white_rand = random(Trunc_gen_t_white,1000,1);
        C_TR_zw = t_Birth_TRAIL_white(i)*mean(C_WE_zw)./t_white_rand;  
        Std_CycleStart_TR_in_white=[Std_CycleStart_TR_in_white, std(C_TR_zw)];
        Mean_CycleStart_TR_in_white=[Mean_CycleStart_TR_in_white, mean(C_TR_zw)];
    end
end





%%%%%%%%%%%% Sceanrio II: TRAIL in SG2M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
filename = 'Dead_Survivor_green.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
generation1 = dataArray{:, 1};
numb=dataArray{:, 2};
whiteendgreenstart1 = dataArray{:, 4};
celldeath=dataArray{:, 8};
lost=dataArray{:, 9};
celldivision = dataArray{:, 7};
survival = dataArray{:, 10};
time = dataArray{:, 11};

survival((celldeath-time)*0.25>max_Exp_time)=1;
celldeath((celldeath-time)*0.25>max_Exp_time)=nan;
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
cd Fig4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% delete all lost cells %%%%%%%%%%%
generation1(isnan(survival) & ~isnan(lost))=[];
numb(isnan(survival) & ~isnan(lost))=[];
whiteendgreenstart1(isnan(survival) & ~isnan(lost))=[];
celldeath(isnan(survival) & ~isnan(lost))=[];
celldivision(isnan(survival) & ~isnan(lost))=[];
time(isnan(survival) & ~isnan(lost))=[];
survival(isnan(survival) & ~isnan(lost))=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_greenStart_to_Trail=[];   %Length of SG2M at TRAIL addition
survival_green=[];
LengthGreen =[];
numb2=0;

for i=1:length(generation1)
    if generation1(i) == 0  && ~isnan(whiteendgreenstart1(i))
        t_greenStart_to_Trail=[t_greenStart_to_Trail, (time(i)-whiteendgreenstart1(i))*0.25];
        survival_green=[survival_green,survival(i)];
        LengthGreen =[LengthGreen,(celldivision(i)-whiteendgreenstart1(i))*0.25];
    elseif generation1(i) == 0  && isnan(whiteendgreenstart1(i))
        numb2=numb2+1;
    end
end

Mean_CycleStart_TR_in_green=[];
Std_CycleStart_TR_in_green=[];



%%%%%%%%%%%%% if length green is changed by TRAIL %%%%%%%%%%%%%%%%%%%%%
for i = 1:length(t_greenStart_to_Trail)
    if t_greenStart_to_Trail(i)==0
        Std_CycleStart_TR_in_green(i)=0;
        Mean_CycleStart_TR_in_green(i)=C_WE_zw;
    elseif LengthGreen(i)>0
        p=-((1-C_WE_zw)/t_greenStart_to_Trail(i) - (1-C_WE_zw)*LengthGreen(i)/(x_SG2M*t_greenStart_to_Trail(i)));
        q=-((1-C_WE_zw)^2)/(x_SG2M*t_greenStart_to_Trail(i));
        g1=-p/2+sqrt((p/2)^2 -q);
        C_TR_zw=C_WE_zw + t_greenStart_to_Trail(i)*g1;  
        Mean_CycleStart_TR_in_green(i)=mean(C_TR_zw);
        Std_CycleStart_TR_in_green(i)=std(C_TR_zw);
    else
        Trunc_gen_t_green=truncate(pdflog_green,t_greenStart_to_Trail(i),Inf);
        t_green_rand = random(Trunc_gen_t_green,1000,1);
        C_TR_zw = mean(C_WE_zw) + (t_greenStart_to_Trail(i)) * ( (1-mean(C_WE_zw) )./(t_green_rand) );  
        Std_CycleStart_TR_in_green(i)=std(C_TR_zw);
        Mean_CycleStart_TR_in_green(i)=mean(C_TR_zw);
    end
end

%%%%%%%%%%%%%%%%%% UNTIL HERE: CYCLE START%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%% cycle start in apoptotic and survivors %%%%%%%%%%%%%%%%%
Mean_CycleStart_TR_in_white_DEAD=[];
Mean_CycleStart_TR_in_green_DEAD=[];
Mean_CycleStart_TR_in_green_MIX=[];   

Mean_CycleStart_TR_in_white_SURVIVORS=[];
Mean_CycleStart_TR_in_green_SURVIVORS=[];
Mean_CycleStart_TR_in_white_MIX=[];   

for i=1:length(Mean_CycleStart_TR_in_white)
    if survival_white(i)==1
        Mean_CycleStart_TR_in_white_SURVIVORS=[Mean_CycleStart_TR_in_white_SURVIVORS,Mean_CycleStart_TR_in_white(i)];
    elseif isnan(survival_white(i))
        Mean_CycleStart_TR_in_white_DEAD=[Mean_CycleStart_TR_in_white_DEAD,Mean_CycleStart_TR_in_white(i)];
    else
        Mean_CycleStart_TR_in_white_MIX=[Mean_CycleStart_TR_in_white_MIX,Mean_CycleStart_TR_in_white(i)];
    end
end


for i=1:length(Mean_CycleStart_TR_in_green)
    if survival_green(i)==1
        Mean_CycleStart_TR_in_green_SURVIVORS=[Mean_CycleStart_TR_in_green_SURVIVORS,Mean_CycleStart_TR_in_green(i)];
    elseif isnan(survival_green(i))
        Mean_CycleStart_TR_in_green_DEAD=[Mean_CycleStart_TR_in_green_DEAD,Mean_CycleStart_TR_in_green(i)];
    else
        Mean_CycleStart_TR_in_green_MIX=[Mean_CycleStart_TR_in_green_MIX,Mean_CycleStart_TR_in_green(i)]; %either one daughter survived or one lost
    end
end


All_green=[Mean_CycleStart_TR_in_green_SURVIVORS,Mean_CycleStart_TR_in_green_DEAD,Mean_CycleStart_TR_in_green_MIX];


[p,h]=ranksum(Mean_CycleStart_TR_in_green_SURVIVORS,Mean_CycleStart_TR_in_green_DEAD)
[p,h]=ranksum(Mean_CycleStart_TR_in_green_DEAD,All_green)
[p,h]=ranksum(Mean_CycleStart_TR_in_green_SURVIVORS,All_green)


All_white=[Mean_CycleStart_TR_in_white_SURVIVORS,Mean_CycleStart_TR_in_white_DEAD,Mean_CycleStart_TR_in_white_MIX];

[p,h]=ranksum(Mean_CycleStart_TR_in_white_SURVIVORS,Mean_CycleStart_TR_in_white_DEAD)
[p,h]=ranksum(Mean_CycleStart_TR_in_white_DEAD,All_white)
[p,h]=ranksum(Mean_CycleStart_TR_in_white_SURVIVORS,All_white)


Numb_Surv=length(Mean_CycleStart_TR_in_white_SURVIVORS)+length(Mean_CycleStart_TR_in_green_SURVIVORS);
Numb_Dead=length(Mean_CycleStart_TR_in_white_DEAD)+length(Mean_CycleStart_TR_in_green_DEAD);
Numb_All=length(All_white)+length(All_green);
Perc_Dead_white=length(Mean_CycleStart_TR_in_white_DEAD)/length(All_white);
Perc_Dead_green=length(Mean_CycleStart_TR_in_green_DEAD)/length(All_green);



