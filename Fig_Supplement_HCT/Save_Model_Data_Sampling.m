clear all;

cd ..
run Import_Cycle_Data.m
cd 7_Model_IQM

load('Data.mat')

Sample_Size=3;    %1: large initial number, 2: same initial number as data

C_WE=0.39;  %0.45 NEW
x_G1=2.5;
x_SG2M=1.95;

%%%%%%%%%%%%%%% TRAIL in SG2M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Sample_Size==2
    CycleStart_reduced_gr=repelem(Mean_CycleStart_TR_in_green,1);
elseif Sample_Size==1
    CycleStart_reduced_gr=repelem(Mean_CycleStart_TR_in_green,83);
elseif Sample_Size==3
    CycleStart_reduced_gr=repelem(Mean_CycleStart_TR_in_green,20);
end

number=length(CycleStart_reduced_gr);

%Para_death_white=lognfit(t_TRAIL_Death(Mean_CycleStart_TR_in_white<0.22 & isnan(t_TRAIL_WhiteEnd_F0)))
Para_death_white(1)=1.7315     %1.9654;    %%%%%%%%%%%%% NEW 2019
Para_death_white(2)=0.4053      %0.5618;
sample_death_white=lognrnd(Para_death_white(1),Para_death_white(2),1,number);
m1=1./sample_death_white;

g_zw=lognrnd(log_SG2M(1),log_SG2M(2),1,number);    
g_zw_NEW=g_zw+x_SG2M;
g_green=(1-C_WE)./g_zw_NEW;

g_zw=lognrnd(log_G1(1),log_G1(2),1,number);    
g_zw_NEW=g_zw+x_G1;                           %prolongation of cell cycle length
g_white_NEW=C_WE./g_zw_NEW; 

%(ones(number,1)*mean(G_zw))';  %0.0656 0.0656 0.0655 0.0656  %from MEMO...
%k=lognrnd(6.25,0.14,1,number);  %green f0 tr gr
%G_zw=(1-C_WE)./(k./60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% TRAIL in G1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Sample_Size==2
    CycleStart_reduced_wh=repelem(Mean_CycleStart_TR_in_white,1);
elseif Sample_Size==1
    CycleStart_reduced_wh=repelem(Mean_CycleStart_TR_in_white,83);
elseif Sample_Size==3
    CycleStart_reduced_wh=repelem(Mean_CycleStart_TR_in_white,20);
end

number1=length(CycleStart_reduced_wh);

g_zw=lognrnd(log_SG2M(1),log_SG2M(2),1,number1);
g_zw_NEW=g_zw+x_SG2M;
g_green_NEW=(1-C_WE)./g_zw_NEW; 

g_zw=lognrnd(log_G1(1),log_G1(2),1,number1);    
g_zw_NEW=g_zw+x_G1;                           %prolongation of cell cycle length
g_white=C_WE./g_zw_NEW; 

sample_death_white=lognrnd(Para_death_white(1),Para_death_white(2),1,number1);
m2=1./sample_death_white;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Sample_Size==1
    save('B0_Data_Model_Sample_1.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
    'g_green_NEW','g_white','number','number1','g_white_NEW');   %NEW g_white_NEW
elseif Sample_Size==2
    save('B0_Data_Model_Sample_LESS.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
    'g_green_NEW','g_white','number','number1','g_white_NEW');  %NEW g_white_NEW
elseif Sample_Size==3
    save('B0_Data_Model_Sample_Med_5.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
    'g_green_NEW','g_white','number','number1','g_white_NEW');  %NEW g_white_NEW
end

