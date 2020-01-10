clear all;

cd ..
run Import_Cycle_Data.m
load('Data.mat')
cd Fig5


Sample_Size=1;

C_WE=0.45;
x_G1=3.6;
x_SG2M=3.4;

%%%%%%%%%%%%%%% TRAIL in SG2M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Sample_Size==2
    CycleStart_reduced_gr=Mean_CycleStart_TR_in_green;
elseif Sample_Size==1
    CycleStart_reduced_gr=repelem(Mean_CycleStart_TR_in_green,83);
end

number=length(CycleStart_reduced_gr);

Para_death_white=lognfit(t_TRAIL_Death(Mean_CycleStart_TR_in_white<0.22)); %death time in early white cells
sample_death_white=lognrnd(Para_death_white(1),Para_death_white(2),1,number);
m1=1./sample_death_white;

g_zw=lognrnd(log_SG2M(1),log_SG2M(2),1,number);    
g_zw_NEW=g_zw+x_SG2M;
g_green=(1-C_WE)./g_zw_NEW;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% TRAIL in G1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Sample_Size==2
    CycleStart_reduced_wh=Mean_CycleStart_TR_in_white;
elseif Sample_Size==1
    CycleStart_reduced_wh=repelem(Mean_CycleStart_TR_in_white,83);
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
    save('Data_Model_Sample_0.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
    'g_green_NEW','g_white','number','number1');
elseif Sample_Size==2
    save('Data_Model_Sample_LESS.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
    'g_green_NEW','g_white','number','number1');
end

