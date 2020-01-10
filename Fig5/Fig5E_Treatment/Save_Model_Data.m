
sync=45;

Start_Gr=0.45;
x_G1=3.6;
x_SG2M=3.4;


if sync==10
    CycleStart_reduced_gr=[];
    CycleStart_reduced_wh=rand(1000,1)*0.1;
elseif sync==45
    CycleStart_reduced_gr=rand(1000,1)*0.1 + Start_Gr;
    CycleStart_reduced_wh=[];
elseif sync==20
    CycleStart_reduced_gr=[];
    number=length(CycleStart_reduced_gr);
    CycleStart_reduced_wh=rand(1000,1)*0.1 +0.1;
elseif sync==30
    CycleStart_reduced_gr=[];
    CycleStart_reduced_wh=rand(1000,1)*0.1 +0.2;
elseif sync==40
    CycleStart_reduced_gr=[];
    CycleStart_reduced_wh=rand(1000,1)*0.1 +0.3;
elseif sync==50
    CycleStart_reduced_gr=rand(1000,1)*0.05 + 0.45;
    CycleStart_reduced_wh=rand(1000,1)*0.05 + 0.4;
elseif sync==60
    CycleStart_reduced_gr=rand(1000,1)*0.1 + 0.5;
    CycleStart_reduced_wh=[];
elseif sync==70
    CycleStart_reduced_gr=rand(1000,1)*0.1 + 0.6;
    CycleStart_reduced_wh=[];
elseif sync==80
    CycleStart_reduced_gr=rand(1000,1)*0.1 + 0.7;
    CycleStart_reduced_wh=[];
elseif sync==90
    CycleStart_reduced_gr=rand(1000,1)*0.1 + 0.8;
    CycleStart_reduced_wh=[];
elseif sync==100
    CycleStart_reduced_gr=rand(1000,1)*0.1 + 0.9;
    CycleStart_reduced_wh=[];
end

number=length(CycleStart_reduced_gr);
number1=length(CycleStart_reduced_wh);

load('Data.mat')



Para_death_white=lognfit(t_TRAIL_Death(Mean_CycleStart_TR_in_white<0.22));
sample_death_white=lognrnd(Para_death_white(1),Para_death_white(2),1,number);
m1=1./sample_death_white;


g_zw=lognrnd(2.1802,0.1915,1,number);
g_zw_NEW=g_zw+x_SG2M;
g_green=(1-Start_Gr)./g_zw_NEW;


%%%%%%%%%%%%%%% TRAIL in G1 %%%%%%%%%%%%%%%%%%%%%%%%%

g_zw=lognrnd(2.1802,0.1915,1,number1);
g_zw_NEW=g_zw+x_SG2M;
g_green_NEW=(1-Start_Gr)./g_zw_NEW; 


g_zw=lognrnd(1.9594,0.2725,1,number1);
g_zw_NEW=g_zw+x_G1;   %prolongation of cell cycle length
g_white=Start_Gr./g_zw_NEW; 

sample_death_white=lognrnd(Para_death_white(1),Para_death_white(2),1,number1);
m2=1./sample_death_white;

CycleStart_reduced_gr=CycleStart_reduced_gr';
CycleStart_reduced_wh=CycleStart_reduced_wh';

if sync==10
    save('Data_Model_Sample_10Per.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
        'g_green_NEW','g_white','number','number1');
elseif sync==45
    save('Data_Model_Sample_S.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
        'g_green_NEW','g_white','number','number1');
elseif sync==20
    save('Data_Model_Sample_20Per.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
        'g_green_NEW','g_white','number','number1');
elseif sync==30
    save('Data_Model_Sample_30Per.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
        'g_green_NEW','g_white','number','number1');
elseif sync==40
    save('Data_Model_Sample_40Per.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
        'g_green_NEW','g_white','number','number1');
elseif sync==50
    save('Data_Model_Sample_50Per.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
        'g_green_NEW','g_white','number','number1')
elseif sync==60
    save('Data_Model_Sample_60Per.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
        'g_green_NEW','g_white','number','number1')
elseif sync==70
    save('Data_Model_Sample_70Per.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
        'g_green_NEW','g_white','number','number1')
elseif sync==80
    save('Data_Model_Sample_80Per.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
        'g_green_NEW','g_white','number','number1')
elseif sync==90
    save('Data_Model_Sample_90Per.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
        'g_green_NEW','g_white','number','number1')
elseif sync==100
    save('Data_Model_Sample_100Per.mat','CycleStart_reduced_gr','CycleStart_reduced_wh','m1','m2','g_green',...
        'g_green_NEW','g_white','number','number1')
end








