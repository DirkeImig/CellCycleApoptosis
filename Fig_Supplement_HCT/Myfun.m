function [y,Cyc_gr,Cyc_wh,Td_gr,Td_wh, d_ph_gr, d_ph_wh] = Myfun(p,PAS,Sample,MyModel)

load('Data.mat')

%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CycleStart_gr= Mean_CycleStart_TR_in_green;
CycleStart_wh= Mean_CycleStart_TR_in_white;
T_death_gr = t_TRAIL_Death_TRinGreen;
T_death_wh = t_TRAIL_Death;  %this is TRinWhite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tmax=60;
step=300;
t=linspace(0,Tmax,step);
method='ode23s';
death_time_wh=[];
death_phase_wh=[];
death_time_gr=[];
death_phase_gr=[];
Cycle_Start_G1=[];
Cycle_Start_SG2M=[];
IC=[0,0];
paramNames = {'m1','p','G','g1','g2','PAD','p1'};
paramvec = [0,0,0,0,0,0,0];
OPTIONS.method='nonstiff';


try
    load(Sample)
catch
    cd .. 
    cd 8_Treatment
    load(Sample)
    cd ..
    cd 7_Model
end



for i=1:length(CycleStart_reduced_wh)
    IC=[0,CycleStart_reduced_wh(i)];
    if PAS<CycleStart_reduced_wh(i)
        p_zw = p;
    else
        p_zw = 0;
    end
    paramvec = [m2(i),p_zw,g_white(i),g_white(i),g_green_NEW(i),PAS,p];
    output = IQMPsimulate(MyModel,t,IC,paramNames,paramvec,OPTIONS); 
    for j=1:length(t)
        if output.statevalues(j,1)>=1
            death_time_wh=[death_time_wh,output.time(j)];
            death_phase_wh=[death_phase_wh,output.statevalues(j,2)];
            Cycle_Start_G1=[Cycle_Start_G1,CycleStart_reduced_wh(i)];
        break
    end
end
end


for i=1:length(CycleStart_reduced_gr)
    IC=[0,CycleStart_reduced_gr(i)];
    if PAS<CycleStart_reduced_gr(i)
        p_zw = p;
    else
        p_zw = 0;
    end
    paramvec = [m1(i),p_zw,g_green(i),g_white_NEW(i),g_green(i),PAS,p];
    output = IQMPsimulate(MyModel,t,IC,paramNames,paramvec,OPTIONS); 
    for j=1:length(t)
        if output.statevalues(j,1)>=1
            death_time_gr=[death_time_gr,output.time(j)];
            death_phase_gr=[death_phase_gr,output.statevalues(j,2)];
            Cycle_Start_SG2M=[Cycle_Start_SG2M,CycleStart_reduced_gr(i)];
        break
    end
end
end

d_phase_wh=zeros(1,length(death_phase_wh));
d_phase_wh(death_phase_wh<0.39)=1;
d_phase_wh(death_phase_wh<1 & death_phase_wh>=0.39)=2;
d_phase_wh(death_phase_wh<1.39 & death_phase_wh>=1)=3;
d_phase_wh(death_phase_wh<2 & death_phase_wh>=1.39)=4;
d_phase_wh(death_phase_wh<2.39 & death_phase_wh>=2)=5;
%irgrendwas stimmt da nicht..

d_phase_gr=zeros(1,length(death_phase_gr));
d_phase_gr(death_phase_gr<0.39)=1;
d_phase_gr(death_phase_gr<1 & death_phase_gr>=0.39)=2;
d_phase_gr(death_phase_gr<1.39 & death_phase_gr>=1)=3;
d_phase_gr(death_phase_gr<2 & death_phase_gr>=1.39)=4;
d_phase_gr(death_phase_gr<2.39 & death_phase_gr>=2)=5;
%else...

%%%%%%%%%%%%%%%%% compare with Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cycle_Start_Model=[Cycle_Start_G1, Cycle_Start_SG2M];
Cycle_Start_Model=round(100.*Cycle_Start_Model)./100;
T_Death_Model=[death_time_wh, death_time_gr];
d_phase=[d_phase_wh, d_phase_gr];

Cycle_Start_Model_NEW=[Cycle_Start_Model,Cycle_Start_Model(d_phase>2)];    %division 
T_Death_Model_NEW=[T_Death_Model,T_Death_Model(d_phase>2)];
d_phase_NEW=[d_phase,d_phase(d_phase>2)];

%%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CycleStart_Data=[CycleStart_wh,CycleStart_gr];
CycleStart_Data=round(100.*CycleStart_Data)./100;
T_Death_Data=[T_death_wh,T_death_gr];


Model_Points=[Cycle_Start_Model_NEW;T_Death_Model_NEW]';
Data_Points=[CycleStart_Data;T_Death_Data]';
A=ksdensity(Model_Points,Data_Points);
A(A==0)=1e-9;
y=sum(-log(A))


Cyc_gr=Cycle_Start_SG2M;
Cyc_wh=Cycle_Start_G1;
Td_gr=death_time_gr;
Td_wh=death_time_wh;
d_ph_gr=d_phase_gr;
d_ph_wh=d_phase_wh;


end




