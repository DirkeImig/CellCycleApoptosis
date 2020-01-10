run Import_Cycle_Data.m;       %cell cycle data of untreated cells

C_WE_zw = 0.45;   %0.39 for HCT
x_G1=3.6;         %2.5 for HCT
x_SG2M=3.4;       %2 for HCT

%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data_white=Function_Import_Data_white('Dead_Survivor_white.csv');
t_Birth_TRAIL_white=Data_white(1,:);
t_TRAIL_WhiteEnd_F0=Data_white(2,:);
t_TRAIL_Div_F0=Data_white(3,:);
t_TRAIL_WhiteEnd_F1=Data_white(4,:);
t_TRAIL_Death=Data_white(5,:);

Data_green=Function_Import_Data_green('Dead_Survivor_green.csv');
t_WhiteEnd_TRAIL_green=Data_green(1,:);
t_TRAIL_Div_F0_TRinGreen=Data_green(2,:);
t_TRAIL_WhiteEnd_F1_TRinGreen=Data_green(3,:);
t_TRAIL_Div_F1_TRinGreen=Data_green(4,:); 
t_TRAIL_Death_TRinGreen=Data_green(5,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pdflog_green=makedist('Lognormal','mu',log_SG2M(1),'sigma',log_SG2M(2));
pdflog_white=makedist('Lognormal','mu',log_G1(1),'sigma',log_G1(2));

t_white_Zw = lognrnd(log_G1(1), log_G1(2), 1, 1000);   %length white phase of untreated cells
t_green_Zw= lognrnd(log_SG2M(1), log_SG2M(2), 1, 1000);  %length green phase of untreated cells



%%%%%%%%% Calculation of C0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Szenario I: TRAIL in G1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Birth_to_WhiteEnd=[];
for i=1:length(t_Birth_TRAIL_white)
    if  ~isnan(t_Birth_TRAIL_white(i)) && t_TRAIL_WhiteEnd_F0(i)>0
        Birth_to_WhiteEnd(i)=t_Birth_TRAIL_white(i)+t_TRAIL_WhiteEnd_F0(i);
    else
        Birth_to_WhiteEnd(i)=nan;
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
        C_TR_zw = t_Birth_TRAIL_white(i)*C_WE_zw./t_white_rand;  
        Std_CycleStart_TR_in_white=[Std_CycleStart_TR_in_white, std(C_TR_zw)];
        Mean_CycleStart_TR_in_white=[Mean_CycleStart_TR_in_white, mean(C_TR_zw)];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%% Szenario II: TRAIL in SG2M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LengthGreen=[];

for i=1:length(t_WhiteEnd_TRAIL_green)
    if ~isnan(t_WhiteEnd_TRAIL_green(i)) && t_TRAIL_Div_F0_TRinGreen(i)>0
        LengthGreen(i)=t_WhiteEnd_TRAIL_green(i)+t_TRAIL_Div_F0_TRinGreen(i);
    else
        LengthGreen(i)=nan;
    end
end

CycleStart_TR_in_green=[];  
Mean_CycleStart_TR_in_green=[];
Std_CycleStart_TR_in_green=[];

for i = 1:length(t_WhiteEnd_TRAIL_green)
    if t_WhiteEnd_TRAIL_green(i)==0
        Std_CycleStart_TR_in_green(i)=0;
        Mean_CycleStart_TR_in_green(i)=C_WE_zw;
    elseif LengthGreen(i)>0
        p=-((1-C_WE_zw)/t_WhiteEnd_TRAIL_green(i) - (1-C_WE_zw)*LengthGreen(i)/(x_SG2M*t_WhiteEnd_TRAIL_green(i)));
        q=-((1-C_WE_zw)^2)/(x_SG2M*t_WhiteEnd_TRAIL_green(i));
        g1=-p/2+sqrt((p/2)^2 -q);
        C_TR_zw=C_WE_zw + t_WhiteEnd_TRAIL_green(i)*g1;  
        Mean_CycleStart_TR_in_green(i)=mean(C_TR_zw);
        Std_CycleStart_TR_in_green(i)=std(C_TR_zw);
    else
        Trunc_gen_t_green=truncate(pdflog_green,t_WhiteEnd_TRAIL_green(i),Inf);
        t_green_rand = random(Trunc_gen_t_green,1000,1);
        C_TR_zw = C_WE_zw + (t_WhiteEnd_TRAIL_green(i) * (1-C_WE_zw)./(t_green_rand) );  
        Std_CycleStart_TR_in_green(i)=std(C_TR_zw);
        Mean_CycleStart_TR_in_green(i)=mean(C_TR_zw);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('Data.mat','t_TRAIL_Death','t_TRAIL_Death_TRinGreen', ...
    'Mean_CycleStart_TR_in_white','Mean_CycleStart_TR_in_green');

