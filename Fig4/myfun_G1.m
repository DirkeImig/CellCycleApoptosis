function y = myfun_G1(z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=500;                  %data points in cycle
n2=1000;                %for normalization
prolong=z;              %'real' prolongation induced by TRAIL

%%%%%%%%%%%%%%%%%%% data unstimulated and stimulated %%%%%%%%%%%%%%%%%%%%%%
cd ..
run Import_Cycle_Data.m;       %cell cycle data of untreated cells
run Mean_Cycle_Start.m;        %for t_death vlaues (below)
cd Fig4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pdflog_green=makedist('Lognormal','mu',log_SG2M(1),'sigma',log_SG2M(2));
pdflog_white=makedist('Lognormal','mu',log_G1(1),'sigma',log_G1(2));
C_WE_zw=0.45;

%t_TRAIL_Death_white_gen0=t_TRAIL_Death(~isnan(t_Birth_TRAIL_white) & t_TRAIL_Div_F0==0);
%t_TRAIL_Death_white_gen1=t_TRAIL_Death(~isnan(t_Birth_TRAIL_white) & t_TRAIL_Div_F0>0);
%t_TRAIL_Death_green_gen0=t_TRAIL_Death_TRinGreen(~isnan(t_WhiteEnd_TRAIL_green) & t_TRAIL_Div_F0_TRinGreen==0);
%t_TRAIL_Death_green_gen1=t_TRAIL_Death_TRinGreen(~isnan(t_WhiteEnd_TRAIL_green) & t_TRAIL_Div_F0_TRinGreen>0);
%t_TRAIL_Death_All_weighted=[t_TRAIL_Death_white_gen0,t_TRAIL_Death_white_gen0,t_TRAIL_Death_green_gen0,t_TRAIL_Death_green_gen0,t_TRAIL_Death_white_gen1,t_TRAIL_Death_green_gen1];
%paracycle3=lognfit(t_TRAIL_Death_All_weighted);

paracycle3(1)=1.1062;   %t_death values for all 
paracycle3(2)=0.5649;



%%%%%%%%%%%%%%%%% Sampling for all scenarios %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_zw=lognrnd(log_G1(1),log_G1(2),n,1);
g=C_WE_zw./g_zw;
g_zw_new=g_zw+prolong;
g_new=C_WE_zw./g_zw_new;
C_WE_zw./g_new - C_WE_zw./g;
C_0=rand(n,1).*0.4499;
T_death = lognrnd(paracycle3(1),paracycle3(2), n,1);
C_end= C_0 + g_new.*T_death;      
t_Birth_to_TR=C_0./g;         
t_Birth_to_end_G1_unchanged = C_WE_zw./g;
t_Birth_to_end_G1_CHANGED = t_Birth_to_TR + (C_WE_zw-C_0)./g_new;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%% No Truncation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_mean_C0=zeros(n,1);
A_C0_Real=zeros(n,1);

for i = 1:length(t_Birth_to_TR)  
    if t_Birth_to_TR(i)==0
        A_mean_C0(i)=0;
        A_C0_Real(i)=C_0(i);
    else
        t_white_rand = random(pdflog_white,n2,1);
        C_TR_zw = t_Birth_to_TR(i)*C_WE_zw./t_white_rand;  
        A_mean_C0(i)=mean(C_TR_zw);
        A_C0_Real(i)=C_0(i);
    end
end
     
compareA=sum((A_mean_C0-A_C0_Real).^2);
compareA2=sum((A_mean_C0-A_C0_Real));

       


%%%%%%%%%%%%% Truncation with trail addition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B_Mean_C0=zeros(n,1);
B_C0_Real=zeros(n,1);

for i = 1:length(t_Birth_to_TR)
    if t_Birth_to_TR(i)==0
        B_Mean_C0(i)=0;
        B_C0_Real(i)=C_0(i);
    else
        Trunc_gen_t_white=truncate(pdflog_white,t_Birth_to_TR(i),Inf);  
        t_white_rand = random(Trunc_gen_t_white,n2,1);
        C_TR_zw = t_Birth_to_TR(i)*C_WE_zw./t_white_rand;              
        B_Mean_C0(i)=mean(C_TR_zw);
        B_C0_Real(i)=C_0(i);
    end
end

compareB=sum((B_Mean_C0-B_C0_Real).^2);
compareB2=sum((B_Mean_C0-B_C0_Real));




%%%%%%%%%%%%%% TRAIL changes phases length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_Mean_C0=zeros(n,1);
C_C0_Real=zeros(n,1);
x=3.6;                  %estimated prolongation
%%%%%%%%%%%%% Truncation with trail addition and tdiv information %%%%%%%%%
for i = 1:length(t_Birth_to_TR) 
    if t_Birth_to_TR(i)==0
        C_Mean_C0(i)=0;
        C_C0_Real(i)=C_0(i);
    elseif C_end(i) <= C_WE_zw    %death before end of G1 phase
        Trunc_gen_t_white=truncate(pdflog_white,t_Birth_to_TR(i),35);  
        t_white_rand = random(Trunc_gen_t_white,n2,1);
        C_TR_zw = t_Birth_to_TR(i)*C_WE_zw./t_white_rand;              
        C_Mean_C0(i)=mean(C_TR_zw);
        C_C0_Real(i)=C_0(i);
    else
        p=-(C_WE_zw/t_Birth_to_TR(i) - t_Birth_to_end_G1_CHANGED(i)*C_WE_zw/(x*t_Birth_to_TR(i)));
        q=-(C_WE_zw^2)/(x*t_Birth_to_TR(i));
        g1=-p/2+sqrt((p/2)^2 -q);
        C_TR_zw=t_Birth_to_TR(i)*g1;   
        C_Mean_C0(i)=mean(C_TR_zw);
        C_C0_Real(i)=C_0(i);
    end
end


compareC=sum((C_Mean_C0-C_C0_Real).^2);
compareC2=sum((C_Mean_C0-C_C0_Real));



y=[compareA,compareA2,compareB, compareB2, compareC,compareC2];

end





