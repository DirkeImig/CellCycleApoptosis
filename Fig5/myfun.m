function [y,Cyc_gr,Cyc_wh,Td_gr,Td_wh, d_ph_gr, d_ph_wh] = myfun(p,PAS,Sample)

C_WE=0.45;

load('Data.mat')


%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CycleStart_gr= Mean_CycleStart_TR_in_green;
CycleStart_wh= Mean_CycleStart_TR_in_white;
T_death_gr = t_TRAIL_Death_TRinGreen;
T_death_wh = t_TRAIL_Death;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    load(Sample)
catch
    cd .. 
    cd Fig5E_Treatment
    load(Sample)
    cd ..
    cd Fig5
end


%%%%%%%%%%%%%%% TRAIL in SG2M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_Death=[];
d_phase_gr=[];
t_PAS=(PAS-CycleStart_reduced_gr)./g_green;
%scenario 1: TRAIL past PAS, m1>p, death after division (f1).
%scenario 2: TRAIL past PAS, m1<=p, death after division (f1).
%sceanrio 3: TRAIL past PAS, death in f0
%scenario 4: TRAIL before PAS, important: Cd2 < PAS: death before PAS
%scenario 5: TRAIL before PAS, death after PAS in f0, m(i)>p
%scenario 6: TRAIL before PAS, death in f1, m(i)>p
%sceanrio 7: TRAIL before PAS, death in f1, m(i)<=p

A_PAS=m1.*t_PAS;
t_Div=(1-CycleStart_reduced_gr)./g_green;             %division time (scenario 1,2 6)
td1_SG2M=1./(m1-p);                                %death time (TRAIL past PAS, death past PAS)  
td2_SG2M=1./m1;                                    %death time (TRAIL before PAS, death before PAS)
td3_SG2M=(1-A_PAS + ((m1-p).*t_PAS))./(m1-p);      %death time (TRAIL before PAS, death past PAS) (scenario 5)
Cd1_SG2M=CycleStart_reduced_gr + (g_green.*td1_SG2M); %cycle at death
Cd2_SG2M=CycleStart_reduced_gr + (g_green.*td2_SG2M); %cycle at death (TRAIL before PAS, death before PAS) (scenario 4)
Cd3_SG2M=CycleStart_reduced_gr + (g_green.*td3_SG2M); %cycle at death (scenario 5)

for i=1:number
    if CycleStart_reduced_gr(i) >= PAS && m1(i)>p  &&  Cd1_SG2M(i)>1                                           %scenario 1: TRAIL past PAS, m1>p, death after division (f1).
        A_Div=(m1(i)-p).*t_Div(i);                                                                             %apoptosis at division (A_Div<1)
        T_Death=[T_Death, (1-A_Div+(m1(i)*t_Div(i)))/m1(i)];                                                   %death time
        d_phase_gr=[d_phase_gr,3];
    elseif CycleStart_reduced_gr(i) >= PAS && p>=m1(i)                                                         %scenario 2: TRAIL past PAS, m1<=p, death after division (f1).
        T_Death=[T_Death, (1+(m1(i)*t_Div(i)))/m1(i)];                                                         %assuming that A_Div non-negativ!
        d_phase_gr=[d_phase_gr,3];
    elseif CycleStart_reduced_gr(i) >= PAS  &&  Cd1_SG2M(i)<=1  && m1(i)>p                                     %sceanrio 3: TRAIL past PAS, death in f0
        T_Death=[T_Death, (1/(m1(i)-p))];
        d_phase_gr=[d_phase_gr,2];
    elseif CycleStart_reduced_gr(i) < PAS && Cd2_SG2M(i)<=PAS                                                     %scenario 4: TRAIL before PAS, important: Cd2 < PAS: death before PAS
        T_Death=[T_Death, 1/m1(i)];
        d_phase_gr=[d_phase_gr,2];
    elseif CycleStart_reduced_gr(i) < PAS && Cd2_SG2M(i)>PAS && Cd3_SG2M(i)<=1 && td3_SG2M(i)>t_PAS(i) && m1(i)>p %scenario 5: TRAIL before PAS, death after PAS in f0, m(i)>p
        T_Death=[T_Death, td3_SG2M(i)];
        d_phase_gr=[d_phase_gr,2];
    elseif CycleStart_reduced_gr(i) < PAS && Cd2_SG2M(i)>PAS && Cd3_SG2M(i)>1 && td3_SG2M(i)>t_PAS(i) && m1(i)>p  %scenario 6: TRAIL before PAS, death in f1, m(i)>p
        A_Div=(m1(i)-p)*t_Div(i) - (m1(i)-p)*t_PAS(i) + A_PAS(i);
        T_Death=[T_Death, (1-A_Div+(m1(i)*t_Div(i)))/m1(i)];
        d_phase_gr=[d_phase_gr,3];
    elseif CycleStart_reduced_gr(i) < PAS && Cd2_SG2M(i)>PAS && p>=m1(i)                                          %sceanrio 7: TRAIL before PAS, death in f1, m(i)<=p
        A_Div=(m1(i)-p)*t_Div(i) - (m1(i)-p)*t_PAS(i) + A_PAS(i);
        T_Death=[T_Death, (1-max(A_Div,0)+(m1(i)*t_Div(i)))/m1(i)];                                               %assuming that A_Div non-negativ!
        d_phase_gr=[d_phase_gr,3];
    else
        i
    end
end


%%%%%%%%%%%%%%% TRAIL in G1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_Death_wh=[];
d_phase_wh=[];
t_CWe=(C_WE-CycleStart_reduced_wh)./g_white;                 %time of white end (for all)
t_DIV=(1-C_WE + (t_CWe.*g_green_NEW))./g_green_NEW;          %tdiv for all

%A: PAS in SG2M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scenario 1: death before PAS in G1
%scenario 2: death in SG2M, before PAS
%scenario 3: death in SG2M, f0, past PAS, m(i)>p
%scenario 4: death in f1, past PAS, m(i)>p
%scenario 5: death in f1, past PAS, m(i)<=p 
t_PAS_1= (PAS-C_WE + (g_green_NEW.*t_CWe))./g_green_NEW;     %time of PAS (for PAS in SG2M)
A_PAS_1= m2.*t_PAS_1;                                        %apoptosis at PAS (for TRAIL before PAS)
td1=1./m2;                                                   %death in G1, f0 (scenario 1)
Cd1=CycleStart_reduced_wh + g_white.*td1;                    %cycle at death (scenario 1)
Cd2=C_WE + td1.*g_green_NEW - t_CWe.*g_green_NEW;            %cycle at death, death before PAS in SG2M (scenario 2)
td3=(1-A_PAS_1 + t_PAS_1.*(m2-p))./(m2-p);                   %time of death, TRAIL in G1, death past PAS in f0, SG2M (scenario 3)
Cd3=C_WE + td3.*g_green_NEW - t_CWe.*g_green_NEW;            
A_Div_4= (t_DIV -t_PAS_1).*(m2-p) + A_PAS_1;                 %apoptosis at division (scenario 4,5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%% B: PAS in G1!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% B1: TRAIL before PAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scenario 6: death in G1, before PAS
%scenario 7: death in G1, past PAS 
%scenario 8: death past PAS in SG2M, 
%scenario 9: death past PAS in f1, m(i)<=p 
%scenario 10: death past PAS in f1, m(i)>p 
t_PAS_7=(PAS-CycleStart_reduced_wh)./g_white;           %time of PAS, PAS in G1 (only for PAS>CycleStart_reduced_wh)
A_PAS_7= m2.*t_PAS_7;                                   %apoptosis at PAS (only for PAS>CycleStart_reduced_wh)
td7= (1- A_PAS_7 + (t_PAS_7.*(m2-p)))./(m2-p);          %time of death, m>p (only for PAS>CycleStart_reduced_wh)
Cd7=CycleStart_reduced_wh + td7.*g_white;               %cycle at time of death, death after PAS in G1 (scenario 7)
A_WE_8= A_PAS_7 + (m2-p).*(t_CWe-t_PAS_7);              %apoptotis white end
td8 =  td7;                                             %time of death for death in SG2M, PAS in G1, TRAIL before PAS
Cd8 = (td8 - t_CWe).*g_green_NEW + C_WE;                %cycle of death 
A_DIV_9=(t_DIV-t_PAS_7).*(m2-p) + A_PAS_7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% B2: TRAIL past PAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scenario 11: death past PAS in G1, m(i)>p
%scenario 12: death in SG2M 
%scenario 13: death in f1, m(i)<=p 
%scenario 14: death in f1, m(i)>p
td11=1./(m2-p);                                         %for death in f0 (if TRAIL past PAS)
Cd11=CycleStart_reduced_wh + g_white.*td11;             %cycle at death, for death in G1 (scenario 11)
td12=td11;                                              %death time (if TRAIL past PAS and death in SG2M) (scenario 12)
Cd12=C_WE + ((td12-t_CWe).*g_green_NEW);                %cycle at death time (scenario 12)
td13=(1+t_DIV.*m2)./m2;                                 %time death if m<p (scenario 13)
A_DIV_14=t_DIV.*(m2-p);                                 %A_Div if m>p  (scenario 14)
td14=(1-A_DIV_14 + (t_DIV.*m2))./m2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i=1:number1
    if PAS>=C_WE                                                                            %A: PAS in SG2M
        if Cd1(i)<=C_WE                                                                     %scenario 1: death before PAS in G1
            T_Death_wh=[T_Death_wh, td1(i)];
            d_phase_wh=[d_phase_wh,1];
        elseif Cd1(i)>C_WE && Cd2(i) <= PAS                                                 %scenario 2: death in SG2M, before PAS
            T_Death_wh=[T_Death_wh, td1(i)];
            d_phase_wh=[d_phase_wh,2];
        elseif Cd2(i) > PAS && Cd3(i) <= 1 && td3(i) > t_PAS_1(i) && m2(i)>p                %scenario 3: death in SG2M, f0, past PAS, m(i)>p
            T_Death_wh=[T_Death_wh, td3(i)];
            d_phase_wh=[d_phase_wh,2];
        elseif Cd2(i) > PAS && Cd3(i) > 1 && td3(i) > t_PAS_1(i) && m2(i)>p                 %scenario 4: death in f1, past PAS, m(i)>p
            T_Death_wh=[T_Death_wh, (1- A_Div_4(i) +(m2(i)*t_DIV(i)))/m2(i)];
            d_phase_wh=[d_phase_wh,3];
        elseif Cd2(i) > PAS && p>=m2(i)                                                     %scenario 5: death in f1, past PAS, m(i)<=p 
            T_Death_wh=[T_Death_wh, (1- max(A_Div_4(i),0) +(m2(i)*t_DIV(i)))/m2(i)];
            d_phase_wh=[d_phase_wh,3];
        else
            i
        end
    else                                                                                         %B: PAS in G1
        if CycleStart_reduced_wh(i) < PAS                                                        %B1 TRAIL before PAS
           if Cd1(i)<= PAS                                                                       %scenario 6: death in G1, before PAS
                T_Death_wh=[T_Death_wh, td1(i)];
                d_phase_wh=[d_phase_wh,1];
            elseif Cd7(i)<=C_WE && Cd1(i)>PAS && m2(i)>p                                         %scenario 7: death in G1, past PAS
                T_Death_wh=[T_Death_wh, td7(i)];
                d_phase_wh=[d_phase_wh,1];
            elseif Cd7(i) > C_WE && Cd8(i)<=1 && m2(i)>p                                         %scenario 8: death past PAS in SG2M
                T_Death_wh=[T_Death_wh, td8(i)];
                d_phase_wh=[d_phase_wh,2];
           elseif Cd1(i)>PAS && m2(i)<=p                                                        %scenario 9: death past PAS in f1, m(i)<=p 
                T_Death_wh=[T_Death_wh, (1- max(A_DIV_9(i),0) + (m2(i)*t_DIV(i)))/m2(i)];
                d_phase_wh=[d_phase_wh,3];
           elseif Cd8(i)>1 && Cd1(i)>PAS &&  m2(i)>p                                            %scenario 10: death past PAS in f1, m(i)>p 
                T_Death_wh=[T_Death_wh, (1- A_DIV_9(i) + (m2(i)*t_DIV(i)))/m2(i)];
                d_phase_wh=[d_phase_wh,3];
           else
                i
           end
        else                                                                                     %B2 TRAIL past PAS
            if Cd11(i)<=C_WE &&  m2(i)>p                                                         %scenario 11: death in G1, m(i)>p
                T_Death_wh=[T_Death_wh, td11(i)];
                d_phase_wh=[d_phase_wh,1];
            elseif  Cd11(i)>C_WE && Cd12(i)<=1 && m2(i)>p                                        %scenario 12: death in SG2M 
                T_Death_wh=[T_Death_wh, td12(i)];
                d_phase_wh=[d_phase_wh,2];
            elseif m2(i)<=p                                                                      %scenario 13: death in f1, m(i)<=p 
                T_Death_wh=[T_Death_wh, td13(i)];
                d_phase_wh=[d_phase_wh,3];
            elseif Cd12(i)>1 && m2(i)>p                                                          %scenario 14: death in f1, m(i)>p
                T_Death_wh=[T_Death_wh, td14(i)];           
                d_phase_wh=[d_phase_wh,3];
            else
                i
            end
        end
    end
end




%%%%%%%%%%%%%%%%% compare with Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cycle_Start_Model=[CycleStart_reduced_wh, CycleStart_reduced_gr];
Cycle_Start_Model=round(100.*Cycle_Start_Model)./100;
T_Death_Model=[T_Death_wh, T_Death];
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
y=sum(-log(A));


Cyc_gr=CycleStart_reduced_gr;
Cyc_wh=CycleStart_reduced_wh;
Td_gr=T_Death;
Td_wh=T_Death_wh;
d_ph_gr=d_phase_gr;
d_ph_wh=d_phase_wh;


end



