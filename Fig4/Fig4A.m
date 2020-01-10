%%%%%%%%%%%%% data normalization example simulations %%%%%%%%%%%%%%%%%%%%%%

prolong1=3.6;    %prolongation of cell cycle in hours
prolong2=3.4; 

n=500;   %data points in cycle
n2=1000; %for normalization

paracycle3(1)=1.1062;   %t_death values for all
paracycle3(2)=0.5649;

%%%%%%%%%%%%%%%%%%% data unstimulated and stimulated %%%%%%%%%%%%%%%%%%%%%%
cd ..
run Import_Cycle_Data.m;       %cell cycle data of untreated cells
cd Fig4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pdflog_green=makedist('Lognormal','mu',log_SG2M(1),'sigma',log_SG2M(2));
pdflog_white=makedist('Lognormal','mu',log_G1(1),'sigma',log_G1(2));
C_WE_zw=0.45;



%%%%%%%%%%%%%%%%% Sampling for all scenarios G1 %%%%%%%%%%%%%%%%%%%%%%%%%%%
g_zw=lognrnd(log_G1(1),log_G1(2),n,1);
g=C_WE_zw./g_zw;
g_zw_new=g_zw+prolong1;
g_new=C_WE_zw./g_zw_new;
C_0=rand(n,1).*0.4499;
T_death = lognrnd(paracycle3(1),paracycle3(2), n,1);
C_end= C_0 + g_new.*T_death;      
t_Birth_to_TR=C_0./g;         
t_Birth_to_end_G1_unchanged = C_WE_zw./g;
t_Birth_to_end_G1_CHANGED = t_Birth_to_TR + (C_WE_zw-C_0)./g_new;

%%%%%%%%%%%%%%%%% Sampling for all scenarios SG2M %%%%%%%%%%%%%%%%%%%%%%%%%
g_zw=lognrnd(log_SG2M(1),log_SG2M(2),n,1);
g=(1-C_WE_zw)./g_zw;
g_zw_new=g_zw+prolong2;
g_new=(1-C_WE_zw)./g_zw_new;
C_0_green=rand(n,1).*0.55+C_WE_zw;
T_death = lognrnd(paracycle3(1),paracycle3(2), n,1);
C_end_green= C_0_green + g_new.*T_death;      
t_StartSG2M_to_TR=(C_0_green-C_WE_zw)./g;        
t_StartSG2M_to_end_SG2M_unchanged = (1-C_WE_zw)./g;
t_StartSG2M_to_end_SG2M_CHANGED = t_StartSG2M_to_TR + (1-C_0_green)./g_new;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%% No Truncation G1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%% No Truncation S/G2/M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_mean_C0_green=zeros(n,1);
A_C0_Real_green=zeros(n,1);

for i = 1:length(t_StartSG2M_to_TR)    
    if t_StartSG2M_to_TR(i)==0
        A_mean_C0_green(i)=0;
        A_C0_Real_green(i)=C_0(i);
    else
        t_green_rand = random(pdflog_green,n2,1);
        C_TR_zw = (t_StartSG2M_to_TR(i)./t_green_rand).*0.55 + C_WE_zw;   
        A_mean_C0_green(i)=mean(C_TR_zw);
        A_C0_Real_green(i)=C_0_green(i);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure()
for i=1:n
    scatter(A_mean_C0(i),A_C0_Real(i),40,[0,0,0],'square','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    hold on;
end
for i=1:n
    scatter(A_mean_C0_green(i),A_C0_Real_green(i),40,[0,0.4,0],'square','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    hold on;
end
plot([0,1],[0,1],'--k','LineWidth',3)
hold on;
plot([0,1.3],[C_WE_zw,C_WE_zw],'k--')
hold on;
plot([C_WE_zw,C_WE_zw],[0,1],'k--')
xlim([0,1.3])
ylim([0,1])
xlabel('$\hat{C}_0$','Interpreter','latex')  
ylabel('$C_0$','Interpreter','latex')   
set(gca,'FontSize',25)    
box on;
grid on;





%%%%%%%%%%%%% Truncation with trail addition G1 %%%%%%%%%%%%%%%%%%%%%%%%%%%
B_mean_C0=zeros(n,1);
B_C0_Real=zeros(n,1);

for i = 1:length(t_Birth_to_TR)
    if t_Birth_to_TR(i)==0
        B_mean_C0(i)=0;
        B_C0_Real(i)=C_0(i);
    else
        Trunc_gen_t_white=truncate(pdflog_white,t_Birth_to_TR(i),Inf);  
        t_white_rand = random(Trunc_gen_t_white,n2,1);
        C_TR_zw = t_Birth_to_TR(i)*C_WE_zw./t_white_rand;              
        B_mean_C0(i)=mean(C_TR_zw);
        B_C0_Real(i)=C_0(i);
    end
end

%%%%%%%%%%%%% Truncation with trail addition S/G2/M %%%%%%%%%%%%%%%%%%%%%%%
B_mean_C0_green=zeros(n,1);
B_C0_Real_green=zeros(n,1);

for i = 1:length(t_StartSG2M_to_TR)  
    if t_StartSG2M_to_TR(i)==0
        B_mean_C0_green(i)=0;
        B_C0_Real_green(i)=C_0(i);
    else
        Trunc_gen_t_green=truncate(pdflog_green,t_StartSG2M_to_TR(i),Inf);  
        t_green_rand = random(Trunc_gen_t_green,n2,1);
        C_TR_zw = (t_StartSG2M_to_TR(i)./t_green_rand).*0.55 + C_WE_zw;            
        B_mean_C0_green(i)=mean(C_TR_zw);
        B_C0_Real_green(i)=C_0_green(i);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
for i=1:n
    scatter(B_mean_C0(i),B_C0_Real(i),40,[0,0,0],'square','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    hold on;
end
for i=1:n
    scatter(B_mean_C0_green(i),B_C0_Real_green(i),40,[0,0.4,0],'square','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    hold on;
end
plot([0,1],[0,1],'--k','LineWidth',3)
hold on;
plot([0,1.3],[C_WE_zw,C_WE_zw],'k--')
hold on;
plot([C_WE_zw,C_WE_zw],[0,1],'k--')
xlim([0,1.3])
xlabel('$\hat{C}_0$','Interpreter','latex')  
ylabel('$C_0$','Interpreter','latex')  
set(gca,'FontSize',25)   
box on;
grid on;




%%%%%%%%%%%%%% TRAIL changes phases length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_Mean_C0=zeros(n,1);
C_C0_Real=zeros(n,1);
x1=3.6;
x2=3.4;
%%%%%%%%%%%%% Truncation with trail addition and tdiv information G1 %%%%%%
for i = 1:length(t_Birth_to_TR) 
    if t_Birth_to_TR(i)==0
        C_Mean_C0(i)=0;
        C_C0_Real(i)=C_0(i);
    elseif C_end(i) <= C_WE_zw
        Trunc_gen_t_white=truncate(pdflog_white,t_Birth_to_TR(i),35);  
        t_white_rand = random(Trunc_gen_t_white,n2,1);
        C_TR_zw = t_Birth_to_TR(i)*C_WE_zw./t_white_rand;              
        C_Mean_C0(i)=mean(C_TR_zw);
        C_C0_Real(i)=C_0(i);
    else
        p=-(C_WE_zw/t_Birth_to_TR(i) - t_Birth_to_end_G1_CHANGED(i)*C_WE_zw/(x1*t_Birth_to_TR(i)));
        q=-(C_WE_zw^2)/(x1*t_Birth_to_TR(i));
        g1=-p/2+sqrt((p/2)^2 -q);
        C_TR_zw=t_Birth_to_TR(i)*g1;   
        C_Mean_C0(i)=mean(C_TR_zw);
        C_C0_Real(i)=C_0(i);
    end
end




%%%%%%%%%%%%%% TRAIL changes phases length S/G2/M %%%%%%%%%%%%%%%%%%%%%%%%%
C_Mean_C0_green=zeros(n,1);
C_C0_Real_green=zeros(n,1);
x2=3.4;
%%%%%%%%%%%%% Truncation with trail addition and tdiv information
for i = 1:length(t_StartSG2M_to_TR)  
    if t_StartSG2M_to_TR(i)==0
        C_Mean_C0_green(i)=0;
        C_C0_Real_green(i)=C_0(i);
    elseif C_end_green(i) <= 1
        Trunc_gen_t_green=truncate(pdflog_green,t_StartSG2M_to_TR(i),Inf);  
        t_green_rand = random(Trunc_gen_t_green,n2,1);
        C_TR_zw = (t_StartSG2M_to_TR(i)./t_green_rand).*0.55 + C_WE_zw;            
        C_Mean_C0_green(i)=mean(C_TR_zw);
        C_C0_Real_green(i)=C_0_green(i);
    else
        p=-(0.55/t_StartSG2M_to_TR(i) - 0.55*t_StartSG2M_to_end_SG2M_CHANGED(i)/(x2*t_StartSG2M_to_TR(i)));
        q=-(0.55^2)/(x2*t_StartSG2M_to_TR(i));
        g1=-p/2+sqrt((p/2)^2 -q);
        C_TR_zw=0.45 + t_StartSG2M_to_TR(i)*g1;   
        C_Mean_C0_green(i)=mean(C_TR_zw);
        C_C0_Real_green(i)=C_0_green(i);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




figure()
for i=1:n
    scatter(C_Mean_C0(i),C_C0_Real(i),40,[0,0,0],'square','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    hold on;
end
for i=1:n
    scatter(C_Mean_C0_green(i),C_C0_Real_green(i),40,[0,0.4,0],'square','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    hold on;
end
plot([0,1],[0,1],'--k','LineWidth',3)
hold on;
plot([0,1.3],[C_WE_zw,C_WE_zw],'k--')
hold on;
plot([C_WE_zw,C_WE_zw],[0,1],'k--')
xlim([0,1.3])
xlabel('$\hat{C}_0$','Interpreter','latex')   
ylabel('$C_0$','Interpreter','latex') 
grid on;









