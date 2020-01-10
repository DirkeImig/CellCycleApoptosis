%%%%%%%%% Fig.2, C,D,E: correlations %%%%%%%%%%%%%%%%%%%%%%%

cd ..
run Import_Cycle_Data.m
cd Fig2

%%%%%%% 0 for no exclusion of data points, 1 for exclusion of outliers %%%%
Excluded=1; 
G1_plot = G1(cens_Cycle==0);
CC_plot = CC(cens_Cycle==0);
SG2M_plot = SG2M(cens_Cycle==0);



%%%%%%%%%%%%%%%%%%%%%%% Correlation of Division and G1 %%%%%%%%%%%%%%%%%%%%
c1=[0:30];
figure
axis([2 25 5 45]);
if Excluded==0
    fit2 = fit(G1_plot,CC_plot,'poly1');   
else
    fit2 = fit(G1_plot,CC_plot,'poly1','Exclude',G1_plot>14 | SG2M_plot>14);
end
ci = confint(fit2,0.95);
conf1=ci(1,1)*c1+ci(1,2);
conf2=ci(2,1)*c1+ci(2,2);
hold on;
scatter(G1_plot,CC_plot,80,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0.4,0.4,0.4]);
p=plot(fit2,'--b');
set(p, 'LineWidth',3);
if Excluded==1
    hold on;
    z=0;
    for i=1:length(G1_plot)
        if G1_plot(i)>14 | SG2M_plot(i)>14
            z=z+1;
            scatter(G1_plot(i),CC_plot(i),80,'MarkerEdgeColor','k','MarkerFaceColor','r');
        end
    end
end
hold on;
p=plot(fit2,'--b');
hold on;
j=plot(c1,conf1,'--');
hold on;
k=plot(c1,conf2,'--k');
set(p, 'LineWidth',3);
set(j, 'LineWidth',2,'Color',[0.8,0.8,0.8]);
set(k, 'LineWidth',2,'Color',[0.8,0.8,0.8]);
grid on;
box on;
xlabel('G_1 [h]');
ylim([0,50]);
ylabel('Cell cycle [h]');
if Excluded==0
    h=legend(['m=' num2str(round(fit2.p1*100)/100),', a=' num2str(round(fit2.p2*100)/100)],['n=' num2str(length(G1_plot))],'Location','northwest');   
else
    h=legend(['n=' num2str(length(G1_plot)-z)],['m=' num2str(round(fit2.p1*100)/100),', a=' num2str(round(fit2.p2*100)/100)],'Location','northwest');   
end
set(h,'FontSize',22);
set(gca,'FontSize',24);
ax = gca;
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);   
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);    

%%%%%%%%% Pearson correlation coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[R,P,RLO,RUP]=corrcoef(G1_plot(G1_plot<14 & SG2M_plot<14),CC_plot(G1_plot<14 & SG2M_plot<14),'alpha',0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%% Correlation of Division and S/G2/M %%%%%%%%%%%%%%%%%%%%%
c1=[0:25];
figure
axis([2 25 5 45]);
if Excluded==0
    fit3 = fit(SG2M_plot,CC_plot,'poly1');  
else
    fit3 = fit(SG2M_plot,CC_plot,'poly1','Exclude',SG2M_plot>14 | G1_plot>14);
end
ci = confint(fit3,0.95);
conf1=ci(1,1)*c1+ci(1,2);
conf2=ci(2,1)*c1+ci(2,2);
hold on;
scatter(SG2M_plot,CC_plot,80,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0.4,0]);
p=plot(fit3,'--b');
set(p, 'LineWidth',3);
hold on;
if Excluded==1
    hold on;
    z=0;
    for i=1:length(SG2M_plot)
        if SG2M_plot(i)>14 | G1_plot(i)>14
            z=z+1;
            scatter(SG2M_plot(i),CC_plot(i),80,'MarkerEdgeColor','k','MarkerFaceColor','r');
        end
    end
end
hold on;
p=plot(fit3,'--b');
j=plot(c1,conf1,'--');
hold on;
k=plot(c1,conf2,'--k');
set(p, 'LineWidth',3);
set(j, 'LineWidth',2,'Color',[0.8,0.8,0.8]);
set(k, 'LineWidth',2,'Color',[0.8,0.8,0.8]);
grid on;
box on;
ylabel('Cell Cycle [h]');
xlabel('S/G_2/M [h]');
ylim([0,50]);
if Excluded==0
    h=legend(['m=' num2str(round(fit3.p1*100)/100),', a=' num2str(round(fit3.p2*100)/100)],['n=' num2str(length(SG2M_plot))],'Location','northwest');   %,['censored=' num2str(length(SG2M(~isnan(SG2M) & cens_Cycle==1)))],'Location','northwest');
else
    h=legend(['n=' num2str(length(SG2M_plot)-z)],['m=' num2str(round(fit3.p1*100)/100),', a=' num2str(round(fit3.p2*100)/100)],'Location','northwest');   %,['censored=' num2str(length(SG2M(~isnan(SG2M) & cens_Cycle==1)))],['excluded=' num2str(z)],'Location','northwest');
end
set(h,'FontSize',22);
set(gca,'FontSize',24);
ax = gca;
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);  
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);    


%%%%%%%%% Pearson correlation coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[R2,P2,RLO2,RUP2]=corrcoef(SG2M_plot(SG2M_plot<14 & G1_plot<14),CC_plot(SG2M_plot<14 & G1_plot<14),'alpha',0.05);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%% Correlation of G1 and S/G2/M %%%%%%%%%%%%%%%%%%%%%%
c2=[0:25];
figure
axis([0 25 0 25]);
if Excluded==0
    fit1 = fit(SG2M_plot,G1_plot,'poly1');  
else
    fit1 = fit(SG2M_plot,G1_plot,'poly1','Exclude',G1_plot>14 | SG2M_plot>14);
end
ci = confint(fit1,0.95);
conf1=ci(1,1)*c2+ci(1,2);
conf2=ci(2,1)*c2+ci(2,2);
hold on;
scatter(SG2M_plot,G1_plot,80,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1]);
hold on;
p=plot(fit1,'--b');
set(p, 'LineWidth',3);
if Excluded~=0
    z=0;
    for i=1:length(SG2M_plot)
        if G1_plot(i)>14 | SG2M_plot(i)>14
            z=z+1;
            scatter(SG2M_plot(i),G1_plot(i),80,'MarkerEdgeColor','k','MarkerFaceColor','r');
        end
    end
end
hold on;
j=plot(c2,conf1,'--');
hold on;
k=plot(c2,conf2,'--k');
hold on;
p=plot(fit1,'--b');
set(p, 'LineWidth',3);
set(j, 'LineWidth',2,'Color',[0.8,0.8,0.8]);
set(k, 'LineWidth',2,'Color',[0.8,0.8,0.8]);
grid on;
box on;
xlabel('S/G_2/M [h]');
ylabel('G_1 [h]');
if Excluded==0
   h=legend(['m=' num2str(round(fit1.p1*100)/100),', a=' num2str(round(fit1.p2*100)/100)],['n=' num2str(length(SG2M_plot))],'Location','northwest');    %,['censored=' num2str(length(SG2M_cens))],'Location','northwest');
else
    h=legend(['n=' num2str(length(SG2M_plot)-z)],['m=' num2str(round(fit1.p1*100)/100),', a=' num2str(round(fit1.p2*100)/100)],'Location','northwest');  %,['censored=' num2str(length(SG2M_cens))],['excluded=' num2str(z)],'Location','northwest');
end
set(h,'FontSize',22);
set(gca,'FontSize',24);
ax = gca;
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);   
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);   

%%%%%%%%% Pearson correlation coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[R3,P3,RLO3,RUP3]=corrcoef(SG2M_plot(SG2M_plot<14 & G1_plot<14),G1_plot(SG2M_plot<14 & G1_plot<14),'alpha',0.05)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





