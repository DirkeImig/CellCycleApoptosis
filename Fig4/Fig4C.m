%%%%%%%%%%%%% C_0 distribution in apoptotic, surviving and all cells %%%%%%

clear all

run Fig4C_probabilities

%%%%%%%%%%%%%%%%%%% I) All G1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=All_white;
R=round(100.*Y)./100;
C=unique(R);

count=zeros(1,length(C));
for i=1:length(C)
    for j=1:length(R) 
        if R(j)==C(i)
            count(i)=count(i)+1;
        end
    end
end
    
figure()       
for i=1:length(C)
    f=-count(i)/2;
    for j=1:length(R)
        if R(j)==C(i)
            scatter(1+f,R(j),'k')  
            hold on;
            f=f+1;
        end
    end
end

plot([-9,9],[median(R),median(R)],'Color',[0.5,0.5,0.5],'LineWidth',4);

ylabel('$\hat{C}_0$','Interpreter','latex');
ax = gca;
set(gca,'xtick',[1]);
set(gca,'xticklabels',({'All'}));
box on;
ax.XLim = [-9.5 9.5];
ax.YLim = [0 0.44];
str = {'n=565'};
text(-8,0.4,str,'FontSize',20,'EdgeColor','k','BackgroundColor','w')
set(gca,'FontSize',20);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 5]);   
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 3 5]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%% II) Survivors G1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=Mean_CycleStart_TR_in_white_SURVIVORS;
R=round(100.*Y)./100;
C=unique(R);

count=zeros(1,length(C));
for i=1:length(C)
    for j=1:length(R) 
        if R(j)==C(i)
            count(i)=count(i)+1;
        end
    end
end
    
figure()   
for i=1:length(C)
    f=-count(i)/2;
    for j=1:length(R)
        if R(j)==C(i)
            scatter(1+f,R(j),'k');
            hold on;
            f=f+1;
        end
    end
end

plot([-9,9],[median(R),median(R)],'Color',[0.5,0.5,0.5],'LineWidth',4);

ylabel('G_1 (0)');
ax = gca;
set(gca,'xtick',[1]);
set(gca,'xticklabels',({'Survivor'}));
box on;
ax.XLim = [-9.5 9.5];
ax.YLim = [0 0.44];
str = {'n=143'};
text(-8,0.4,str,'FontSize',20,'EdgeColor','k','BackgroundColor','w')
set(gca,'FontSize',20);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 5]);   
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 3 5]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%% III) Apoptotic G1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=Mean_CycleStart_TR_in_white_DEAD;
R=round(100.*Y)./100;
C=unique(R);

count=zeros(1,length(C));
for i=1:length(C)
    for j=1:length(R) 
        if R(j)==C(i)
            count(i)=count(i)+1;
        end
    end
end
    
figure()   
for i=1:length(C)
    f=-count(i)/2;
    for j=1:length(R)
        if R(j)==C(i)
            scatter(1+f,R(j),'k');
            hold on;
            f=f+1;
        end
    end
end

plot([-9,9],[median(R),median(R)],'Color',[0.5,0.5,0.5],'LineWidth',4)

ylabel('G_1 (0)');
ax = gca;
set(gca,'xtick',[1]);
set(gca,'xticklabels',({'Apoptotic'}));
box on;
ax.XLim = [-9.5 9.5];
ax.YLim = [0 0.44];
str = {'n=417'};
text(-8,0.4,str,'FontSize',20,'EdgeColor','k','BackgroundColor','w')
set(gca,'FontSize',20);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 5]);   
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 3 5]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%% IV) All SG2M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=All_green;
R=round(100.*Y)./100;
C=unique(R);

count=zeros(1,length(C));
for i=1:length(C)
    for j=1:length(R) 
        if R(j)==C(i)
            count(i)=count(i)+1;
        end
    end
end
    
figure()   
for i=1:length(C)
    f=-count(i)/2;
    for j=1:length(R)
        if R(j)==C(i)
            scatter(1+f,R(j),'k');
            hold on;
            f=f+1;
        end
    end
end

plot([-9,9],[median(R),median(R)],'Color',[0.5,0.5,0.5],'LineWidth',4);

ylabel('$\hat{C}_0$','Interpreter','latex');
ax = gca;
set(gca,'xtick',[1]);
set(gca,'xticklabels',({'All'}));
box on;
ax.XLim = [-9.5 9.5];
ax.YLim = [0.44 1];
str = {'n=490'};
text(-8,0.95,str,'FontSize',20,'EdgeColor','k','BackgroundColor','w')
set(gca,'FontSize',20);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 5]);  
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 3 5]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%% V) Survivor SG2M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=Mean_CycleStart_TR_in_green_SURVIVORS;
R=round(100.*Y)./100;
C=unique(R);

count=zeros(1,length(C));
for i=1:length(C)
    for j=1:length(R) 
        if R(j)==C(i)
            count(i)=count(i)+1;
        end
    end
end
    
figure()   
for i=1:length(C)
    f=-count(i)/2;
    for j=1:length(R)
        if R(j)==C(i)
            scatter(1+f,R(j),'k');
            hold on;
            f=f+1;
        end
    end
end

plot([-9,9],[median(R),median(R)],'Color',[0.5,0.5,0.5],'LineWidth',4)

ylabel('S/G_2/M (0)');
ax = gca;
set(gca,'xtick',[1]);
set(gca,'xticklabels',({'Survivor'}));
box on;
ax.XLim = [-9.5 9.5];
ax.YLim = [0.44  1];
str = {'n=140'};
text(-8,0.95,str,'FontSize',20,'EdgeColor','k','BackgroundColor','w')
set(gca,'FontSize',20);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 5]);   
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 3 5]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%% VI) Apoptotic SG2M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=Mean_CycleStart_TR_in_green_DEAD;
R=round(100.*Y)./100;
C=unique(R);

count=zeros(1,length(C));
for i=1:length(C)
    for j=1:length(R) 
        if R(j)==C(i)
            count(i)=count(i)+1;
        end
    end
end
    
figure()   
for i=1:length(C)
    f=-count(i)/2;
    for j=1:length(R)
        if R(j)==C(i)
            scatter(1+f,R(j),'k');
            hold on;
            f=f+1;
        end
    end
end

plot([-9,9],[median(R),median(R)],'Color',[0.5,0.5,0.5],'LineWidth',4);

ylabel('S/G_2/M (0)');
ax = gca;
set(gca,'xtick',[1]);
set(gca,'xticklabels',({'Apoptotic'}));
box on;
ax.XLim = [-9.5 9.5];
ax.YLim = [0.44 1];
str = {'n=321'};
text(-8,0.95,str,'FontSize',20,'EdgeColor','k','BackgroundColor','w')
set(gca,'FontSize',20);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 5]);   
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 3 5]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




