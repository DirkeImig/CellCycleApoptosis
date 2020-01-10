%%%%%%%%%%%%%%% box plots of death times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

cd ..
run Mean_Cycle_Start.m
cd Fig2

%%%%%%%%%% G1 %%%%%%%%%%%%%
figure()
boxplot(t_TRAIL_Death)
ax = gca;
set(gca,'xtick',[1]);
set(gca,'xticklabels',({''}));      
ylabel('t_{death} [h]');
set(gca,'ytick',[0 2.5 5 7.5 10 12.5 15 17.5]);
xlabel('G_1');
ax.YLim = [0 19.5];
set(gca,'FontSize',16);


%%%%%%%%%% SG2M %%%%%%%%%%%%%
figure()
boxplot(t_TRAIL_Death_TRinGreen)
ax = gca;
set(gca,'xtick',[1]);
set(gca,'xticklabels',({''}));    
set(gca,'ytick',[0 2.5 5 7.5 10 12.5 15 17.5]);
ylabel('t_{death} [h]');
xlabel('S/G_2/M');
ax.YLim = [0 19.5];
set(gca,'FontSize',16);


[h,p]=ranksum(t_TRAIL_Death,t_TRAIL_Death_TRinGreen)







