%%%%%%%%%%%%%%%%%% synchronization strategies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Sync_10.mat')
T_Death_10=T_Death_All;
PerAllMitosis_10=PerAllMitosis;

load('Sync_20.mat')
T_Death_20=T_Death_All;
PerAllMitosis_20=PerAllMitosis;

load('Sync_30.mat')
T_Death_30=T_Death_All;
PerAllMitosis_30=PerAllMitosis;

load('Sync_40.mat')
T_Death_40=T_Death_All;
PerAllMitosis_40=PerAllMitosis;

load('Sync_50.mat')
T_Death_50=T_Death_All;
PerAllMitosis_50=PerAllMitosis;

load('Sync_60.mat')
T_Death_60=T_Death_All;
PerAllMitosis_60=PerAllMitosis;

load('Sync_70.mat')
T_Death_70=T_Death_All;
PerAllMitosis_70=PerAllMitosis;

load('Sync_80.mat')
T_Death_80=T_Death_All;
PerAllMitosis_80=PerAllMitosis;

load('Sync_90.mat')
T_Death_90=T_Death_All;
PerAllMitosis_90=PerAllMitosis;

load('Sync_100.mat')
T_Death_100=T_Death_All;
PerAllMitosis_100=PerAllMitosis;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

a=[mean(T_Death_10) mean(T_Death_20) mean(T_Death_30) mean(T_Death_40) mean(T_Death_50) mean(T_Death_60) mean(T_Death_70) mean(T_Death_80) mean(T_Death_90) mean(T_Death_100)];
b=[mean(PerAllMitosis_10) mean(PerAllMitosis_20) mean(PerAllMitosis_30) mean(PerAllMitosis_40) mean(PerAllMitosis_50) mean(PerAllMitosis_60) mean(PerAllMitosis_70) mean(PerAllMitosis_80) mean(PerAllMitosis_90) mean(PerAllMitosis_100)];

error_a=[mean(T_Death_10) mean(T_Death_20) mean(T_Death_30) mean(T_Death_40) mean(T_Death_50) mean(T_Death_60) mean(T_Death_70) mean(T_Death_80) mean(T_Death_90) mean(T_Death_100)];

fig=figure
left_color = [0,0,0];
right_color = [0.4,0.6,1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
x_a=[1:length(a)]
a1=bar(x_a,a)
set(a1,'FaceColor',[0.7,0.7,0.7],'BarWidth',0.3)
ylabel('t_{death}[h]')
xlabel(' ')
hold on;
errorbar([mean(T_Death_10) mean(T_Death_20) mean(T_Death_30) mean(T_Death_40) mean(T_Death_50) mean(T_Death_60) mean(T_Death_70) mean(T_Death_80) mean(T_Death_90) mean(T_Death_100)],[std(T_Death_10) std(T_Death_20) std(T_Death_30) std(T_Death_40) std(T_Death_50) std(T_Death_60) std(T_Death_70) std(T_Death_80) std(T_Death_90) std(T_Death_100)],'r.','CapSize',6);
box off;
yyaxis right
x_b=[1:length(b)]+0.4
b1=bar(x_b,b)
set(b1,'FaceColor',[0.4,0.6,1],'BarWidth',0.3)
ylabel('% divisions')
set(gca,'XTick',([1 2 3 4 5 6 7 8 9 10]+0.2),'xlim',([0.5 11]))
set(gca,'XTickLabel',{'0-0.1','0.1-0.2','0.2-0.3','0.3-0.4','0.4-0.5','0.5-0.6','0.6-0.7','0.7-0.8','0.8-0.9','0.9-1'})    %,{'f0', 'f1'})
xlabel('synchronization in initial cell cycle position C_0')









