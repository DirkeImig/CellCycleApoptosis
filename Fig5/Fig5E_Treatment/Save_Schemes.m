clear all;
p=0.15;    %best 
PAS=0.52;  %best

sync=45      %synchronization in cell cycle

if sync==10
    Sample=('Data_Model_Sample_10Per.mat')
elseif sync==45
    Sample=('Data_Model_Sample_S.mat')
elseif sync==20
    Sample=('Data_Model_Sample_20Per.mat')
elseif sync==30
    Sample=('Data_Model_Sample_30Per.mat')
elseif sync==40
    Sample=('Data_Model_Sample_40Per.mat')
elseif sync==50
    Sample=('Data_Model_Sample_50Per.mat')
elseif sync==60
    Sample=('Data_Model_Sample_60Per.mat')
elseif sync==70
    Sample=('Data_Model_Sample_70Per.mat')
elseif sync==80
    Sample=('Data_Model_Sample_80Per.mat')
elseif sync==90
    Sample=('Data_Model_Sample_90Per.mat')
elseif sync==100
    Sample=('Data_Model_Sample_100Per.mat')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
[~,Cyc_gr,Cyc_wh,T_Death,T_Death_wh,d_phase_gr,d_phase_wh]=myfun(p,PAS,Sample);
cd ..
cd Fig5E_Treatment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mitosis=length(T_Death(d_phase_gr>2)) + length(T_Death_wh(d_phase_wh>2));
T_Death_All=[T_Death,T_Death_wh];


%%%%%%%%% Figure: time of death %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
b1=bar(mean([T_Death,T_Death_wh]),'BarWidth', 0.7)
set(b1,'FaceColor',[0.7,0.7,0.7])
ylabel('t_{death}[h]')
xlabel(' ')
hold on;
errorbar(mean([T_Death,T_Death_wh]),std([T_Death,T_Death_wh]),'r.','CapSize',6);
set(gca,'XTickLabel',{''})    
box off;
set(gca,'FontSize',20,'ylim',([0 12]))
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [2 2.5]);   
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 1.8 2.5]);   
if sync==10
    print('Scheme_10Per','-dpdf')
elseif sync==45
    print('Scheme_S','-dpdf')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% Number of cells that undergo mitosis %%%%%%%%%%%%%%%%%%%%%%%%%%%%
DeadMitosis=mitosis/length([T_Death,T_Death_wh])  %number of cells that die after mitosis     
if sync==10
    PerAllMitosis=100*(DeadMitosis*0.77) + 23     %including survivors
elseif sync==45
    PerAllMitosis=100*(DeadMitosis*0.73)  + 27 
elseif sync==20
    PerAllMitosis=100*(DeadMitosis*0.77) + 23
elseif sync==30
    PerAllMitosis=100*(DeadMitosis*0.77) + 23
elseif sync==40
    PerAllMitosis=100*(DeadMitosis*0.77) + 23
elseif sync==50
    PerAllMitosis=100*(DeadMitosis*0.77) + 23
elseif sync==60
    PerAllMitosis=100*(DeadMitosis*0.73) + 27
elseif sync==70
    PerAllMitosis=100*(DeadMitosis*0.73) + 27
elseif sync==80
    PerAllMitosis=100*(DeadMitosis*0.73)  + 27
elseif sync==90
    PerAllMitosis=100*(DeadMitosis*0.73)  + 27
elseif sync==100
    PerAllMitosis=100*(DeadMitosis*0.73)  + 27
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Figure mitosis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
b1=bar(PerAllMitosis,'BarWidth', 0.7)
set(b1,'FaceColor',[0.4,0.6,1])
ylabel('% divisions')
xlabel(' ')
set(gca,'XTickLabel',{''})   
box off;
set(gca,'FontSize',20,'ylim',([0 100]))
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [2 2.5]);   
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 1.8 2.5]);   
if sync==10
    print('Mitosis_10Per','-dpdf')
elseif sync==45
    print('Mitosis_S','-dpdf')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% save data files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sync==10
    save('Sync_10','T_Death_All','PerAllMitosis')
elseif sync==45
    save('Sync_45','T_Death_All','PerAllMitosis')
elseif sync==20
    save('Sync_20','T_Death_All','PerAllMitosis')
elseif sync==30
    save('Sync_30','T_Death_All','PerAllMitosis')
elseif sync==40
    save('Sync_40','T_Death_All','PerAllMitosis')
elseif sync==50
    save('Sync_50','T_Death_All','PerAllMitosis')
elseif sync==60
    save('Sync_60','T_Death_All','PerAllMitosis')
elseif sync==70
    save('Sync_70','T_Death_All','PerAllMitosis')
elseif sync==80
    save('Sync_80','T_Death_All','PerAllMitosis')
elseif sync==90
    save('Sync_90','T_Death_All','PerAllMitosis')
elseif sync==100
    save('Sync_100','T_Death_All','PerAllMitosis')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
