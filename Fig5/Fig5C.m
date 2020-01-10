%%%%%%%%%%%%%%%%% exemplary dotplot model data (t_death over C_0) %%%%%%%%%

p=0.15;  %best 
PAS=0.52;
Sample='Data_Model_Sample_LESS.mat';

[~,Cyc_gr,Cyc_wh,Td_gr,Td_wh,d_phase_gr,d_phase_wh]=myfun(p,PAS,Sample);

Cycle_Start_Model=[Cyc_wh, Cyc_gr];
Cycle_Start_Model=round(100.*Cycle_Start_Model)./100;
T_Death_Model=[Td_wh, Td_gr];
d_phase=[d_phase_wh, d_phase_gr];
Cycle_Start_Model_NEW=[Cycle_Start_Model,Cycle_Start_Model(d_phase>2)];  %dividing cells
T_Death_Model_NEW=[T_Death_Model,T_Death_Model(d_phase>2)];
d_phase_NEW=[d_phase,d_phase(d_phase>2)];


figure()
x=[0.45   0        0          1 ];
y=[0    7.41      16.43      0  ];
p1=patch(x,y,-1*ones(size(x)),[0 0.6 0],'Linestyle','none');
p1.FaceAlpha = 0.1 ;
p1.FaceVertexAlphaData=[];   
hold on;
x2=[1          0           0          1      ];
y2=[7.41    23.84      32.86         16.43   ];
p1=patch(x2,y2,-1*ones(size(x2)),[0 0.6 0],'Linestyle','none');
p1.FaceAlpha = 0.1 ;
p1.FaceVertexAlphaData=[];
hold on;

y=[16.43,0];
x=[0,1]
plot(x,y,'Color',[0.6,0.8,0.6],'LineWidth',3,'LineStyle','--');
hold on;

y=[32.86,16.43];
x=[0,1]
p4=plot(x,y,'Color',[0.6,0.8,0.6],'LineWidth',3,'LineStyle','--');
hold on;


for i=1:length(Cycle_Start_Model_NEW)
    if d_phase_NEW(i)==1
        scatter(Cycle_Start_Model_NEW(i),T_Death_Model_NEW(i),40,'^','MarkerFaceColor',[1,1,1],'MarkerEdgeColor',[0,0,0]);   %death in G1 f0
        hold on;
    elseif d_phase_NEW(i)==2      
        scatter(Cycle_Start_Model_NEW(i),T_Death_Model_NEW(i),40,'^','MarkerFaceColor',[0,0.4,0.2],'MarkerEdgeColor',[0,0,0]);   %death in SG2M f0
        hold on;
    elseif d_phase_NEW(i)==3
        scatter(Cycle_Start_Model_NEW(i),T_Death_Model_NEW(i),40,'square','MarkerFaceColor',[1,1,1],'MarkerEdgeColor',[0,0,0]); %death in G1 f1
        hold on;
    elseif d_phase_NEW(i)==4
        scatter(Cycle_Start_Model_NEW(i),T_Death_Model_NEW(i),40,'square','MarkerFaceColor',[0,0.4,0.2],'MarkerEdgeColor',[0,0,0]); %death in SG2M f1
        hold on;
    end
end

y=[0,30];
x=[0.45,0.45]
plot(x,y,'k--');
xlabel(' ');
ylabel('t_{death}[h]');
ax = gca;
set(gca,'xtick',[]);
box on;
set(gca,'FontSize',20);
ax.YLim = [0 19.5];
set(gca,'ytick',[0,2,4,6,8,10,12,14,16,18]);
set(gca,'xticklabels',({'0,2,4,6,8,10,12,14,16,18'})); 












