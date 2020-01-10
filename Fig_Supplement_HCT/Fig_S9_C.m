%%%%%% IntiQuan IQM Tools is needed: https://iqmtools.intiquan.com %%%%%%%%
%%%%%% exemplary Model Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=0.12; 
PAS=0.32;
Sample='Data_Model_Sample.mat';

model = IQMmodel('my_model.txt');
IQMmakeMEXmodel(model,'MEXmodel1');
MyModel='MEXmodel1';

[~,Cyc_gr,Cyc_wh,Td_gr,Td_wh,d_phase_gr,d_phase_wh]=Myfun(p,PAS,Sample,MyModel);

Cycle_Start_Model=[Cyc_wh, Cyc_gr];
Cycle_Start_Model=round(100.*Cycle_Start_Model)./100;
T_Death_Model=[Td_wh, Td_gr];
d_phase=[d_phase_wh, d_phase_gr];
Cycle_Start_Model_NEW=[Cycle_Start_Model,Cycle_Start_Model(d_phase>2)];  
T_Death_Model_NEW=[T_Death_Model,T_Death_Model(d_phase>2)];
d_phase_NEW=[d_phase,d_phase(d_phase>2)];

w=6.9;
g=10.7;
CWE=w/(w+g);


figure()
x=[CWE   0        0          1 ];
y=[0    w      (w+g)      0  ];
p1=patch(x,y,-1*ones(size(x)),[0 0.6 0],'Linestyle','none');
p1.FaceAlpha = 0.1 ;
p1.FaceVertexAlphaData=[];   


hold on;
x2=[1          0           0          1      ];
y2=[w       (w+g+w)      (w+g+w+g)         (w+g)];
p1=patch(x2,y2,-1*ones(size(x2)),[0 0.6 0],'Linestyle','none');
p1.FaceAlpha = 0.1 ;
p1.FaceVertexAlphaData=[];
hold on;

y=[(w+g),0];
x=[0,1];
plot(x,y,'Color',[0.6,0.8,0.6],'LineWidth',3,'LineStyle','--');
hold on;

y=[2*(w+g),(w+g)];
x=[0,1];
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
    elseif d_phase_NEW(i)==5
        scatter(Cycle_Start_Model_NEW(i),T_Death_Model_NEW(i),45,[0,0,0],'d','MarkerFaceColor',[1,1,1],'MarkerEdgeColor',[0,0,0]);
        hold on;
    end
end

    
y=[0,30];
x=[CWE,CWE]
plot(x,y,'k--');
xlabel(' ');
ylabel('t_{death}[h]');
ax = gca;
set(gca,'xtick',[]);
box on;
set(gca,'FontSize',20);
ax.YLim = [0 24.75];
set(gca,'ytick',[0,4,8,12,16,20,24]);
set(gca,'xticklabels',({'0,4,8,12,16,20,24'})); 

