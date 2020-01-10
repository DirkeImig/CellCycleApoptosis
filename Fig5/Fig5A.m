%%%%%%%%%%%%%%%%% Dotplot experimental data (t_death over C_0) %%%%%%%%%%%%
clear all

cd ..
run Mean_Cycle_Start.m
cd Fig5

le=0;

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

phase_Data_wh=[];  %1: death in G1 f0, 2: death in SG2M f0, 3: death in G1 f1

for j = 1: length(t_TRAIL_Div_F0)
    if t_TRAIL_WhiteEnd_F0(j)>0 && isnan(t_TRAIL_Div_F0(j))                                   %cells died in SG2M, f0
        p3=scatter(Mean_CycleStart_TR_in_white(j),t_TRAIL_Death(j),40,[0,0.4,0.2],'^','MarkerFaceColor',[0,0.4,0.2],'MarkerEdgeColor',[0,0,0]);  
        phase_Data_wh=[phase_Data_wh,2];
        hold on;
    elseif isnan(t_TRAIL_WhiteEnd_F0(j)) && (t_TRAIL_Death(j) + t_Birth_TRAIL_white(j))<25   %cells died in G1, f0 (and OVERALL G1 phase smaller than X h)
        phase_Data_wh=[phase_Data_wh,1];
        p4=scatter(Mean_CycleStart_TR_in_white(j),t_TRAIL_Death(j),40,[0,0,0],'^','MarkerFaceColor',[1,1,1],'MarkerEdgeColor',[0,0,0]);
        hold on;
    elseif isnan(t_TRAIL_WhiteEnd_F1(j)) && t_TRAIL_WhiteEnd_F0(j)>0 && t_TRAIL_Div_F0(j)>0   %cells died in G1, f1
        p5=scatter(Mean_CycleStart_TR_in_white(j),t_TRAIL_Death(j),45,[0,0,0],'square','MarkerFaceColor',[1,1,1],'MarkerEdgeColor',[0,0,0]);
        phase_Data_wh=[phase_Data_wh,3];
        hold on;
    elseif t_TRAIL_WhiteEnd_F1(j)>0 && t_TRAIL_WhiteEnd_F0(j)>0 && t_TRAIL_Div_F0(j)>0        %cells died in SG2M ,f1
        o=scatter(Mean_CycleStart_TR_in_white(j),t_TRAIL_Death(j),45,[0,0.4,0.2],'square','MarkerFaceColor',[0,0.4,0.2],'MarkerEdgeColor',[0,0,0]) ;
        le=le+1;
        phase_Data_wh=[phase_Data_wh,4];
        hold on;
    end
end


       
hold on;
phase_Data_gr=[];  %1: death in G1, f0, 2: SG2M, f0, 3: G1 f1

for i = 1: length(t_TRAIL_Div_F0_TRinGreen)
    if t_TRAIL_Div_F0_TRinGreen(i)>0 && isnan(t_TRAIL_WhiteEnd_F1_TRinGreen(i))           %cells died in G1, f1
        p1=scatter(Mean_CycleStart_TR_in_green(i),t_TRAIL_Death_TRinGreen(i),45,[0,0,0],'square','MarkerFaceColor',[1,1,1],'MarkerEdgeColor',[0,0,0]);
        phase_Data_gr=[phase_Data_gr,3];
        hold on;
    elseif isnan(t_TRAIL_Div_F0_TRinGreen(i))                                             %cells died in SG2M, f0 
        p2=scatter(Mean_CycleStart_TR_in_green(i),t_TRAIL_Death_TRinGreen(i),40,[0,0.4,0.2],'^','MarkerFaceColor',[0,0.4,0.2],'MarkerEdgeColor',[0,0,0]);
        phase_Data_gr=[phase_Data_gr,2];
        hold on;
    elseif t_TRAIL_Div_F0_TRinGreen(i)>0 & t_TRAIL_WhiteEnd_F1_TRinGreen(i)>0             %SG2M, f1
        o=scatter(Mean_CycleStart_TR_in_green(i),t_TRAIL_Death_TRinGreen(i),45,[0,0.4,0.2],'square','MarkerFaceColor',[0,0.4,0.2],'MarkerEdgeColor',[0,0,0]);
        phase_Data_gr=[phase_Data_gr,4];
        hold on;
    else
        o=scatter(Mean_CycleStart_TR_in_green(i),t_TRAIL_Death_TRinGreen(i),45,[0.5,0.5,0.5],'square') ;
        phase_Data_gr=[phase_Data_gr,5];
        hold on;
    end
end


hold on;
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


[p,h]=ranksum(t_TRAIL_Death,t_TRAIL_Death_TRinGreen)




%%%%%%%%%%%%%%%%%%% Plot Outliers G1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=t_TRAIL_Death(isnan(t_TRAIL_WhiteEnd_F0) & (t_TRAIL_Death + t_Birth_TRAIL_white)>25);
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
            scatter(1+f,R(j),'^','MarkerFaceColor',[1,1,1],'MarkerEdgeColor',[0,0,0]);
            hold on;
            f=f+1;
        end
    end
end

plot([-9,9],[mean(R),mean(R)],'k','LineWidth',4);

ylabel(' ');
set(gca,'yticklabels',({''}));
ax = gca;
set(gca,'xtick',[1]);
set(gca,'xticklabels',({''}));      %Non-dividing
box on;
ax.XLim = [-9.5 9.5];
ax.YLim = [0 19.5];
set(gca,'FontSize',20);

