%%%%%%%%%%%%%%%%% exemplary plot parameter estimation %%%%%%%%%%%%%%%%%%%%%

n=12;
p=linspace(0.1,0.18,n);   
PAS=linspace(0.41,0.58,n);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('ModelFit_0.mat')

C(C==min(C))
Opti_p_1=Plot_p(C==min(C))
Opti_PAS_1=Plot_PAS(C==min(C))


Opti=min(C)
Area_PAS=[];
Area_p=[];

for i=1:length(C)
    if C(i)<Opti+(5.99/2)
        Area_PAS=[Area_PAS,Plot_PAS(i)];
        Area_p=[Area_p,Plot_p(i)];
    end
end

Min_PAS_1=min(Area_PAS)
Max_PAS_1=max(Area_PAS)
Min_p_1=min(Area_p)
Max_p_1=max(Area_p)
Opti_1=min(C)

Z=[C(1:12);...
    C(13:24);...
    C(25:36);...
    C(37:48);...
    C(49:60);...
    C(61:72);...
    C(73:84);...
    C(85:96);...
    C(97:108);...
    C(109:120);...
    C(121:132);...
    C(133:144)
    ];

v=[Opti+0.001,Opti+(5.99/2),Opti+(9.21/2),Opti+(13.82/2)];

figure()
contourf(PAS,p,Z,v)
ylabel('p')
xlabel('PAS')
hold on;
for i=1:length(Plot_p)
    if C(i)<Opti+0.001
        scatter(Plot_PAS(i),Plot_p(i),150,'+k','MarkerFaceColor',[0,0,0],'LineWidth',2.5)
        hold on;
    end
end
hold on;
map=[1,1,1
    0.6,0.8,1
    0.4,0.6,1
    0.1,0.3,1
    ]
colormap(map)
h=colorbar('Ticks',[1.841e+03,1.8428e+03,1.8445e+03,1.8462e+03],'TickLabels',{'0.05','0.01','0.001','<0.001'})
set(h, 'FontSize',20, 'location','northoutside');
ax = gca;
set(gca,'FontSize',20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


