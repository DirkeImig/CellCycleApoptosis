%%%%%%%%%%%  average apoptotic progression over age %%%%%%%%%%%%%%%%%%%%%%%%

Zyklus=16.3
g=1/Zyklus    %growth rate 
tPAD=0.52/g   %time of PAD
m=0.4         %slope before PAD
mp=0.25       %slope after PAD: m minus p


%%%%%% TRAIL before PAD %%%%%%%%%
t0=linspace(0,tPAD,4)      
a_tPAD=m*tPAD - t0.*m
a0=0

a_t_Div=a_tPAD + mp*Zyklus - tPAD*mp     
t_Div=16.3

tend=21
a_tend=a_t_Div + m*tend - t_Div*m



t0_d=linspace(10.8667,Zyklus,3)
a_t_Div_d=mp*Zyklus - t0_d.*mp
a0_d=0
a_tend_d=a_t_Div_d + m*tend - t_Div*m



figure()
for i=1:length(t0)
   k=plot([t0(i),tPAD,   t_Div,  tend], [a0, a_tPAD(i), a_t_Div(i), a_tend(i)],'LineWidth',4,'Color',[0,0,0])
   hold on;
end
hold on;
for i=1:length(t0_d)
   l=plot([t0_d(i), t_Div, tend], [a0_d, a_t_Div_d(i), a_tend_d(i)],'LineWidth',4,'Color',[0,0,0])
   hold on;
end
x = [tPAD tPAD Zyklus Zyklus tPAD];
y = [0 1 1 0 0];
patch(x, y, -1 * ones(size(x)), [0.7 0.7 0.7], 'LineStyle', 'none')
hold on;
plot([tPAD,tPAD],[0,1],'Color',[0.6,0.6,0.6],'LineWidth',1)
hold on;
plot([Zyklus,Zyklus],[0,1],'Color',[0.6,0.6,0.6],'LineWidth',1)
ylim([0,1])
xlim([0,19])
yticks([0,1])
ylabel('apoptosis')
xlabel('age [h]')
set(gca,'FontSize',22)


