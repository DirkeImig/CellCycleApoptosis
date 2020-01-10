%%%%%%%%%% Fig.2, A,B: histograms %%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..
run Import_Cycle_Data.m
cd Fig2


%%%%%%%%%%%%%%%%%%%%%%%%%%%% G1 PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% #cells G1 phase (without phase > 38h) %%%%%%%%%%%%%%%%%%%%%

figure()
subplot('Position',[0.11 0.1790 0.65 0.7460])   
[fa,xa]=hist(G1(cens_G1==0), 70);   % without censored data
h=bar(xa,fa)
set(h,'Facecolor',[106,104,117]./ 255,'EdgeColor','k');
hold on;
xt = 0:0.1:30;
plot(xt,trapz(xa,fa)*lognpdf(xt,log_G1(1),log_G1(2)),'Color',[0.2,0.2,0.2],'LineWidth',2)
xlabel('G_1 [h]', 'FontSize',14)
ylabel('# cells', 'FontSize',14)
set(gca,'FontSize',16,'ylim',([0,max(fa)]),'xlim',([0,max(G1(cens_G1==0))+1]))
legend(['n= ' num2str(length(G1(cens_G1==0)))],'Location','northeast');

subplot('Position',[0.8 0.1790 0.1 0.7460])    
[fa1,xa1]=hist(G1(cens_G1==1), 1);   %only censored data (in red)
b=bar(xa1,fa1)
xticks([41])
xticklabels({'>38'})
yticks([ ])
set(b,'Facecolor',[1,0,0],'EdgeColor','k');
set(gca,'FontSize',16,'ylim',([0,max(fa)]))  



%%%%%%%%%%%%%%%%%%%%%%%%%%%% SG2M PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% #cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
%subplot('Position',[0.11 0.1790 0.65 0.7460])   %%% only for HCT %%%%%%%%%
[fa,xa]=hist(SG2M(cens_SG2M==0), 70);   % without censored data
h=bar(xa,fa)
hold on;
[fa1,xa1]=hist(SG2M(cens_SG2M==1 & SG2M<=30), 70);   %only censored data (in red)
b=bar(xa1,fa1)
set(h,'Facecolor',[47,167,10]./ 255,'EdgeColor','k');
set(b,'Facecolor',[1,0,0],'EdgeColor','k');
hold on;
xt = 0:0.1:25;
plot(xt,trapz(xa,fa)*lognpdf(xt,log_SG2M(1),log_SG2M(2)),'Color',[0,0.4,0],'LineWidth',2)
xlabel('S/G_2/M [h]', 'FontSize',14)
ylabel('# cells', 'FontSize',14)
set(gca,'FontSize',16,'ylim',([0,max(fa)]),'xlim',([0,max(SG2M)+1]))
legend(['n= ' num2str(length(SG2M(cens_SG2M==0)))],'Location','northeast');

%%%%%%%% only for HCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
subplot('Position',[0.8 0.1790 0.1 0.7460])   
[fa1,xa1]=hist(SG2M(cens_SG2M==1 & SG2M>30), 1);     
b=bar(xa1,fa1)
xticks([41])
xticklabels({'>30'})
yticks([ ])
set(b,'Facecolor',[1,0,0],'EdgeColor','k');
set(gca,'FontSize',16,'ylim',([0,max(fa)]))  
%}



