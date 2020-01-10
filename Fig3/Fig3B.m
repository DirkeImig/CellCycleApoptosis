%%%%%%%%% histogram of S/G2/M phase lengths in treated cells %%%%%%%%%%%%%%%%%%

clear all

%%%%%%%%%%%%%%%%%%% data  unstimulated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
run Import_Cycle_Data.m;            %cell cycle data of untreated cells
cens=0                              %0 without censored data

Para_green_MODEL_Apop=[2.1463, 0.1909];           %only apoptotic
Para_green_MODEL_Surv_AND_Apop=[2.1602, 0.1918];  %survivors as in data



%%%%%%%%%%%%%%%%%%%%%% data  stimulated in S/G2/M %%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'Dead_Survivor_green.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
generation = dataArray{:, 1};
cell = dataArray{:, 2};
birthwhitestart= dataArray{:, 3};
whiteendgreenstart = dataArray{:, 4};
greenend = dataArray{:, 5};
roundstartMstart = dataArray{:, 6};
celldivision = dataArray{:, 7};
celldeath = dataArray{:, 8};
celllost = dataArray{:, 9};
survived = dataArray{:, 10};
time = dataArray{:, 11};
videoend = dataArray{:, 12};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
cd Fig3

%%%%%% 1: cell division: end of phase is known
Uncens_Length_SG2M_F0 = (celldivision(~isnan(celldivision) & generation==0 & ~isnan(whiteendgreenstart)) - whiteendgreenstart(~isnan(celldivision) & generation==0 & ~isnan(whiteendgreenstart)) )*0.25;

%%%%%% 2: cell death as (censored) end of phase
A_Cens_Length_SG2M_F0 = (celldeath(~isnan(celldeath) & generation==0 & ~isnan(whiteendgreenstart)) - whiteendgreenstart(~isnan(celldeath) & generation==0 & ~isnan(whiteendgreenstart)) )*0.25;

%%%%% 3: no death, no division, cell survived --> end of video as censoring
B_Cens_Length_SG2M_F0 = (videoend(isnan(celldivision) & survived==1 & generation==0 & ~isnan(whiteendgreenstart)) - whiteendgreenstart(isnan(celldivision) & survived==1 & generation==0 & ~isnan(whiteendgreenstart)) )*0.25;
Data=[Uncens_Length_SG2M_F0', A_Cens_Length_SG2M_F0', B_Cens_Length_SG2M_F0'];
uncens=zeros(1,length(Uncens_Length_SG2M_F0));
censor=ones(1,length([A_Cens_Length_SG2M_F0', B_Cens_Length_SG2M_F0']));
vect=[uncens, censor];
alpha=0.05; %confidence interval

if cens==0
    New_Para_green=lognfit(Uncens_Length_SG2M_F0);
else
    New_Para_green=lognfit(Data, alpha, vect);
end


mean_Data=exp(New_Para_green(1) + (New_Para_green(2)^2)/2)
mean_Model_Surv_AND_Apop=exp(Para_green_MODEL_Surv_AND_Apop(1) + (Para_green_MODEL_Surv_AND_Apop(2)^2)/2 )


figure()
[fa,xa]=hist(Uncens_Length_SG2M_F0, 70);   %without censored data
h=bar(xa,fa)
set(h,'Facecolor',[47,167,10]./ 255,'EdgeColor','k');
hold on;
xt = 0:0.1:32;
plot(xt,trapz(xa,fa)*lognpdf(xt,New_Para_green(1),New_Para_green(2)),'Color',[0,0.4,0],'LineWidth',3);
hold on;
plot(xt,trapz(xa,fa)*lognpdf(xt,Para_green_MODEL_Surv_AND_Apop(1),Para_green_MODEL_Surv_AND_Apop(2)),'Color',[0.2,0.2,0.2],'LineWidth',3,'LineStyle','--');
set(gca,'FontSize',20);
leg{1}=['uncensored data, n=' num2str(length(Uncens_Length_SG2M_F0))];
leg{2}=['fitted to data (mean= ' num2str(round(100*mean_Data)/100) 'h)'];  
leg{3}=['expected (mean= ' num2str(round(100*mean_Model_Surv_AND_Apop)/100) 'h)']; 
h=legend(leg, 'Location','northeast');
set(h,'FontSize',14);
xlabel('S/G_2/M [h]');
ylabel('# cells');
xlim([0,max(Uncens_Length_SG2M_F0)+1]);
grid on;


[prob1,hyp1]=ranksum(Uncens_Length_SG2M_F0',SG2M(cens_SG2M==0)')    %nullhypothesis







