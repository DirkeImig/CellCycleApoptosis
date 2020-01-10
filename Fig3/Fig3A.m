%%%%%%%%% histogram of G1 phase lengths in treated cells %%%%%%%%%%%%%%%%%%

clear all

%%%%%%%%%%%%%%%%%%% data  unstimulated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
run Import_Cycle_Data.m;                          %cell cycle data of untreated cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cens=0                                            %0 without taking censored data into account

max_Exp_time=19.25;
paracycle1(1)=1.1062;                             %t_death G1
paracycle1(2)=0.5649;
C_WE=0.45;                                        %0.39 for HCT                                   
Para_white_MODEL_Apop=[1.8938, 0.2687];           %lognormal distr. phase length of apoptotic, uncensored
Para_white_MODEL_Surv_AND_Apop=[1.9259, 0.2725];  %lognormal distr. phase length including survivors (data)



%%%%%%%%%%%%%%%%%%% data  stimulated G1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'Dead_Survivor_white.csv';
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

survived((celldeath-time)*0.25>max_Exp_time)=1;
celldeath((celldeath-time)*0.25>max_Exp_time)=nan;

clearvars filename delimiter startRow formatSpec fileID dataArray ans;
cd Fig3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% 1.: end of G1 phase is known
Uncens_Length_G1_F0 = (whiteendgreenstart(~isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart)) - birthwhitestart(~isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart)) )*0.25;

%%%%%% 2.: cell death as (censored) end of phase
A_Cens_Length_G1_F0 = (celldeath(~isnan(celldeath) & isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart)) - birthwhitestart(~isnan(celldeath) & isnan(whiteendgreenstart) & generation==0 & ~isnan(birthwhitestart)) )*0.25;

%%%%% 3.: no death, no end of G1, cell survived --> end of video as censoring
B_Cens_Length_G1_F0 = (videoend(isnan(whiteendgreenstart) & survived==1 & generation==0 & ~isnan(birthwhitestart)) - birthwhitestart(isnan(whiteendgreenstart) & survived==1 & generation==0 & ~isnan(birthwhitestart)) )*0.25;
Data=[Uncens_Length_G1_F0', A_Cens_Length_G1_F0', B_Cens_Length_G1_F0'];
uncens=zeros(1,length(Uncens_Length_G1_F0));
censor=ones(1,length([A_Cens_Length_G1_F0', B_Cens_Length_G1_F0']));
vect=[uncens, censor];
alpha=0.05; %confidence interval

if cens==0
    New_Para_G1=lognfit(Uncens_Length_G1_F0);
else
    New_Para_G1=lognfit(Data, alpha, vect);
end

mean_Data=exp(New_Para_G1(1) + (New_Para_G1(2)^2)/2 );                   %Data
mean_Model_Surv_AND_Apop=exp(Para_white_MODEL_Surv_AND_Apop(1) + (Para_white_MODEL_Surv_AND_Apop(2)^2)/2 );

figure()
[fa,xa]=hist(Uncens_Length_G1_F0, 70);                                   %without censored data
h=bar(xa,fa);
set(h,'Facecolor',[0.4,0.4,0.4],'EdgeColor','k');
hold on;
xt = 0:0.1:50;
plot(xt,trapz(xa,fa)*lognpdf(xt,New_Para_G1(1),New_Para_G1(2)),'Color',[0.6,0.6,0.6],'LineWidth',3);
hold on;
plot(xt,trapz(xa,fa)*lognpdf(xt,Para_white_MODEL_Surv_AND_Apop(1),Para_white_MODEL_Surv_AND_Apop(2)),'Color',[0,0,0],'LineWidth',3,'LineStyle','--');
set(gca,'FontSize',26);
leg{1}=['uncensored data, n=' num2str(length(Uncens_Length_G1_F0))];
leg{2}=['fitted to data (mean= ' num2str(round(100*mean_Data)/100) 'h)']; 
leg{3}=['expected, (mean= ' num2str(round(100*mean_Model_Surv_AND_Apop)/100) 'h)'];
h=legend(leg, 'Location','northeast');
set(h,'FontSize',18);
xlabel('G_1');
ylabel('# cells');
xlim([0,25]);
ylim([0,15]);
grid on;
   




