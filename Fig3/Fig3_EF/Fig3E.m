%%%%%%%%%%% parameter estimation - exemplary %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% EXAMPLE: Calculation of Likelihood %%%%%%%%%%%%%%%%%%%%%%%%%%
%{
load('Data_Model_G1_0.mat') %0-6 times, from save data files %%%%%%%
load('Data_real_G1.mat')    %Time_TRAIL_G1 and Uncens_Length_G1_F0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uni=unique(Time_TRAIL_G1);   %from data
numb=20;
z=linspace(0,10,numb);       %from save Data

%%%%%% compare the y value (length G1) for each time of TRAIL addition (model and data)
negloglike=zeros(length(uni),numb);
for k=1:length(z)
    for i=1:length(uni)
        my_field = strcat('z',num2str(round(z(k)*100)));
        t_0_without_censored=t_0_model.(my_field);
        t_G1_new_without_censored=t_G1_model.(my_field);
        Fit_Model=fitdist(t_G1_new_without_censored(t_0_without_censored==uni(i)),'Lognormal');
        x_values = Uncens_Length_G1_F0(Time_TRAIL_G1==uni(i));
        yv = pdf(Fit_Model,x_values);
        negloglike(i,k)=-log(prod(yv));
    end
end
Summe_0=sum(negloglike,1);
%}

%save('Summe_G1.mat','Summe_0','Summe_1','Summe_2','Summe_3','Summe_4','Summe_5')
%only example here, Summe_G1 consists of 6 repeats.

load('Summe_G1.mat')
load('Data_real_G1.mat')    %Time_TRAIL_G1 and Uncens_Length_G1_F0
numb=20;
z=linspace(0,10,numb);   

mean_All=zeros(numb,1);
std_All=zeros(numb,1);
for i=1:numb
    mean_All(i)=mean([Summe_0(i),Summe_1(i),Summe_2(i),Summe_3(i),Summe_4(i),Summe_5(i)]);
    std_All(i)=std([Summe_0(i),Summe_1(i),Summe_2(i),Summe_3(i),Summe_4(i),Summe_5(i)]);
end

a=fit(z',mean_All,'poly5');
c=coeffvalues(a);
cd=polyder(c);
roots(cd);
z_a=linspace(0,10,100);   

x=3.6434;      %2.5 for HCT
mini_y= a.p1*x^5 + a.p2*x^4 + a.p3*x^3 + a.p4*x^2 + a.p5*x + a.p6;


figure()
plot(z_a,a(z_a),'-b', 'LineWidth',3)
hold on;
errorbar(z',mean_All,std_All,'ko')
grid on;
set(gca,'FontSize',20);
xlabel('z');
ylabel('-ln($\mathcal{L}$)','Interpreter','Latex');
xlim([2,6]);
h=legend(['exp=6, n=3e7'],'Location','northwest'); 


