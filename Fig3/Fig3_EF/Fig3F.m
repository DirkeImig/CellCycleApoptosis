%%%%%%%%%%% parameter estimation - exemplary %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% EXAMPLE: Calculation of Likelihood %%%%%%%%%%%%%%%%%%%%%%%%%%
%{
load('Data_Model_SG2M_0.mat')  
load('Data_real_SG2M.mat')   %Time_TRAIL_SG2M and Uncens_Length_SG2M_F0

uni=unique(Time_TRAIL_SG2M);    %from data
numb=20;
z=linspace(0,10,numb); %from save Data


%%%%%% compare the y value (length SG2M) for each time of TRAIL addition (model and data)
negloglike=zeros(length(uni),numb);

for k=1:length(z)
    for i=1:length(uni)
        my_field = strcat('z',num2str(round(z(k)*100)));
        t_0_without_censored=t_0_model.(my_field);
        t_SG2M_new_without_censored=t_SG2M_model.(my_field);
        Fit_Model=fitdist(t_SG2M_new_without_censored(t_0_without_censored==uni(i)),'Lognormal');
        x_values = Uncens_Length_SG2M_F0(Time_TRAIL_SG2M==uni(i));
        yv = pdf(Fit_Model,x_values);
        negloglike(i,k)=-log(prod(yv));
    end
end
Summe_0=sum(negloglike,1);
%}

%save('Summe_SG2M.mat','Summe_0','Summe_1','Summe_2','Summe_3','Summe_4','Summe_5')
load('Data_real_SG2M.mat')   %Time_TRAIL_SG2M and Uncens_Length_SG2M_F0
load('Summe_SG2M.mat')
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

x=3.4086;   %1.9508 for HCT
mini_y= a.p1*x^5 + a.p2*x^4 + a.p3*x^3 + a.p4*x^2 + a.p5*x + a.p6;


figure()
plot(z_a,a(z_a),'-b', 'LineWidth',3)
hold on;
errorbar(z',mean_All,std_All,'ko')
grid on;
set(gca,'FontSize',20);
xlabel('z');
ylabel('-ln($\mathcal{L}$)','Interpreter','Latex');
xlim([2,5]);
h=legend(['exp=6, n=3e7'],'Location','northwest'); 


