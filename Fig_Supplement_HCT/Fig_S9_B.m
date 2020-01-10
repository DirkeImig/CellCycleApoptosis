%%%%%% IntiQuan IQM Tools is needed: https://iqmtools.intiquan.com %%%%%%%%
%%%%%% exemplary parameter estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=12;

p=linspace(0.05,0.15,n);   
PAS=linspace(0.1,0.5,n);  


Sample_Data='Data_Model_Sample.mat';


model = IQMmodel('my_model.txt');
IQMmakeMEXmodel(model,'MEXmodel1');
MyModel='MEXmodel1';

Plot_p=[];
Plot_PAS=[];
C=[];
tic
for i=1:n
    for k=1:n
        Plot_p=[Plot_p,p(i)];
        Plot_PAS=[Plot_PAS,PAS(k)];
        C=[C,Myfun(p(i),PAS(k),Sample_Data,MyModel)];  
    end
end
toc


save('B2_ModelFit_0.mat','Plot_p','Plot_PAS','C')


