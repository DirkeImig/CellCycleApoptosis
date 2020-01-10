n=12;

p=linspace(0.1,0.18,n);   
PAS=linspace(0.41,0.58,n);  

Sample_Data='Data_Model_Sample_0.mat';

Plot_p=[];
Plot_PAS=[];
C=[];
tic
for i=1:n
    for k=1:n
        Plot_p=[Plot_p,p(i)];
        Plot_PAS=[Plot_PAS,PAS(k)];
        C=[C,myfun(p(i),PAS(k),Sample_Data)];  
    end
end
toc


%save('ModelFit_0.mat','Plot_p','Plot_PAS','C')




