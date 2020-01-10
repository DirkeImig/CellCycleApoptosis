%%%%%%%%%%%%% create heatmap for different approaches and z-values %%%%%%%%
%%%%%%%%%%%%%%%%%%%%% S/G2/M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=6;
n=100;

%%%%%%%%%%%%%%% z=0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z=0;
Z0=zeros(m,2);
B=zeros(m,n);
for j=1:n
    B(:,j)=myfun_SG2M(z);
end

for i=1:m
    Z0(i,1)=mean(B(i,:));
    Z0(i,2)=std(B(i,:));
end



%%%%%%%%%%%%%%% z=2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z=2;
Z2=zeros(m,2);
B=zeros(m,n);
for j=1:n
    B(:,j)=myfun_SG2M(z);
end

for i=1:m
    Z2(i,1)=mean(B(i,:));
    Z2(i,2)=std(B(i,:));
end




%%%%%%%%%%%%%%% z=4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z=4;
Z4=zeros(m,2);
B=zeros(m,n);
for j=1:n
    B(:,j)=myfun_SG2M(z);
end

for i=1:m
    Z4(i,1)=mean(B(i,:));
    Z4(i,2)=std(B(i,:));
end
     



%%%%%%%%%%%%%%% z=6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z=6;
Z6=zeros(m,2);
B=zeros(m,n);
for j=1:n
    B(:,j)=myfun_SG2M(z);
end

for i=1:m
    Z6(i,1)=mean(B(i,:));
    Z6(i,2)=std(B(i,:));
end



%%%%%%%%%%%%%%% z=8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z=8;
Z8=zeros(m,2);
B=zeros(m,n);
for j=1:n
    B(:,j)=myfun_SG2M(z);
end

for i=1:m
    Z8(i,1)=mean(B(i,:));
    Z8(i,2)=std(B(i,:));
end


Plot_A=[Z0(1,1), Z2(1,1),Z4(1,1),Z6(1,1),Z8(1,1); ...
    Z0(3,1),Z2(3,1),Z4(3,1),Z6(3,1),Z8(3,1);...
    Z0(5,1), Z2(5,1),Z4(5,1),Z6(5,1),Z8(5,1)];
Plot_A=round(100.*Plot_A)./100;


map = [0.8,1,1
    0.5,0.7,1
    0.4 0.6 1
    0.3 0.5 1
    0.2,0.4,1
    ];

figure()
xvalues = {'0','2','4','6','8'};
yvalues = {'I','II','III'};
h=heatmap(xvalues,yvalues,Plot_A, 'FontSize',26,'colormap', map);
