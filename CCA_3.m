%
clear all; close all; 
k=3;% this means the iteration equal 2 

%[tr_inp1,tr_out1,te_inp1,lab] = Gen_Data_NN(1);
% tr_inp1 is represent the training data input (size 7738 * 30)
% tr_out1 is represent the training data output(size 7738 * 1)
%te_inp1_1 is represent the testing data input( size 3317 * 30)
%te_out1 is represent the testing data output( size 3317 * 1)

A=xlsread('D:\orgenal\ds1\ds1_1.xlsx','Sheet1'); %real
[m,n] = size(A) ;
P = 0.80 ;
idx = randperm(m)  ;
Training = A(idx(1:round(P*m)),:) ; 
Testing = A(idx(round(P*m)+1:end),:) ;
tr_inp1=Training(:,1:end-1);
te_inp1=Testing(:,1:end-1)
tr_out1=Training(:,end);
lab=Testing(:,end);


x=1;y=1;
% in this loop we will divided the training datasets to Legitimate and phishing
for i=1:size(tr_inp1,1)% we will read the matrix for the input data from 1 to the end 
    if tr_out1(i,1)==1;% we will read and make examination of 
%the first number in the output matrix to classify the phishing from Legitimate
     data_org1(x,:)= tr_inp1(i,:); x=x+1; % we classify the legitimate result here
    else
     data_org2(y,:)= tr_inp1(i,:); y=y+1;% we classify the phishing result here
    end;
end;
%here we will get the size of the tow matrix



dim1= size(data_org1,1);% legitimate matrix   
dim2= size(data_org2,1);% phishing matrix

% here we will start our model where we will start first iteration 
%-------------------------------------------------------------------
%-------------------------------------------------------------------
%-------------------------------------------------------------------
%-------------------------------------------------------------------
%-------------------------------------------------------------------
%-------------------------------------------------------------------
%-------------------------------------------------------------------

for ite=1:2
    % here we will get the tow matrices into the Kmaens cluster to get K
    % cluster as needed
    
[IDX1] = kmeans (data_org1,k); % kmedoids    kmeans
[IDX2] = kmeans (data_org2,k); 

x=1;y=1;z=1;
C1_k1=[]; C1_k2=[]; C1_k3=[];% this virible represent the Class {1} and Cluster {1,2,..}
C2_k1=[]; C2_k2=[]; C2_k3=[];% this virible represent the Class {2} and Cluster {1,2,..}
%set matrices for the clusteres 

data_org11=[];data_org12=[];data_org13=[];
data_org21=[];data_org22=[];data_org23=[];

%this loop is for divide the clusteres into it's matriex for the legitimate class 
for i=1:size(data_org1,1)
    if IDX1(i,1)==1;
     data_org11(x,:)= data_org1(i,:); C1_k1(x,3)=i;  x=x+1;
    elseif IDX1(i,1)==2;
     data_org12(y,:)= data_org1(i,:); C1_k2(y,3)=i;  y=y+1;
    else
     data_org13(z,:)= data_org1(i,:); C1_k3(z,3)=i;  z=z+1;
  
    end;
end;


x=1;y=1;z=1;
%this loop is for divide the clusteres into it's matriex for the phishing class 
for i=1:size(data_org2,1)
    if IDX2(i,1)==1;
       data_org21(x,:)= data_org2(i,:); C2_k1(x,3)=i; x=x+1;
    elseif IDX2(i,1)==2; 
       data_org22(y,:)= data_org2(i,:); C2_k2(y,3)=i; y=y+1;
    else
       data_org23(z,:)= data_org2(i,:); C2_k3(z,3)=i; z=z+1; 
    end;
end;

% here we will use the function name Dis_F where its giving back the distance matrix for each 
%cluster
DM11=Dis_F(data_org11);DM12=Dis_F(data_org12);DM13=Dis_F(data_org13);  
DM21=Dis_F(data_org21);DM22=Dis_F(data_org22);DM23=Dis_F(data_org23);

% here we will get the size for each cluster 
dim11= size(data_org11,1); dim12= size(data_org12,1);dim13= size(data_org13,1);  
dim21= size(data_org21,1); dim22= size(data_org22,1);dim23= size(data_org23,1);

% in this loop we will get the sumation for each ROw in it's cluster for
% Legitimate class and it's cluster
for i=1:dim11;
   C1_k1(i,1)=sum(DM11(i,:));
 end;
 for i=1:dim12;
   C1_k2(i,1)=sum(DM12(i,:));
 end;
 for i=1:dim13;
   C1_k3(i,1)=sum(DM13(i,:));
 end;
 
 
 for i=1:dim21;
   C2_k1(i,1)=sum(DM21(i,:));
 end;
 for i=1:dim22;
   C2_k2(i,1)=sum(DM22(i,:));
 end;
 for i=1:dim23;
   C2_k3(i,1)=sum(DM23(i,:));
 end;
 
 % in this step we will used MIN function to get the smallest  
[min11_c,pos11_c] = min(C1_k1(:,1));%center
[min12_c,pos12_c] = min(C1_k2(:,1));%center
[min13_c,pos13_c] = min(C1_k3(:,1));%center




[min21_c,pos21_c] = min(C2_k1(:,1));%center
[min22_c,pos22_c] = min(C2_k2(:,1));%center
[min23_c,pos23_c] = min(C2_k3(:,1));%center


for i=1:size(te_inp1,1)
 %center
 dis1(1,1)  =   pdist2(te_inp1(i,:),data_org1(C1_k1(pos11_c,3),:),'euclidean'); 
 dis1(2,1)  =   pdist2(te_inp1(i,:),data_org1(C1_k2(pos12_c,3),:),'euclidean');
 dis1(3,1)  =   pdist2(te_inp1(i,:),data_org1(C1_k3(pos13_c,3),:),'euclidean');

 
 
 dis1(4,1)  =   pdist2(te_inp1(i,:),data_org2(C2_k1(pos21_c,3),:),'euclidean'); 
 dis1(5,1)  =   pdist2(te_inp1(i,:),data_org2(C2_k2(pos22_c,3),:),'euclidean');
 dis1(6,1)  =   pdist2(te_inp1(i,:),data_org2(C2_k3(pos23_c,3),:),'euclidean');

 
 [min1,pos1]= min(dis1(:,1));%GC
 if pos1<4    
   lab(i,2)=1;
 else
   lab(i,2)=0;
 end;
  end;
  %{
  [min11_c,pos11_c] = min(C1_k1(:,1));%center
[min12_c,pos12_c] = min(C1_k2(:,1));%center
[min21_c,pos21_c] = min(C2_k1(:,1));%center
[min22_c,pos22_c] = min(C2_k2(:,1));%center
  accu_ce(1,2)=
    accu_ce(1,2)=
  accu_ce(1,2)=
  accu_ce(1,2)=
%}
 diff_ce(:,1)=abs(lab(:,1)-lab(:,2)) ; 
num1_ce=sum(diff_ce==0); num2_ce=sum(diff_ce==1);
accu_ce(ite,1)=num1_ce/size(lab,1);
end;
acc1_ce=max(accu_ce(:,1));

predicted = lab(:,1);
actual =lab(:,2);
confMatx = confusionmat(actual,predicted);
confMatx = confMatx';
diagonal = diag(confMatx);
sum_of_rows = sum(confMatx , 2);
precision = diagonal ./ sum_of_rows ;
overall_precision = mean ( precision);
sum_of_columns = sum ( confMatx , 1);
recall = diagonal ./ sum_of_columns';
overall_recall = mean(recall);
f1_score = 2*((overall_precision * overall_recall)/(overall_precision + overall_recall));

%
figure;    hold on 
for i=1:size(data_org11)
     plot(data_org11(i, 1),     data_org11(i,2),'b-O');
end;
for i=1:size(data_org12)
     plot(data_org12(i, 1),     data_org12(i,2),'g-O');
end;
for i=1:size(data_org21)
     plot(data_org21(i, 1),     data_org21(i,2),'r-O');
end;
for i=1:size(data_org22)
     plot(data_org22(i, 1),     data_org22(i,2),'m-O');
end;

plot(data_org1(C1_k1(pos11_c,3),1), data_org1(C1_k1(pos11_c,3),2),'r-+');
plot(data_org1(C1_k2(pos12_c,3),1), data_org1(C1_k2(pos12_c,3),2),'r-+');
plot(data_org2(C2_k1(pos21_c,3),1), data_org2(C2_k1(pos21_c,3),2),'k-*');
plot(data_org2(C2_k2(pos22_c,3),1), data_org2(C2_k2(pos22_c,3),2),'k-*');


%}

