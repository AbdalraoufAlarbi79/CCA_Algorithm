clear all; close all; 

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
 
x=1;y=1;result1=[];result2=[];
for i=1:size(tr_inp1,1)
    if tr_out1(i,1)==1;
     DATA_ORG1(x,:)= tr_inp1(i,:); result1(x,3)=i; x=x+1;
    else
     DATA_ORG2(y,:)= tr_inp1(i,:); result2(y,3)=i;  y=y+1;
    end;
end;
DM1=Dis_F(DATA_ORG1);  
dim1= size(DATA_ORG1,1);  
DM2=Dis_F(DATA_ORG2); 
dim2= size(DATA_ORG2,1);  

for i=1:dim1;
 
   result1(i,1)=sum(DM1(i,:));
 end;
for i=1:dim2;
  
   result2(i,1)=sum(DM2(i,:));
end;
[min11,pos11]   = min(result1(:,1));%center
result1(pos11,1)= max(result1(:,1))+1;
[min12,pos12]   = min(result1(:,1));%center
result1(pos12,1)= max(result1(:,1))+1;
[min13,pos13]   = min(result1(:,1));%center
result1(pos13,1)= max(result1(:,1))+1;
[min14,pos14]   = min(result1(:,1));%center
result1(pos14,1)= max(result1(:,1))+1;
[min15,pos15]    = min(result1(:,1));%center


[min21,pos21]   = min(result2(:,1));%center
result2(pos21,1)= max(result2(:,1))+1;
[min22,pos22]   = min(result2(:,1));%center
result2(pos22,1)= max(result2(:,1))+1;
[min23,pos23]   = min(result2(:,1));%center
result2(pos23,1)= max(result2(:,1))+1;
[min24,pos24]   = min(result2(:,1));%center
result2(pos24,1)= max(result2(:,1))+1;
[min25,pos25]   = min(result2(:,1));%center

for i=1:size(te_inp1,1)
 %center
 dis11 = pdist2(te_inp1(i,:),tr_inp1(result1(pos11,3),:),'euclidean'); 
 dis12 = pdist2(te_inp1(i,:),tr_inp1(result1(pos12,3),:),'euclidean'); 
 dis13 = pdist2(te_inp1(i,:),tr_inp1(result1(pos13,3),:),'euclidean'); 
 dis14 = pdist2(te_inp1(i,:),tr_inp1(result1(pos14,3),:),'euclidean'); 
 dis15 = pdist2(te_inp1(i,:),tr_inp1(result1(pos15,3),:),'euclidean'); 

 dis21 = pdist2(te_inp1(i,:),tr_inp1(result2(pos21,3),:),'euclidean');
 dis22 = pdist2(te_inp1(i,:),tr_inp1(result2(pos22,3),:),'euclidean');
 dis23 = pdist2(te_inp1(i,:),tr_inp1(result2(pos23,3),:),'euclidean');
 dis24 = pdist2(te_inp1(i,:),tr_inp1(result2(pos24,3),:),'euclidean');
 dis25 = pdist2(te_inp1(i,:),tr_inp1(result2(pos25,3),:),'euclidean');
 dis1=(dis11+dis12+dis13+dis14+dis15)/5
 dis2=(dis21+dis22+dis23+dis24+dis25)/5
 
 if dis1<dis2
    lab(i,2)=1;
 else
      lab(i,2)=0;
 end;
 end
diff_ce(:,1)=abs(lab(:,1)-lab(:,2)) ; 
num1_ce=sum(diff_ce==0); num2_ce=sum(diff_ce==1);
accu_ce=num1_ce/size(lab,1);

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
for i=1:size(DATA_ORG1)
     plot(DATA_ORG1(i, 1),     DATA_ORG1(i,2),'b-o');
end;
for i=1:size(DATA_ORG2)
     plot(DATA_ORG2(i, 1),     DATA_ORG2(i,2),'m-O');
end;

plot(tr_inp1(result1(pos11,3),1), tr_inp1(result1(pos11,3),2),'r-*');
plot(tr_inp1(result2(pos21,3),1), tr_inp1(result2(pos21,3),2),'b-*');

%plot(tr_inp1(result1(pos12,3),1), tr_inp1(result1(pos12,3),2),'g-..');
%plot(tr_inp1(result2(pos22,3),1), tr_inp1(result2(pos22,3),2),'g-..');

%}