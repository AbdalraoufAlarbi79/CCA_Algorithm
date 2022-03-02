clear all; close all; 
k1=6;
%[tr_inp1,tr_out1,te_inp1,lab] = Gen_Data_NN(1);

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
     data_org1(x,:)= tr_inp1(i,:); x=x+1; 
    else
     data_org2(y,:)= tr_inp1(i,:); y=y+1;
    end;
end;
dim1= size(data_org1,1);  
dim2= size(data_org2,1); 
for ite=1:2
[IDX1] = kmedoids (data_org1,k1); % kmedoids    kmeans
[IDX2] = kmedoids (data_org2,k1); 
x=1;y=1;z=1; w=1;h=1;k=1;
C1_k1=[]; C1_k2=[];  C1_k3=[]; C1_k4=[]; C1_k5=[];  C1_k6=[]; 

C2_k1=[]; C2_k2=[];  C2_k3=[]; C2_k4=[]; C2_k5=[];  C2_k6=[];
data_org11=[];data_org12=[]; data_org13=[]; data_org14=[];data_org15=[]; data_org16=[];

data_org21=[];data_org22=[]; data_org23=[];  data_org24=[];data_org25=[]; data_org26=[];

for i=1:size(data_org1,1)
   if IDX1(i,1)==1;
     data_org11(x,:)= data_org1(i,:); C1_k1(x,3)=i;  x=x+1;
   elseif IDX1(i,1)==2;
     data_org12(y,:)= data_org1(i,:); C1_k2(y,3)=i;  y=y+1;
   elseif IDX1(i,1)==3;
     data_org13(z,:)= data_org1(i,:); C1_k3(z,3)=i;  z=z+1;
  elseif IDX1(i,1)==4;
    data_org14(w,:)= data_org1(i,:); C1_k4(w,3)=i;  w=w+1;
  elseif IDX1(i,1)==5;
    data_org15(h,:)= data_org1(i,:); C1_k5(h,3)=i;  h=h+1;
  elseif IDX1(i,1)==6;
    data_org16(k,:)= data_org1(i,:); C1_k6(k,3)=i;  k=k+1;
  end;
end;
x=1;y=1;z=1; w=1;h=1;k=1;

for i=1:size(data_org2,1)
    if IDX2(i,1)==1;
       data_org21(x,:)= data_org2(i,:); C2_k1(x,3)=i; x=x+1;
    elseif IDX2(i,1)==2;
       data_org22(y,:)= data_org2(i,:); C2_k2(y,3)=i; y=y+1;
    elseif IDX2(i,1)==3;       
       data_org23(z,:)= data_org2(i,:); C2_k3(z,3)=i; z=z+1;        
    elseif IDX2(i,1)==4;
       data_org24(w,:)= data_org2(i,:); C2_k4(w,3)=i; w=w+1;
    elseif IDX2(i,1)==5;
       data_org25(h,:)= data_org2(i,:); C2_k5(h,3)=i; h=h+1;
    elseif IDX2(i,1)==6;       
       data_org26(k,:)= data_org2(i,:); C2_k6(k,3)=i; k=k+1;        
    end;
end;
DM11=Dis_F(data_org11);DM12=Dis_F(data_org12); DM13=Dis_F(data_org13); DM14=Dis_F(data_org14);DM15=Dis_F(data_org15); DM16=Dis_F(data_org16);  
  
DM21=Dis_F(data_org21);DM22=Dis_F(data_org22); DM23=Dis_F(data_org23); DM24=Dis_F(data_org24);DM25=Dis_F(data_org25); DM26=Dis_F(data_org26);  

dim11= size(data_org11,1); dim12= size(data_org12,1);  dim13= size(data_org13,1); dim14= size(data_org14,1); dim15= size(data_org15,1);  dim16= size(data_org16,1);

dim21= size(data_org21,1); dim22= size(data_org22,1);  dim23= size(data_org23,1); dim24= size(data_org24,1); dim25= size(data_org25,1);  dim26= size(data_org26,1); 
for i=1:dim11;
    C1_k1(i,1)=sum(DM11(i,:));
 end;
 for i=1:dim12;
    C1_k2(i,1)=sum(DM12(i,:));
 end;
 for i=1:dim13;
    C1_k3(i,1)=sum(DM13(i,:));
 end;
  for i=1:dim14;
    C1_k4(i,1)=sum(DM14(i,:));
 end;
 for i=1:dim15;
    C1_k5(i,1)=sum(DM15(i,:));
 end;
 for i=1:dim16;
    C1_k6(i,1)=sum(DM16(i,:));
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
 
 for i=1:dim24;
    C2_k4(i,1)=sum(DM24(i,:));
 end;
 for i=1:dim25;
    C2_k5(i,1)=sum(DM25(i,:));
 end;
 for i=1:dim26;
    C2_k6(i,1)=sum(DM26(i,:));
 end;
 
[min11_c,pos11_c] = min(C1_k1(:,1));%center
[min12_c,pos12_c] = min(C1_k2(:,1));%center
[min13_c,pos13_c] = min(C1_k3(:,1));%center
[min14_c,pos14_c] = min(C1_k4(:,1));%center
[min15_c,pos15_c] = min(C1_k5(:,1));%center
[min16_c,pos16_c] = min(C1_k6(:,1));%center

[min21_c,pos21_c] = min(C2_k1(:,1));%center
[min22_c,pos22_c] = min(C2_k2(:,1));%center
[min23_c,pos23_c] = min(C2_k3(:,1));%center
[min24_c,pos24_c] = min(C2_k4(:,1));%center
[min25_c,pos25_c] = min(C2_k5(:,1));%center
[min26_c,pos26_c] = min(C2_k6(:,1));%center

for i=1:size(te_inp1,1)
 %center
 dis1(1,1)  =   pdist2(te_inp1(i,:),data_org1(C1_k1(pos11_c,3),:),'euclidean'); 
 dis1(2,1)  =   pdist2(te_inp1(i,:),data_org1(C1_k2(pos12_c,3),:),'euclidean');
 dis1(3,1)  =   pdist2(te_inp1(i,:),data_org1(C1_k3(pos13_c,3),:),'euclidean');
 dis1(4,1)  =   pdist2(te_inp1(i,:),data_org1(C1_k4(pos14_c,3),:),'euclidean'); 
 dis1(5,1)  =   pdist2(te_inp1(i,:),data_org1(C1_k5(pos15_c,3),:),'euclidean');
 dis1(6,1)  =   pdist2(te_inp1(i,:),data_org1(C1_k6(pos16_c,3),:),'euclidean');
 
 dis1(7,1)  =   pdist2(te_inp1(i,:),data_org2(C2_k1(pos21_c,3),:),'euclidean'); 
 dis1(8,1)  =   pdist2(te_inp1(i,:),data_org2(C2_k2(pos22_c,3),:),'euclidean');
 dis1(9,1)  =   pdist2(te_inp1(i,:),data_org2(C2_k3(pos23_c,3),:),'euclidean');
 dis1(10,1)  =  pdist2(te_inp1(i,:),data_org2(C2_k4(pos24_c,3),:),'euclidean'); 
 dis1(11,1)  =  pdist2(te_inp1(i,:),data_org2(C2_k5(pos25_c,3),:),'euclidean');
 dis1(12,1)  =  pdist2(te_inp1(i,:),data_org2(C2_k6(pos26_c,3),:),'euclidean');
 [min1,pos1]= min(dis1(:,1));%GC
 if pos1<7    
   lab(i,2)=1;
 else
   lab(i,2)=0;
 end;
 
 diff_ce(:,1)=abs(lab(:,1)-lab(:,2)) ; 
num1_ce=sum(diff_ce==0); num2_ce=sum(diff_ce==1);
accu_ce(ite,1)=num1_ce/size(lab,1);
end;
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

figure;    hold on 
%plot(data_org11,'bo')
for i=1:size(data_org11)
     plot(data_org11(i, 1),     data_org11(i,2),'g*');
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
