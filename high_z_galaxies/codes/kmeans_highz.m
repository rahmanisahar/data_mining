clear all, clc
load ~/Desktop/project/data_mining/high_z_galaxies/data/kinney_n.txt
cat=kinney_n(:,1:end);

load ~/Desktop/project/data_mining/high_z_galaxies/data/model_n.txt
catv=model_n(:,1:end);

load ~/Desktop/project/data_mining/nearby_galaxies/H_paper_2012_som/xx.txt
x_ax=xx;
%>>>>---------------------------------------------------------------------------------------------------
%>>>>  Inverse  'cat' to N x M  for normalization and other processes
kinney=cat;
annm=catv';
catv=catv';
cat=cat';        

%>>>>---------------------------------------------------------------------------------------------------
%>>>> fix NaN by takining avarage and assign a number to NaN and flag the NaN
%>>>> Network can not accept NaN or NULL

cat_fix=(cat);
catv_fix=(catv);

%>>>>---------------------------------------------------------------------------------------------------
%>>>> normalization data, select only one  
%>>>> mapminmax: mormalization between -1 and 1.
%>>>> mapstd: Gaussian normalization to sigma=1 and median=0
%>>>> cat_fix_norm= (cat_fix)  DO Nothing 

cat_fix_norm= (cat_fix);
%cat_fix_norm= mapminmax(cat_fix);
%cat_fix_norm= mapstd(cat_fix);

catv_fix_norm= (catv_fix);
%catv_fix_norm= mapminmax(catv_fix);
%catv_fix_norm= mapstd(catv_fix);

%>>>>---------------------------------------------------------------------------------------------------
%>>>> a  name change only to introduce t network
annt=cat_fix_norm';
annv=catv_fix_norm';
annt=annt';

[idx,C] = kmeans(annt,3);
% opts = statset('Display','final');
% [idx,C] = kmeans(annt,2,'Distance','cityblock',...
%     'Replicates',5,'Options',opts);
annt_tr=annt';
cluster1 = median(annt(idx==1,:));
cluster2 = median(annt(idx==2,:));
cluster3 = median(annt(idx==3,:));

figure;
plot(annt_tr(idx==1,:),annt_tr(idx==1,:),'r.','MarkerSize',12)
hold on
plot(annt_tr(idx==2,1),annt_tr(idx==2,2),'b.','MarkerSize',12)
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3)
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
title 'Cluster Assignments and Centroids'
hold off

figure(2);
plot(xx,cluster1,'r','MarkerSize',12)
hold on
plot(xx,cluster2,'b','MarkerSize',12)
hold on
plot(xx,cluster3,'g','MarkerSize',12)
% plot(median(C(1,:)),median(C(2,:)),'kx',...
%      'MarkerSize',15,'LineWidth',3)
% legend('Cluster 1','Cluster 2','Centroids',...
%        'Location','NW')
ylabel('Normalized flux density')
xlabel('Wavelength($$\AA$$)','Interpreter', 'Latex' ) 
title 'Cluster Assignments and Centroids'
hold off