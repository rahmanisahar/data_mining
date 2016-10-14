clear, close all,  clc
dir ='~/Desktop/project/data_mining/high_z_galaxies/results_SED/after_referee/kmeans/';
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

annv=annv';
%>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------
%>>>> KMeans
%>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------

k =4;
%[idx,C] = kmeans(annt,k);
opts = statset('Display','final');
[idx,C] = kmeans(annt,k,'Distance','cityblock',...
    'Replicates',5,'Options',opts);

if (size (annv(idx==1))==1) 
    cluster1 = annv(idx==1,:);
else
    cluster1 = median(annv(idx==1,:));
end
if (size (annv(idx==2))==1) 
    cluster2 = annv(idx==2,:);
else
    cluster2 = median(annv(idx==2,:));
end

if (size (annv(idx==3))==1) 
    cluster3 = annv(idx==3,:);
else
    cluster3 = median(annv(idx==3,:));
end


if (size (annv(idx==4))==1) 
    cluster4 = annv(idx==4,:);
else
    cluster4 = median(annv(idx==4,:));
end

SOM_Res_for_comp=[1,1,1,1,2,3,4,4,4,4,3,3];

if (size (annt(SOM_Res_for_comp==1))==1) 
    cluster_som1 = annt(SOM_Res_for_comp==1,:);
else
    cluster_som1 = median(annt(SOM_Res_for_comp==1,:));
end
if (size (annt(SOM_Res_for_comp==2))==1) 
    cluster_som2 = annt(SOM_Res_for_comp==2,:);
else
    cluster_som2 = median(annt(SOM_Res_for_comp==2,:));
end

if (size (annt(SOM_Res_for_comp==3))==1) 
    cluster_som3 = annt(SOM_Res_for_comp==3,:);
else
    cluster_som3 = median(annt(SOM_Res_for_comp==3,:));
end


if (size (annt(SOM_Res_for_comp==4))==1) 
    cluster_som4 = annt(SOM_Res_for_comp==4,:);
else
    cluster_som4 = median(annt(SOM_Res_for_comp==4,:));
end
SOM_in2=[1,1,1,1,1,2,2,2,2,2,2,2];


if (size (annt(SOM_in2==1))==1) 
    cluster_som1_in2 = annt(SOM_in2==1,:);
else
    cluster_som1_in2  = median(annt(SOM_in2==1,:));
end
if (size (annt(SOM_in2==2))==1) 
    cluster_som2_in2 = annt(SOM_in2==2,:);
else
    cluster_som2_in2 = median(annt(SOM_in2==2,:));
end



%cluster2 = median(annt(idx==2,:));
%cluster3 = median(annt(idx==3,:));
%cluster4 = median(annt(idx==4,:));
%>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------
%>>>>plotting
%>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------



figure(1);
plot(xx, annv)
ylabel('Normalized flux density')
xlabel('Wavelength($$\AA$$)','Interpreter', 'Latex' ) 

rep = strcat(dir,'all_spectrums_in_one_plot_142.fig');
 saveas(figure(1),rep,'fig')
 rep = strcat(dir,'all_spectrums_in_one_plot_142.pdf');
 saveas(figure(1),rep,'pdf')
 
 
% figure(2);
% plot(annt(idx==1,1),annt(idx==1,2),'r.','MarkerSize',12)
% hold on
% plot(annt(idx==2,1),annt(idx==2,2),'b.','MarkerSize',12)
% hold on
% plot(annt(idx==3,1),annt(idx==3,2),'y.','MarkerSize',12)
% plot(C(:,1),C(:,2),'kx',...
%      'MarkerSize',15,'LineWidth',3)
% legend('Cluster 1','Cluster 2','Centroids',...
%        'Location','NW')
% title 'Cluster Assignments and Centroids'
% hold off
% 
% rep = strcat(dir,'spec_vs_spec.fig');
%  saveas(figure(2),rep,'fig')
%  rep = strcat(dir,'spec_vs_spec.pdf');
%  saveas(figure(2),rep,'pdf')

figure(2);
plot(xx,cluster1,'r','MarkerSize',12)
hold on
plot(xx,cluster2,'b','MarkerSize',12)
hold on
plot(xx,cluster3,'g','MarkerSize',12)
hold on
plot(xx,cluster4,'Color',[0,0,0],'MarkerSize',12)
% plot(median(C(1,:)),median(C(2,:)),'kx',...
%      'MarkerSize',15,'LineWidth',3)
legend('Cluster 1','Cluster 2','Cluster 3',...
       'Cluster 4','NW')
ylabel('Normalized flux density')
xlabel('Wavelength($$\AA$$)','Interpreter', 'Latex' ) 
%title 'Cluster Assignments and Centroids'
hold off


rep = strcat(dir,'classified_group_in_',int2str(k),'cluster_142.fig');
saveas(figure(2),rep,'fig')
rep = strcat(dir,'classified_group_in_',int2str(k),'cluster_142.pdf');
saveas(figure(2),rep,'pdf')



figure(3);
plot(xx,cluster_som1,'r','MarkerSize',12)
hold on
plot(xx,cluster_som2,'b','MarkerSize',12)
hold on
plot(xx,cluster_som3,'g','MarkerSize',12)
hold on
plot(xx,cluster_som4,'Color',[0,0,0],'MarkerSize',12)
% plot(median(C(1,:)),median(C(2,:)),'kx',...
%      'MarkerSize',15,'LineWidth',3)
legend('SOM cluster 1','SOM cluster 2','SOM cluster 3',...
       'SOM cluster 4','NW')
ylabel('Normalized flux density')
xlabel('Wavelength($$\AA$$)','Interpreter', 'Latex' ) 
%title 'SOM cluster Assignments and Centroids'
hold off


rep = strcat(dir,'classified_group_in_',int2str(k),'cluster_som.fig');
saveas(figure(3),rep,'fig')
rep = strcat(dir,'classified_group_in_',int2str(k),'cluster_som.pdf');
saveas(figure(3),rep,'pdf')

figure(4);
plot(xx,cluster_som1_in2,'r','MarkerSize',12)
hold on
plot(xx,cluster_som2_in2,'b','MarkerSize',12)
% plot(median(C(1,:)),median(C(2,:)),'kx',...
%      'MarkerSize',15,'LineWidth',3)
legend('SOM cluster 1','SOM cluster 2','SOM cluster 3',...
       'SOM cluster 4','NW')
ylabel('Normalized flux density')
xlabel('Wavelength($$\AA$$)','Interpreter', 'Latex' ) 
%title 'SOM cluster Assignments and Centroids'
hold off


rep = strcat(dir,'classified_group_in_2cluster_som.fig');
saveas(figure(4),rep,'fig')
rep = strcat(dir,'classified_group_in_2cluster_som.pdf');
saveas(figure(4),rep,'pdf')

figure(10)
m1=0;
type=zeros(1,12);
% text(-10,10.2,'Wavelength($$\AA$$)','Interpreter', 'Latex') 

for h1=1:1:n_1
 for   h2=1:1:n_2
     check_s=CATv_1{h1,h2};
               size_ch=size(check_s);
               if (isnan(catm{h1,h2}) ~= 1) 
                if (size_ch(2) > 1)
    m1=m1+1;
    chi_s=zeros(1,12);
    a = catm{h1,h2};
    for n1=1:1:12
       chi = 0;
       for i = 1:1:630
           chi = chi+(a(i)-kinney(i,n1))^2/kinney(i,n1);
       end
       chi_s(n1)=chi;
    end   
    [dum,num_kinney]=min(chi_s(:));
   index=find(min(chi_s));
    m1=m1+1;
    subplot(n_1,n_2,m1)
    plot(xx,catm{h1,h2},'K')%,xlim([1 630])
    type(m1)=num_kinney;
    plot(xx,kinney(1:630,num_kinney),'r')%,xlim([1 630])
    strmax = ['Type = ',num2str(num_kinney)];
    text(3500,2,strmax,'HorizontalAlignment','right');
   ylabel('flux (arbitary unit)')
   xlabel('Wavelength($$\AA$$)','Interpreter', 'Latex' ) 
     end    
               end
 end
end
%text(-10,10.2,'Wavelength(\AA)')
SEDim = strcat(dir,'SED_total',n1st,'by',n2st,'.fig');

