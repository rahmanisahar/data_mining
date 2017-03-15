clear, close all,  clc

dir ='~/Desktop/project/data_mining/high_z_galaxies/results_SED/after_referee/after_referee2/kmeans/';

load ~/Desktop/project/data_mining/high_z_galaxies/data/kinney_n.txt
cat=kinney_n(:,1:end)'; %%%>>>>Kinny_data

load ~/Desktop/project/data_mining/high_z_galaxies/data/model_n.txt
catv=model_n(:,1:end)'; %>>>>> H_paper_data

load ~/Desktop/project/data_mining/nearby_galaxies/H_paper_2012_som/xx.txt
x_ax=xx;

%>>>>---------------------------------------------------------------------------------------------------
%>>>> normalization data, select only one  
%>>>> mapminmax: mormalization between -1 and 1.
%>>>> mapstd: Gaussian normalization to sigma=1 and median=0
%>>>> cat_fix_norm= (cat_fix)  DO Nothing 

cat_fix_norm= (cat); %%%>>>>Kinny_data
%cat_fix_norm= mapminmax(cat_fix);
%cat_fix_norm= mapstd(cat_fix);

catv_fix_norm= (catv); %>>>>> H_paper_data
%catv_fix_norm= mapminmax(catv_fix);
%catv_fix_norm= mapstd(catv_fix);

%>>>>---------------------------------------------------------------------------------------------------
%>>>> a  name change only to introduce t network
annt=cat_fix_norm'; %%%>>>>Kinny_data
annv=catv_fix_norm'; %>>>>> H_paper_data
annt=annt'; 

annv=annv';  %>>>>> H_paper_data
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
[idx,Centroid] = kmeans(annt,k,'Distance','cityblock',...
    'Replicates',5,'Options',opts);
csvwrite(strcat(dir,'indexes_in_',num2str(k),'_clusters.csv'),idx)
csvwrite(strcat(dir,'centroid_in_',num2str(k),'_clusters.csv'),Centroid)
%>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------
%medians_for_plotting_kmeans
%>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------
%%%%if it is not equal Centroid, there is a problem!!

 sz_cnt =size(Centroid);
% medians_for_plotting=zeros(sz_cnt);
% 
% for h=1:k
%   if (size (annt(idx==h))==1) 
%       medians_for_plotting(h,:) = annt(idx==h,:);
%   else
%       medians_for_plotting(h,:) = median(annt(idx==h,:));
%   end
% end
% %>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------
%medians_for_plotting_SOM
%>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------
%>>>>---------------------------------------------------------------------------------------------------


SOM_Res_for_comp=[1,1,1,1,2,3,4,4,4,4,3,3]';

cluster_som_in4=zeros(4,sz_cnt(2));

for h=1:k
  if (size (annt(SOM_Res_for_comp==h))==1) 
      cluster_som_in4(h,:) = annt(SOM_Res_for_comp==h,:);
  else
      cluster_som_in4(h,:) = median(annt(SOM_Res_for_comp==h,:));
  end
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



% figure(1);
% plot(xx, annv)
% ylabel('Normalized flux density')
% xlabel('Wavelength($$\AA$$)','Interpreter', 'Latex' ) 
% 
% rep = strcat(dir,'all_spectrums_in_one_plot_142.fig');
%  saveas(figure(1),rep,'fig')
%  rep = strcat(dir,'all_spectrums_in_one_plot_142.pdf');
%  saveas(figure(1),rep,'pdf')
%>>>>> Pllotting_original_cluster
 colorVec = hsv(k);
figure(1);
  for h=1:k
    plot(xx,Centroid(h,:),'Color',colorVec(h,:),'MarkerSize',12)
    hold on
  end
  legend('K-means cluster 1','K-means cluster 2','K-means cluster 3',...
       'K-means cluster 4','NW')
  ylabel('Normalized flux density')
  xlabel('Wavelength($$\AA$$)','Interpreter', 'Latex' ) 
  hold off

  rep = strcat(dir,'classified_group_in_',int2str(k),'cluster_kinney.fig');
  saveas(figure(1),rep,'fig')
  rep = strcat(dir,'classified_group_in_',int2str(k),'cluster_kinney.pdf');
  saveas(figure(1),rep,'pdf')
 

% colorVec = hsv(k);
% figure(2);
%   for h=1:k
%     plot(xx,medians_for_plotting(h),'Color',colorVec(h,:),'MarkerSize',12)
%     hold on
%   end
%   legend('Cluster 1','Cluster 2','Cluster 3',...
%        'Cluster 4','NW')
%   ylabel('Normalized flux density')
%   xlabel('Wavelength($$\AA$$)','Interpreter', 'Latex' ) 
%   hold off

%   rep = strcat(dir,'classified_group_in_',int2str(k),'cluster_142.fig');
%   saveas(figure(2),rep,'fig')
%   rep = strcat(dir,'classified_group_in_',int2str(k),'cluster_142.pdf');
%   saveas(figure(2),rep,'pdf')



figure(3);
  for h=1:k
    plot(xx,cluster_som_in4(h,:),'Color',colorVec(h,:),'MarkerSize',12)
    hold on
  end
  legend('SOM cluster 1','SOM cluster 2','SOM cluster 3',...
         'SOM cluster 4','NW')
  ylabel('Normalized flux density')
    xlabel('Wavelength($$\AA$$)','Interpreter', 'Latex' ) 
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
  hold off

  rep = strcat(dir,'classified_group_in_2cluster_som.fig');
  saveas(figure(4),rep,'fig')
  rep = strcat(dir,'classified_group_in_2cluster_som.pdf');
  saveas(figure(4),rep,'pdf')


