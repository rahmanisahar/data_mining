%%%%% Classifying galaxies  using idx and centroid from K-means:
clc
%%%% Loading data:
dir ='~/Desktop/project/data_mining/high_z_galaxies/results_SED/after_referee/after_referee2/kmeans/';

load ~/Desktop/project/data_mining/high_z_galaxies/data/model_n.txt
catv=model_n(:,1:end)';

load ~/Desktop/project/data_mining/nearby_galaxies/H_paper_2012_som/xx.txt
x_ax=xx;

idx=csvread ...
    ('~/Desktop/project/data_mining/high_z_galaxies/results_SED/after_referee/after_referee2/kmeans/indexes_in_4_clusters.csv');
centroid=csvread ...
    ('~/Desktop/project/data_mining/high_z_galaxies/results_SED/after_referee/after_referee2/kmeans/centroid_in_4_clusters.csv');

sz_cat=size(catv);
sz_cnt =size(centroid);

dist_min=zeros(sz_cnt(1),1);
closest_cluster=zeros(sz_cat(1),1);

 for n=1:sz_cat(1)
     for h=1:sz_cnt(1)
        dist_min(h)=sqrt(sum((centroid(h,:)-catv(n,:)).^2));
     end
     closest_cluster(n)=find(dist_min==min(dist_min));
 end
 
medians_for_plotting=zeros(sz_cnt);
 
for h=1:sz_cnt(1)
	if (size (catv(closest_cluster==h))==1) 
    	medians_for_plotting(h,:) = catv(closest_cluster==h,:);
	else
    	medians_for_plotting(h,:) = median(catv(closest_cluster==h,:));
	end
end


csvwrite(strcat(dir,'closest_cluster_for_142_in_',num2str(sz_cnt(1)),'_k_means_clusters.csv'),closest_cluster)
csvwrite(strcat(dir,'medians_of_142_galaxies_in_each_cluster_in_',num2str(sz_cnt(1)),'_k_means_clusters.csv'),medians_for_plotting)

 
colorVec = hsv(sz_cnt(1));
figure(1);
  for h=1:sz_cnt(1)
    plot(xx,medians_for_plotting(h,:),'Color',colorVec(h,:),'MarkerSize',12)
    hold on
  end
  legend('Cluster 1','Cluster 2','Cluster 3',...
       'Cluster 4','NW')
  ylabel('Normalized flux density')
  xlabel('Wavelength($$\AA$$)','Interpreter', 'Latex' ) 
  hold off

  rep = strcat(dir,'classified_group_in_',int2str(sz_cnt(1)),'cluster_142.fig');
  saveas(figure(1),rep,'fig')
  rep = strcat(dir,'classified_group_in_',int2str(sz_cnt(1)),'cluster_142.pdf');
  saveas(figure(1),rep,'pdf')
 
 

