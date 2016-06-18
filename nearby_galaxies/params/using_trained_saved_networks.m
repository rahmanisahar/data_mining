close all , clear, clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This program reads the previously generated networks. 
% Applies them to data in order to create plots in grey cycle clours.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reading data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 n_1=10; %number of raws in network
 n_2=10; %number of columns in network
%%% NETWORK NAME
net10by10=load('~/Desktop/project/data_mining/nearby_galaxies/SOM/org_colours2d/subsets/without_reg1/net_without_reg1_all.mat'); %Load network
net=net10by10.net;

%%%% Trining file
cat_t = csvread('~/Desktop/project/data_mining/nearby_galaxies/m31/ascii_tables/derived_ones_with_mean_per_arcsec.csv',2,1,[2,1,10,26]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preparing the training file for the net
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
catt=cat_t';
catt_fix=fixunknowns(catt);
y_min = -1;
y_max = 1;
sz = size(catt);
catt_min = min(catt')' * ones(1,sz(2));
catt_max = max(catt')' * ones(1,sz(2));
catt_fix_norm = (y_max - y_min) * (catt - catt_min) ./ (catt_max - catt_min) + y_min;
annt=catt_fix_norm; %changing namme to introduce to network
sz_t=size(annt); %finding size of original data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Giving data to our network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_t=sim(net, annt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Saving informations for training set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k1=1:n_1*n_2
 at{k1}=find(sim_t(k1,:)==1);
end


for k1=1:n_1*n_2
inpt{k1}=annt(:,at{k1});
end
 


m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
    m1=m1+1;
    
    Tabt_1{h1,h2}=at{m1};
    
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting and saving tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
    plotsomnd_shar(net,annt) %MATLAB som built-in SOM plots; shows distance between each neuran' Neighbours
figure(2)
    plotsomhits_sahar(net,annt) %MATLAB som built-in SOM plots; shows density of each neurans
figure(3)
    plotsomplanes_sahar(net) %MATLAB som built-in SOM plots; shows density of each neurans
figure(4)
    plotsom_sahar(net,annt) %MATLAB som built-in SOM plots; shows distance between each neuran' Neighbours


    saveas(figure(1),'~/Desktop/project/data_mining/nearby_galaxies/SOM/grey_cycle/2d/subset1_dist.fig','fig')
    saveas(figure(1),'~/Desktop/project/data_mining/nearby_galaxies/SOM/grey_cycle/2d/subset1_dist.png','png')
    saveas(figure(4),'~/Desktop/project/data_mining/nearby_galaxies/SOM/grey_cycle/2d/subset1_dist_with_hits_t.fig','fig')
    saveas(figure(4),'~/Desktop/project/data_mining/nearby_galaxies/SOM/grey_cycle/2d/subset1_dist_with_hits_t.png','png')
    saveas(figure(3),'~/Desktop/project/data_mining/nearby_galaxies/SOM/grey_cycle/2d/subset1_planes.fig','fig')
    saveas(figure(3),'~/Desktop/project/data_mining/nearby_galaxies/SOM/grey_cycle/2d/subset1_planes.png','png')
   close all