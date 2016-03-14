

net10by10=load('~/Desktop/project/data_mining/SOM/derived/per_arcsec_sq/subsets/net_subset9.mat'); %Load network
net=net10by10.net;
%catv = csvread('~/Desktop/project/data_mining/m31/ascii_tables/derived_ones_with_mean_per_arcsec.csv',1,3,[1,3,10,26]); %load data
catv = csvread('~/Desktop/project/data_mining/m101/ascii_table/m101_total_with_mean_per_arcsecsq_same_as_derived_ones.csv',1,1);
%load m101_va.txt
%catv=m101_va(:,1:end);
catv=catv';
catv_fix=fixunknowns(catv);
[catv1_fix_norm,~]= mapminmax(catv);
y_min = -1;
y_max = 1;
sz = size(catv);
catv_min = min(catv')' * ones(1,sz(2));
catv_max = max(catv')' * ones(1,sz(2));
catv_fix_norm = (y_max - y_min) * (catv - catv_min) ./ (catv_max - catv_min) + y_min;
annv=catv_fix_norm; %changing namme to introduce to network
sz=size(annv); %finding size of original data
 n_1=10;
 n_2=10;


%%%%%%%%%%%%%
%%%%%%%%%%%%%
%Giving data to our network
%%%%%%%%%%%%%
%%%%%%%%%%%%%
sim_v=sim(net, annv);

  
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%Saving informations
%%%%%%%%%%%%%
%%%%%%%%%%%%% 

for k1=1:n_1*n_2
 av{k1}=find(sim_v(k1,:)==1);
end


for k1=1:n_1*n_2
inpv{k1}=annv(:,av{k1});
end
 


m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
    m1=m1+1;
    
    Tabv_1{h1,h2}=av{m1};
    
end
end



%%%%%%%%%%%%%
%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%
%%%%%%%%%%%%%
figure(1)
    plotsomnd(net,annv) %MATLAB som built-in SOM plots; shows distance between each neuran' Neighbours
figure(2)
    plotsomhits(net,annv) %MATLAB som built-in SOM plots; shows density of each neurans
figure(3)
    plotsomplanes_sahar(net) %MATLAB som built-in SOM plots; shows density of each neurans

    saveas(figure(2),'~/Desktop/v_with_m101.jpeg','jpeg')
     saveas(figure(3),'~/Desktop/planes_weight.jpeg','jpeg')
   

   