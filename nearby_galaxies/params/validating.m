
close all, clear, clc
%net10by10=load('~/Desktop/project/data_mining/nearby_galaxies/SOM/grey_cycle/2d/subsets/without_reg1/sfr/net_subset19.mat'); %Load network
net10by10=load('~/Desktop/net_same_as_m101.mat'); %Load network
net=net10by10.net;
%catv = csvread('~/Desktop/project/data_mining/nearby_galaxies/m31/ascii_tables/subsets/subset10.csv',1,3);%col1,[1,col1,10,col2]); %load data
catv = csvread('~/Desktop/project/data_mining/nearby_galaxies/m101/ascii_table/m101_total_with_mean_per_arcsecsq_same_as_derived_ones_2.csv',1,1);
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
annv=catv_fix_norm;%(:,1); %changing namme to introduce to network
annt=catv_fix_norm;%(:,2:10);
sz=size(annv); %finding size of original data
 n_1=10;
 n_2=10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Giving data to our network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_v=sim(net, annv);
sim_t=sim(net, annt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Saving informations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for k1=1:n_1*n_2
 at{k1}=find(sim_t(k1,:)==1);
 sz=size(at{k1});
    if (sz(2) == 1)
 	at_2(k1)=at{k1}(1);
end

for k1=1:n_1*n_2
 av{k1}=find(sim_v(k1,:)==1);

 end
end

for k1=1:n_1*n_2
inpv{k1}=annv(:,av{k1});
end
 for k1=1:n_1*n_2
inpt{k1}=annt(:,at{k1});
end


m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
    m1=m1+1;
    
    Tabv_1{h1,h2}=av{m1};
    
end
end

m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
    m1=m1+1;
    
    Tabt_1{h1,h2}=at{m1};
    
end
end

for h1=n_1:-1:1
 for   h2=1:1:n_2
    sz=size(Tabt_1{h1,h2});
    if (sz(2) == 1)
    
    Tabt(h1,h2)=Tabt_1{h1,h2}(1);
end
    
end
end


for h1=n_1:-1:1
 for   h2=1:1:n_2
    sz=size(Tabv_1{h1,h2});
    if (sz(2) == 1)
    
    Tabv(h1,h2)=Tabv_1{h1,h2}(1);

end
    
end
end

% [nn1, mm1]=find(Tabv==1);
% nonem_at=find(~cellfun(@isempty,at));





% at_m=sparse(Tabt);
%[d pred]=shortest_paths(at_m,[nn1,mm1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting and saving information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%csvwrite('~/Desktop/levels_without_stellar_mass_2.csv',(1-levels)*100); 
figure(2)
    plotsom_sahar(net,annv) %MATLAB som built-in SOM plots; shows density of each neurans

    %saveas(figure(2),'~/Desktop/project/data_mining/nearby_galaxies/SOM/grey_cycle/2d/subsets/without_reg1/where_is_reg_1_without_stellar_mass.png','png')






