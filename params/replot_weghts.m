close all, clear, clc
net10by10=load('~/Desktop/project/data_mining/SOM/derived/per_arcsec_sq/diff_dimention_10by10/nei1/nets/net10by10_data_between_cols3and26.mat');
net=net10by10.net;
catv = csvread('~/Desktop/project/data_mining/m31/ascii_tables/derived_ones_with_mean_per_arcsec.csv',1,3,[1,3,10,26]);
catv=catv';
catv_fix=fixunknowns(catv);
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

names=char('PAH6.2flx','PAH7.7flx','PAH8.3flx','PAH8.6flx','PAH11.3flx','PAH12.0flx', ...
           'PAH12.7flx','PAH17.0flx','H\alpha','OIII continuum','SII continumm','IRAC 5.7 \mum', ...
           'PACS 100 \mum','SPIRE 250 \mum','SPIRE 350 \mum','SPIRE 500 \mum','Dust luminosity','Dust Mass', ...
           'SFR', 'Stellar Mass', 'TIR', 'Total gas mass', 'RHI','Metallicity');
net.inputs{1}.userdata = names; 
figure(1)
plotsomplanes_sahar(net)

saveas(figure(1),'~/Desktop/project/data_mining/SOM/derived/per_arcsec_sq/subsets/all.jpeg','jpeg')
