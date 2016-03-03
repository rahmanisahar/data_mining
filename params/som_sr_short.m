function som_sr_short(infile,col1,col2,n_1,n_2,dir)

clc
%This function creats a som for available data in m31 (54 raws for 10
%regions) and plots som plots plus F(FUV)/F(NUV) vs RHI and PAH line in
%lambda vs all the other parameters in M31. 
% Map size (n_1 x n_2) of the Network parametrs 
%dir: results' file directory 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Load all necessarily data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%>>>> Loading data in "txt" format
%>>>> Make sure there is no "NAN" values in this file; replace them with proper # or remove them
%load ~/Desktop/project/data_mining/m31/ascii_tables/m31_table_without_UBVRIJHKs_nohd_noNAN.txt
%cat = m31_table_without_UBVRIJHKs_nohd_noNAN; % new name 'cat':  M x N  (M=regions, N= parameters; this is the format of original file)


%>>>> Loading data in "CSV" format; When first raw is header and first
%three column are region name, RA and DEC, respectively.
cat = csvread(infile,1,col1,[1,col1,10,col2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Making file ready for Network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat=cat'; %Inverse  'cat' to N x M

%>>>> If there is any "NAN", fix it by takining avarage and assign a number to NaN and flag the NaN
%>>>> Network can not accept NaN or NULL
 cat_fix=fixunknowns(cat);

%>>>> normalization data, select only one  
%>>>> mapminmax: mormalization between -1 and 1.
%>>>> mapstd: Gaussian normalization to sigma=1 and mean=0
%>>>> cat_fix_norm= (cat_fix)  DO Nothing 
 %cat_fix_norm= (cat_fix);
 %cat_fix_norm= mapminmax(cat_fix);
y_min = -1;
y_max = 1;
sz = size(cat_fix);
catt_min = min(cat_fix')' * ones(1,sz(2));
catt_max = max(cat_fix')' * ones(1,sz(2));
cat_fix_norm = (y_max - y_min) * (cat_fix - catt_min) ./ (catt_max - catt_min) + y_min;
 
 %cat_fix_norm= mapstd(cat_fix);

annt=cat_fix_norm; %changing namme to introduce to network
sz=size(annt); %finding size of original data
nums=sz(2); % #of regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Create a SOM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%n_cen: Ordering phase steps
%n_nei: Tuning phase neighborhood distance, default = 1.
n_cen=500;
%>>>> MATLAB NETWORK
 net = newsom(annt,[n_2,n_1],'hextop','linkdist',0.9,n_cen,0.02,1);
 net.trainParam.epochs = 2000;
 net.trainParam.showWindow =false;
 net.trainParam.showCommandLine=true;
 net.trainParam.show=100;
 net = train(net,annt);
 sim_t=sim(net,annt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Generating plots and tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %for nameing
  n1st=int2str(n_1);
  n2st=int2str(n_2);
  col1st=int2str(col1);
  col2st=int2str(col2);
%  neist=int2str(n_nei);
 % n_censt=int2str(n_cen);
for k1=1:n_1*n_2
 at{k1}=find(sim_t(k1,:)==1); %create 1X(n1*n2) cell and show regions which are going to same neuran 
end 
for k1=1:n_1*n_2
inpt{k1}=annt(:,at{k1}); % initial data in for each of the regions in each "at" cell
end

m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
    m1=m1+1;
    TAB_1{h1,h2}=at{m1}'; % create n1 X n2 cell; shows which regions goes to which Neuran 
end
end

m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
    m1=m1+1;
    
    CAT_1{h1,h2}=cat(:,at{m1}); % initial data in for each of the regions in each neuran
    
end
end

%>>>>Following creats matrix version of TAB_1 an to be able to save them as
%csv files
  max_n = 0;
  for j=1:n_1*n_2
    size_tab_1 = size(TAB_1{j}); 
    if size_tab_1(1) > max_n 
       max_n = size_tab_1(1);
    end
  end 
 
   Mtx_TAB_1 = zeros(n_1*n_2,max_n);
  for j=1:n_1*n_2
    size_temp = size(TAB_1{j});
    Mtx_TAB_1(1:size_temp(1),j) = TAB_1{j};
  end
  Mtx_TAB_1=Mtx_TAB_1';

  
  names=char('names','RADEg','DECDEG','PAH6.2flx','PAH7.7flx','PAH8.3flx','PAH8.6flx','PAH11.3flx','PAH12.0flx', ...
           'PAH12.7flx','PAH17.0flx','halpha','OIII continuum','SII continumm','IRAC3', ...
           'PACS100','SPIRE250','SPIRE350','SPIRE500','Dust luminosity','Dust Mass', ...
           'SFR','Stellar Mass', 'TIR', 'Total gas mass', 'RHI','X12plogOH');
diff=col1-1;
for K = col1:1:col2
    
    new_names(K-diff,:)=names(K+1,:);
   
end

net.inputs{1}.userdata = new_names;

%>>>>end

figure(1)
    plotsomnd(net,annt) %MATLAB som built-in SOM plots; shows distance between each neuran' Neighbours
 
figure(2)
    plotsomhits(net,annt) %MATLAB som built-in SOM plots; shows density of each neurans
%figure(3)
 %  plotsomplanes(net,annt)
figure(4)
    plotsomplanes_sahar(net) %MATLAB som built-in SOM plots with small changes; shows weights for each input


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Saving Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fig1 = strcat(dir,'dist',n1st,'by',n2st,'_data_between_cols',col1st,'and',col2st,'.jpeg');
  fig2 = strcat(dir,'hits',n1st,'by',n2st,'_data_between_cols',col1st,'and',col2st,'.jpeg');
  %fig3 = strcat(dir,'weigth',n1st,'by',n2st,'_data_between_cols',col1st,'and',col2st,'.jpeg');
  fig4 = strcat(dir,'weigth_planes',n1st,'by',n2st,'_data_between_cols',col1st,'and',col2st,'.jpeg');
  %pers= strcat(dir,'pers',n1st,'by',n2st,'_data_between_cols',col1st,'and',col2st,'.csv');
  pos = strcat(dir,'pos',n1st,'by',n2st,'_data_between_cols',col1st,'and',col2st,'.csv');
  nett= strcat(dir,'net',n1st,'by',n2st,'_data_between_cols',col1st,'and',col2st);
save(nett, 'net')
%table = cell2table(pers_result);
%writetable(table,pers);
csvwrite(pos,Mtx_TAB_1); 
saveas(figure(1),fig1,'jpeg')
saveas(figure(2),fig2,'jpeg')
%saveas(figure(3),fig3,'jpeg')
saveas(figure(4),fig4,'jpeg')
close all
end