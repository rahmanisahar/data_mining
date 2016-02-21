function som_sr_short(filename,n_1,n_2,dir)

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
cat = csvread(filename,1,4);
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
 cat_fix_norm= (cat_fix);
 %cat_fix_norm= mapminmax(cat_fix);
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
n_cen=1000;
%>>>> MATLAB NETWORK
 net = newsom(annt,[n_2,n_1],'hextop','linkdist',0.9,n_cen,0.02,3);
 net.trainParam.epochs = n_cen*3;
 net.trainParam.showWindow =false;
 net.trainParam.showCommandLine=true;
 net.trainParam.show=1000;
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
    dummy=size(TAB_1{h1,h2});
    pers_result{h1,h2}=dummy(1)*100/nums; %%Shows percentage of regions in each neuran 
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
  max = 0;
  for j=1:n_1*n_2
    size_tab_1 = size(TAB_1{j}); 
    if size_tab_1(1) > max 
       max = size_tab_1(1);
    end
  end 
 
   Mtx_TAB_1 = zeros(n_1*n_2,max);
  for j=1:n_1*n_2
    size_temp = size(TAB_1{j});
    Mtx_TAB_1(1:size_temp(1),j) = TAB_1{j};
  end
  Mtx_TAB_1=Mtx_TAB_1';

  
 
%>>>>end

figure(1)
    plotsomnd(net,annt) %MATLAB som built-in SOM plots; shows distance between each neuran' Neighbours
 
figure(2)
    plotsomhits(net,annt) %MATLAB som built-in SOM plots; shows density of each neurans
figure(3)
   plotsompos(net,annt)
figure(4)
    plotsomnd_s(net,annt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Saving Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fig1 = strcat(dir,'dist',n1st,'by',n2st,'.jpeg');
  fig2 = strcat(dir,'hits',n1st,'by',n2st,'.jpeg');
  fig3 = strcat(dir,'weigth',n1st,'by',n2st,'.jpeg');
  pers= strcat(dir,'pers',n1st,'by',n2st,'.csv');
  pos = strcat(dir,'pos',n1st,'by',n2st,'.csv');
  nett= strcat(dir,'net',n1st,'by',n2st);
save(nett, 'net')
table = cell2table(pers_result);
writetable(table,pers);
csvwrite(pos,Mtx_TAB_1); 
saveas(figure(1),fig1,'jpeg')
saveas(figure(2),fig2,'jpeg')
saveas(figure(3),fig3,'jpeg')
%close all
end