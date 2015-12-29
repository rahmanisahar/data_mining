clear
clc

dir= '~/Desktop/project/data_minning/new_paper_2012/'; %results' file directory 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Load all necessarily data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%>>>> Loading data in "txt" format
%>>>> Make sure there is no "NAN" values in this file; replace them with proper # or remove them
load m31.txt
cat = m31; % new name 'cat':  M x N  (M=regions, N= parameters; this is the format of original file)
%>>>> If there is any "NAN", fix it by takining avarage and assign a number to NaN and flag the NaN
%>>>> Network can not accept NaN or NULL
 cat_fix=fixunknowns(cat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Making file ready for Network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat=cat'; %Inverse  'cat' to N x M

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
%>>>> Map size (n_1 x n_2) of the Network parametrs 
 n_1=1;
 n_2=12;

n_cen=100; % number of training steps; the smaler n_cen the more separated groups (more covering space)
n_nei=2;  % #of Neighbours; each neuran can be connected with (n_nei) nth Neighbours

%>>>> MATLAB NETWORK
 net = newsom(annt,[n_2,n_1],'hextop','linkdist',n_cen,n_nei);
 net.trainParam.epochs = 200;
 net.trainParam.showWindow =false;
 net.trainParam.showCommandLine=true;
 net.trainParam.show=10;
 net = train(net,annt);
 sim_t=sim(net,annt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Generating plots and tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  Mtx_TAB_1 = zeros(max,n_1*n_2);
  for j=1:n_1*n_2
    size_temp = size(TAB_1{j});
    Mtx_TAB_1(1:size_temp(1),j) = TAB_1{j};
  end
%>>>>end
figure(1)
    plotsomnd(net,annt) %MATLAB som built-in SOM plots; shows distance between each neuran' Neighbours
 
figure(2)
    plotsomhits(net,annt) %MATLAB som built-in SOM plots; shows density of each neurans
 
figure(3) %test plots
    
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Saving Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for nameing
  n1st=int2str(n_1);
  n2st=int2str(n_2);

