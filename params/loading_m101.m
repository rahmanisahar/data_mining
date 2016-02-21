

load ~/Desktop/project/data_mining/SOM/derived/per_arcsec_sq/same_as_M101/net10by10.mat %Load network

cat_new = csvread('~/Desktop/project/data_mining/m101/ascii_table/m101_total_per_arcsecsq.csv',1,2); %load data
cat_new=cat_new./10;
cat_new=cat_new';
cat_fix=fixunknowns(cat_new);
cat_fix_norm= (cat_fix);
annt=cat_fix_norm; %changing namme to introduce to network
sz=size(annt); %finding size of original data
 n_1=10;
 n_2=10;


%%%%%%%%%%%%%
%%%%%%%%%%%%%
%Giving data to our network
%%%%%%%%%%%%%
%%%%%%%%%%%%%
sim_n=sim(net, annt);


%%%%%%%%%%%%%
%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%
%%%%%%%%%%%%%
figure(1)
    plotsomnd(net,annt) %MATLAB som built-in SOM plots; shows distance between each neuran' Neighbours
figure(2)
    plotsomhits(net,annt) %MATLAB som built-in SOM plots; shows density of each neurans
figure(3)
   plotsompos(net,annt)
   
  
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%Saving informations
%%%%%%%%%%%%%
%%%%%%%%%%%%% 
for k1=1:n_1*n_2
 at{k1}=find(sim_n(k1,:)==1); %create 1X(n1*n2) cell and show regions which are going to same neuran 
end 
for k1=1:n_1*n_2
inpt{k1}=annt(:,at{k1}); % initial data in for each of the regions in each "at" cell
end
nums=sz(2); 
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
    
    CAT_1{h1,h2}=cat_new(:,at{m1}); % initial data in for each of the regions in each neuran
    
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


   