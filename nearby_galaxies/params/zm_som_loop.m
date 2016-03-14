clear
clc

dir= '~/Dropbox/for_SOM/Hossein_s_Data/';
  
%>>>>  load your input dara in the Excel format
%>>>> new name 'cat':  M x N  (M=regions, N= parameters; this is the format of original file)
%>>>> by this command 'Null' is converted to NaN which is readable by Matlab 

%[cat] = xlsread('m31.csv');  
%>>>>---------------------------------------------------------------------------------------------------
%>>>> if all components are numbers then you can us 'txt' format then run this to load the data and inactive the above commend by % 

load ~/Downloads/sent-Oct14/inp_scat.txt
cat=inp_scat;
%>>>>---------------------------------------------------------------------------------------------------
%>>>>  Inverse  'cat' to N x M  for normalization and other processes

 cat=cat';        
%>>>>---------------------------------------------------------------------------------------------------
%>>>> fix NaN by takining avarage and assign a number to NaN and flag the NaN
%>>>> Network can not accept NaN or NULL

cat_fix=fixunknowns(cat);
%>>>>---------------------------------------------------------------------------------------------------
%>>>> normalization data, select only one  
%>>>> mapminmax: mormalization between -1 and 1.
%>>>> mapstd: Gaussian normalization to sigma=1 and mean=0
%>>>> cat_fix_norm= (cat_fix)  DO Nothing 

cat_fix_norm= (cat_fix);
%cat_fix_norm= mapminmax(cat_fix);
%cat_fix_norm= mapstd(cat_fix);

%>>>>---------------------------------------------------------------------------------------------------
%>>>> a  name change only to introduce t network
annt=cat_fix_norm;
sz=size(annt); %finding size of original data
nums=sz(2); % #of galaxies (regions)
%>>>>---------------------------------------------------------------------------------------------------
%>>>> Map size (n_1 x n_2) of the Network parametrs 
for n_1=5:5:50
   for n_2=n_1:5:50
     n1st=int2str(n_1);
     n2st=int2str(n_2);
%n_1=5
%n_2=5
%>>>>---------------------------------------------------------------------------------------------------
%>>>>  Parameters of Neighbours (n_nei) and number of training steps (n_cen) 
%>>>> the smaler n_cen the more separated groups (more covering space) 
%>>>> each neuran can be connected with (n_nei) nth Neighbours

n_cen=1;
n_nei=3;
%>>>>---------------------------------------------------------------------------------------------------
%>>>> MATLAB NETWORK  ; should not change
 net = newsom(annt,[n_2,n_1],'hextop','linkdist',n_cen,n_nei);
 net.trainParam.epochs = 200;
 net.trainParam.showWindow =false;
 net.trainParam.showCommandLine=true;
 net.trainParam.show=10;
 net = train(net,annt);
  sim_t=sim(net,annt);

%>>>>---------------------------------------------------------------------------------------------------
%>>>> Generating plots and tables 
%>>>> should not change


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
    
    TAB_1{h1,h2}=at{m1}';
    dummy=size(TAB_1{h1,h2});
    pers_result{h1,h2}=dummy(1)*100/nums; %%Shows percentage of regions in each neuran 
    
end
end


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


m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
    m1=m1+1;
    
    CAT_1{h1,h2}=cat(:,at{m1});
    
end
end

colorim = strcat(dir,'dist',n1st,'by',n2st,'.pdf');
resim = strcat(dir,'mat',n1st,'by',n2st,'.pdf');
pers= strcat(dir,'pers',n1st,'by',n2st,'.csv');
pos = strcat(dir,'pos',n1st,'by',n2st,'.csv');

figure(1)
 plotsomnd(net,annt);
saveas(figure(1),colorim,'pdf')
figure(2)
plotsomhits(net,annt);
saveas(figure(2),resim,'pdf')


table = cell2table(pers_result);
writetable(table,pers);

csvwrite(pos,Mtx_TAB_1);



 end
end


%>>>>---------------------------------------------------------------------------------------------------



