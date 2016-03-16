clear
clc

dir= '~/Desktop/';
%>>>>  load your input dara in the Excel format
%>>>> new name 'cat':  M x N  (M=regions, N= parameters; this is the format of original file)
%>>>> by this command 'Null' is converted to NaN which is readable by Matlab 

%[cat] = xlsread('/Users/Andromeda/Dropbox/for_SOM/original_data_m31_same_as_m101 Andromeda/m31_secondry_per_pc_nohd.csv');  
%>>>>---------------------------------------------------------------------------------------------------
%>>>> if all components are numbers then you can us 'txt' format then run this to load the data and inactive the above commend by % 

%load m31.txt
load ~/Desktop/project/data_mining/m31/acsii_tables/total_m31_table_nohd_noNAN.txt
cat=total_m31_table_nohd_noNAN;
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
%>>>> a  name change only to introduce to network
annt=cat_fix_norm;
sz=size(annt); %finding size of original data
nums=sz(2); % #of galaxies (regions)
%>>>>---------------------------------------------------------------------------------------------------
%>>>> Map size (n_1 x n_2) of the Network parametrs 

n_1=2;
n_2=2;
n1st=int2str(n_1);
n2st=int2str(n_2);
%>>>>---------------------------------------------------------------------------------------------------
%>>>>  Parameters of Neighbours (n_nei) and number of training steps (n_cen) 
%>>>> the smaler n_cen the more separated groups (more covering space) 
%>>>> each neuran can be connected with (n_nei) nth Neighbours

n_cen=100;    
n_nei=1;  
%>>>>---------------------------------------------------------------------------------------------------
%>>>> MATLAB NETWORK  ; should not change
 net = newsom(annt,[n_2,n_1],'hextop','linkdist',n_cen,n_nei);
 net.trainParam.epochs = 500;
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

colorim=strcat(dir,'dist',n1st,'by',n2st,'.pdf');
resim=strcat(dir,'mat',n1st,'by',n2st,'.pdf');
pers=strcat(dir,'pers',n1st,'by',n2st,'.csv');
pos=strcat(dir,'pos',n1st,'by',n2st,'.csv');

figure(1)
 plotsomnd(net,annt);
saveas(figure(1),'~/Desktop/dist_only_observed_values.pdf','pdf')

%figure(2)
 %plotsomnc(net,annt);
%saveas(figure(2),'~/Desktop/project/data_minning/new_paper_2012/connections_only_observed_values.pdf','pdf')

figure(3)
plotsomhits(net,annt);
saveas(figure(3),'~/Desktop/hits_only_observed_values.pdf','pdf')

%figure(4)
%plotsomplanes(net,annt);
%saveas(figure(4),'~/Desktop/project/data_minning/new_paper_2012/weight_planes_only_observed_values.pdf','pdf')

%figure(5)
%plotsompos(net,annt);
%saveas(figure(5),'~/Desktop/project/data_minning/new_paper_2012/weight_positions_only_observed_values.pdf','pdf')


%figure(6)
%plotsomtop(net,annt);
%saveas(figure(6),'~/Desktop/project/data_minning/new_paper_2012/topology_only_observed_values.pdf','pdf')



load ~/Desktop/project/data_minning/new_paper_2012/ks_n.txt
 model=ks_n;

load ~/Desktop/project/data_minning/new_paper_2012/b21s_n.txt
sed=b21s_n;


load ~/Desktop/project/data_minning/new_paper_2012/xx.txt
x_ax=xx;


m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
    m1=m1+1;
    
    sep_sed{h1,h2}=sed(:,TAB_1{h1,h2});
    
end
end

m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
    m1=m1+1;
    ch=sep_sed{h1,h2};
    size_ch=size(ch);
    if (size_ch(2) > 1) 
         trans=transpose(ch); 
         catm{h1,h2}=mean(trans);
    else
        trans=ch; 
        catm{h1,h2}=trans;
    end
    
    
end
end



figure(10)
m1=0;
type=zeros(1,12);
for h1=1:1:n_1
 for   h2=1:1:n_2
    m1=m1+1;
    chi_s=zeros(1,12);
    a = catm{h1,h2};
    for n1=1:1:12
       chi = 0;
       for i = 1:1:630
           chi = chi+(a(i)-model(i,n1))^2/model(i,n1);
       end
       chi_s(n1)=chi;
    end   
    [dum,num_model]=min(chi_s(:));
    %index=find(min(chi_s));
    subplot(n_1,n_2,m1)
    plot(xx,catm{h1,h2},'K');%,xlim([1 630])
    type(m1)=num_model;
    hold on
    plot(xx,model(1:630,num_model),'r');%,xlim([1 630])
    strmax = ['type = ',num2str(num_model)];
    text(3500,2,strmax,'HorizontalAlignment','right');
    
    
end
end

saveas(figure(10),'~/Desktop/project/data_minning/new_paper_2012/SED_only_observed_values.png','png')



figure(11);
m1=0;
for h1=1:1:n_1
 for   h2=1:1:n_2
    m1=m1+1;
    
subplot(n_1,n_2,m1)
    plot(xx,model(1:630,m1),'r')%,xlim([1 630])
    
end
end


load ~/Desktop/project/data_minning/new_paper_2012/bay_g21.txt
params=bay_g21;
%SSFR = params(:,33)-params(:,29);

for h1=n_1:-1:1
 for   h2=1:1:n_2
    
    stellar_mass{h1,h2}=params(TAB_1{h1,h2},31);
    D4000{h1,h2}=params(TAB_1{h1,h2},37);
    ssfr{h1,h2}=params(TAB_1{h1,h2},33)-params(TAB_1{h1,h2},29)-9;
    
    
end
end

figure(12)

subplot(3,1,1)
m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
     
     m1=m1+1;
     if (type(m1) == 1) 
       scatter(ssfr{h1,h2},D4000{h1,h2},'s','m')
     elseif (type(m1)==2) || (type(m1)==3)
        scatter(ssfr{h1,h2},D4000{h1,h2},'s', 'r')
     elseif (type(m1)==4) || (type(m1)==5) || (type(m1)==6)
         scatter(ssfr{h1,h2},D4000{h1,h2},'s', 'y')
     elseif (type(m1)==7) || (type(m1)==8) || (type(m1)==9)
         scatter(ssfr{h1,h2},D4000{h1,h2},'s', 'b')
     else
         scatter(ssfr{h1,h2},D4000{h1,h2},'s', 'g')
     end
   hold on
   xlabel('log \phi [Gyr^{-1}]')
   ylabel('log ^{t}D4000 [Gyr]')
 end
end
subplot(3,1,2)
m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
     
     m1=m1+1;
     if (type(m1) == 1) 
       scatter(stellar_mass{h1,h2},ssfr{h1,h2},'s','m')
     elseif (type(m1)==2) || (type(m1)==3)
        scatter(stellar_mass{h1,h2},ssfr{h1,h2},'s', 'r')
     elseif (type(m1)==4) || (type(m1)==5) || (type(m1)==6)
         scatter(stellar_mass{h1,h2},ssfr{h1,h2},'s', 'y')
     elseif (type(m1)==7) || (type(m1)==8) || (type(m1)==9)
         scatter(stellar_mass{h1,h2},ssfr{h1,h2},'s', 'b')
     else
         scatter(stellar_mass{h1,h2},ssfr{h1,h2},'s', 'g')
     end
   hold on 
   xlabel('log M_{star} [M_{Sun}]')
   ylabel('log \phi [Gyr^{-1}]')
 end
end
subplot(3,1,3)
m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
     
     m1=m1+1;
     if (type(m1) == 1) 
       scatter(stellar_mass{h1,h2},D4000{h1,h2},'s','m')
     elseif (type(m1)==2) || (type(m1)==3)
        scatter(stellar_mass{h1,h2},D4000{h1,h2},'s', 'r')
     elseif (type(m1)==4) || (type(m1)==5) || (type(m1)==6)
         scatter(stellar_mass{h1,h2},D4000{h1,h2},'s', 'y')
     elseif (type(m1)==7) || (type(m1)==8) || (type(m1)==9)
         scatter(stellar_mass{h1,h2},D4000{h1,h2},'s', 'b')
     elseif (type(m1)==10) || (type(m1)==11) || (type(m1)==12)
         scatter(stellar_mass{h1,h2},D4000{h1,h2},'s', 'g')
     end
   hold on  
   xlabel('log M_{star} [M_{Sun}]')
   ylabel('log ^{t}D4000 [Gyr]')
 end
end
   
 saveas(figure(12),'~/Desktop/project/data_minning/new_paper_2012/rep_fig_16_only_observed_values.png','png')  
   
   
   %>>>>---------------------------------------------------------------------------------------------------



