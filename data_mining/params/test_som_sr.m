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
  
  
for h1=n_1:-1:1
    for   h2=1:1:n_2
        
    end
end
%>>>>end

figure(1)
    plotsomnd(net,annt) %MATLAB som built-in SOM plots; shows distance between each neuran' Neighbours
 
figure(2)
    plotsomhits(net,annt) %MATLAB som built-in SOM plots; shows density of each neurans

figure(3)  
    for h1=n_1:-1:1
        m1=0;
            for h2=1:1:n_2
                m1=m1+1;
                check_s=CAT_1{h1,h2};
                size_ch=size(check_s);
                if (size_ch(2) > 0) 
                    params=CAT_1{h1,h2};
                    s=scatter(params(24,:)./params(25,:),params(51,:));
                    s.MarkerFaceColor = [m1/10 h1/10 h2/10];
                    hold on
                    xlabel('F(FUV)/F(NUV)')
                    ylabel('RHI')
                end
            end
    end
%>>>> test plots
count=2;
for i=1:17
    i_name=int2str(i);
    m1=0;
    count=count+2;
    figure(count)
    for ind=18:33
        m1=m1+1;
        subplot(4,4,m1)
        for h1=n_1:-1:1
            for h2=1:1:n_2
                check_s=CAT_1{h1,h2};
                size_ch=size(check_s);
                if (size_ch(2) > 0) 
                    params=CAT_1{h1,h2};
                    s=scatter(params(ind,:),params(i,:));
                    s.MarkerFaceColor = [m1/20 h1/10 h2/10];
                    hold on
                end
            end
        end
    end
    name1 = strcat(dir,i_name,'_vs_raws_21_to_36_for_',n1st,'by',n2st,'.pdf');
    saveas(figure(count),name1,'pdf')    
    m1=0;
    figure(count+1)
    for ind=34:49
        m1=m1+1;
        subplot(4,4,m1)
        for h1=n_1:-1:1
            for h2=1:1:n_2
                check_s=CAT_1{h1,h2};
                size_ch=size(check_s);
                if (size_ch(2) > 0) 
                    params=CAT_1{h1,h2};
                    s=scatter(params(ind,:),params(i,:));
                    s.MarkerFaceColor = [m1/20 h1/10 h2/10];
                    hold on
                end
            end
        end
    end
    name2 = strcat(dir,i_name,'_vs_raws_37_to_53_for_',n1st,'by',n2st,'.pdf');
    saveas(figure(count+1),name2,'pdf')    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Saving Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fig1 = strcat(dir,'dist',n1st,'by',n2st,'.pdf');
  fig2 = strcat(dir,'hits',n1st,'by',n2st,'.pdf');
  fig3 = strcat(dir,'plot_radiation_hardness_index_vs_ratio_of_galex',n1st,'by',n2st,'.pdf');
  
  
saveas(figure(1),fig1,'pdf')
saveas(figure(2),fig2,'pdf')
saveas(figure(3),fig3,'pdf')
