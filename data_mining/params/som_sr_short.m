function som_sr_short(n_1,n_2,n_cen,n_nei,dir)
clc
%This function creats a som for available data in m31 (54 raws for 10
%regions) and plots som plots plus F(FUV)/F(NUV) vs RHI and PAH lin in
%lambda vs all the other parameters in M31. 
% Map size (n_1 x n_2) of the Network parametrs 
%n_cen number of training steps; the smaler n_cen the more separated groups
%(more covering space) (suggestion between 1 to 200)
%n_nei of Neighbours; each neuran can be connected with (n_nei) nth
%Neighbours (suggestion between 1 to 5)
%dir: results' file directory 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Load all necessarily data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%>>>> Loading data in "txt" format
%>>>> Make sure there is no "NAN" values in this file; replace them with proper # or remove them
load ~/Desktop/project/data_mining/m31/acsii_tables/m31_table_without_UBVRIJHKs_nohd_noNAN.txt
cat = m31_table_without_UBVRIJHKs_nohd_noNAN; % new name 'cat':  M x N  (M=regions, N= parameters; this is the format of original file)

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
   %for nameing
  n1st=int2str(n_1);
  n2st=int2str(n_2);
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
  Mtx_TAB_1=Mtx_TAB_1'

  
  
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
   plotsompos(net,annt)
% figure(3)  
%     for h1=n_1:-1:1
%         m1=0;
%             for h2=1:1:n_2
%                 m1=m1+1;
%                 check_s=CAT_1{h1,h2};
%                 size_ch=size(check_s);
%                 if (size_ch(2) > 0) 
%                     params=CAT_1{h1,h2};
%                     s=scatter(params(42,:),params(24,:)./params(25,:));
%                     s.MarkerFaceColor = [0 h1/20 h2/20];
%                     hold on
%                     ylabel('F(FUV)/F(NUV)')
%                     xlabel('RHI')
%                 end
%             end
%     end
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
                    s.MarkerFaceColor = [0 h1/20 h2/20];
                    hold on
                end
            end
        end
    end
    name1 = strcat(dir,i_name,'vs_raws_18_to_33_for_',n1st,'by',n2st,'.pdf');
    saveas(figure(count),name1,'pdf')    
    m1=0;
    figure(count+1)
    for ind=34:43
        m1=m1+1;
        subplot(2,5,m1)
        for h1=n_1:-1:1
            for h2=1:1:n_2
                check_s=CAT_1{h1,h2};
                size_ch=size(check_s);
                if (size_ch(2) > 0) 
                    params=CAT_1{h1,h2};
                    s=scatter(params(ind,:),params(i,:));
                    s.MarkerFaceColor = [0 h1/20 h2/20];
                    hold on
                end
            end
        end
    end
    name2 = strcat(dir,i_name,'vs_raws_34_to_49_for_',n1st,'by',n2st,'.pdf');
    saveas(figure(count+1),name2,'pdf')    
end

for i=1:17
    for ind=18:43
    i_name=int2str(i);
    j_name=int2str(ind);
    corr_mat=zeros(n_1,n_2)+0.000000006;
        for h1=n_1:-1:1
            for h2=1:1:n_2
                check_s=CAT_1{h1,h2};
                size_ch=size(check_s);
                if (size_ch(2) > 2) 
                    params=CAT_1{h1,h2};
                    corr_mat(h1,h2)=corr(params(i,:)',params(ind,:)');
                end
            end
        end
        if (abs(corr_mat) > 0.5)
          name1 = strcat(dir,i_name,'_',j_name,n1st,'_','by',n2st,'.csv');   
          % table1 = cell2table(corr_mat);
         csvwrite(name1,corr_mat);
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Saving Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fig1 = strcat(dir,'dist',n1st,'by',n2st,'.pdf');
  fig2 = strcat(dir,'hits',n1st,'by',n2st,'.pdf');
  fig3 = strcat(dir,'weigth',n1st,'by',n2st,'.pdf');
  pers= strcat(dir,'pers',n1st,'by',n2st,'.csv');
  pos = strcat(dir,'pos',n1st,'by',n2st,'.csv');

table = cell2table(pers_result);
writetable(table,pers);
csvwrite(pos,Mtx_TAB_1); 
saveas(figure(1),fig1,'pdf')
saveas(figure(2),fig2,'pdf')
saveas(figure(3),fig3,'pdf')
close all
end