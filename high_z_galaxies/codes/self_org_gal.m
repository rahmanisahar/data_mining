close all, clear, clc
dir='~/Desktop/';
load ~/Desktop/project/data_mining/high_z_galaxies/data/info_model.txt
params = info_model;
%>>>>  load your input dara in the Excel format
%>>>> new name 'cat':  M x N  (M=regions, N= parameters; this is the format of original file)
%>>>> by this command 'Null' is converted to NaN which is readable by Matlab 

%[cat] = xlsread('m31.csv');  
%>>>>---------------------------------------------------------------------------------------------------
%>>>> if all components are numbers then you can us 'txt' format then run this to load the data and inactive the above commend by % 

load ~/Desktop/project/data_mining/high_z_galaxies/data/kinney_n.txt
cat=kinney_n(:,1:end);

load ~/Desktop/project/data_mining/high_z_galaxies/data/model_n.txt
catv=model_n(:,1:end);

load ~/Desktop/project/data_mining/nearby_galaxies/H_paper_2012_som/xx.txt
x_ax=xx;
%>>>>---------------------------------------------------------------------------------------------------
%>>>>  Inverse  'cat' to N x M  for normalization and other processes
kinney=cat;
annm=catv';
catv=catv';
cat=cat';        

%>>>>---------------------------------------------------------------------------------------------------
%>>>> fix NaN by takining avarage and assign a number to NaN and flag the NaN
%>>>> Network can not accept NaN or NULL

cat_fix=(cat);
catv_fix=(catv);

%>>>>---------------------------------------------------------------------------------------------------
%>>>> normalization data, select only one  
%>>>> mapminmax: mormalization between -1 and 1.
%>>>> mapstd: Gaussian normalization to sigma=1 and mean=0
%>>>> cat_fix_norm= (cat_fix)  DO Nothing 

cat_fix_norm= (cat_fix);
%cat_fix_norm= mapminmax(cat_fix);
%cat_fix_norm= mapstd(cat_fix);

catv_fix_norm= (catv_fix);
%catv_fix_norm= mapminmax(catv_fix);
%catv_fix_norm= mapstd(catv_fix);

%>>>>---------------------------------------------------------------------------------------------------
%>>>> a  name change only to introduce t network
annt=cat_fix_norm';
annv=catv_fix_norm';

%>>>>---------------------------------------------------------------------------------------------------
%>>>> Map size (n_1 x n_2) of the Network parametrs 


%>>>>---------------------------------------------------------------------------------------------------
%>>>>  Parameters of Neighbours (n_nei) and number of training steps (n_cen) 
%>>>> the smaler n_cen the more separated groups (more covering space) 
%>>>> each neuran can be connected with (n_nei) nth Neighbours
n_1=12;
n_2=12;
 n1st=int2str(n_1);
 n2st=int2str(n_2);


ost=1000;
tnd=1;

orl=.9;
tlr=0.02;


n_epoch =ost*2;
%>>>>---------------------------------------------------------------------------------------------------
%>>>> MATLAB NETWORK  ; should not change
%net = newsom(PR,[D1,D2,...],TFCN,DFCN,OLR,OSTEPS,TLR,TND)
 net = selforgmap([n_2,n_1],ost,3,'hextop','linkdist');
 %net = newsom(annt,[n_2,n_1]);
 net.trainParam.epochs = n_epoch;
 net.trainParam.showWindow =false;
 net.trainParam.showCommandLine=true;
net.trainParam.show=10;
 net = train(net,annt);
  sim_t=sim(net,annt);
  sim_v=sim(net,annv);
 
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
    Tab_1{h1,h2}=at{m1};
 
end
end

%elseif
colorim = strcat(dir,'dist_',n1st,'_by_self_org_res_',n2st,'.png');
hitim_t= strcat(dir,'hit_t_',n1st,'_by_self_org_res_',n2st,'.png');
 figure(1)
 plotsomnd_shar(net,annt)
figure(2)
plotsomhits_sahar(net,annt)
 %figure(20)
 %plotsom_sahar(net,annt)
saveas(figure(1),colorim,'png')
saveas(figure(2),hitim_t,'png')
%-------------------

for k1=1:n_1*n_2
 av{k1}=find(sim_v(k1,:)==1);
end


for k1=1:n_1*n_2
inpv{k1}=annv(:,av{k1});
end
 


m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
    m1=m1+1;
    
    
    Tabv_1{h1,h2}=av{m1};

    
end
end
hitim_v= strcat(dir,'hit_v_',n1st,'_by_self_org_res_',n2st,'.png');
%  figure(3)
%  plotsomnd(net,annv)
figure(3)
plotsomhits(net,annv)
saveas(figure(3),hitim_v,'png')

if n_1==1

    m1=0;
    for h1=1:1:n_1
        for   h2=1:1:n_2
            m1=m1+1;
    
    
            MEDv_mass(h1,h2)=median(info_model(Tabv_1{h1,h2},31));
            MEDv_ssfr(h1,h2)=median(info_model(Tabv_1{h1,h2},33)-info_model(Tabv_1{h1,h2},31)-9);
            MEDv_age(h1,h2)=median(info_model(Tabv_1{h1,h2},37));
            MEDv_frac(h1,h2)=median(info_model(Tabv_1{h1,h2},49));

    
        end
    end

    x= 1:n_2;
    figure(5)
    subplot(2,2,1)
    plot(x,MEDv_mass, '-O')
    set(gca,'xtick',1:n_2);
    title('Stellar Mass')
    ylabel('log M_{star} [M_{Sun}]')
    subplot(2,2,2)
    plot(x,MEDv_age, '-O')
    set(gca,'xtick',1:n_2);
    title('Age')
    ylabel('log ^{t}D4000 [Gyr]')
    subplot(2,2,3)
    plot(x,MEDv_ssfr, '-O')
    set(gca,'xtick',1:n_2);
    title('SSFR')
    ylabel('log \phi [Gyr^{-1}]')
    subplot(2,2,4)
    plot(x,MEDv_frac, '-O')
    set(gca,'xtick',1:n_2);
    title('Fraction')
prop= strcat(dir,'prop',n1st,'_by_self_org_res_',n2st,'.png');
saveas(figure(5),prop,'png')
end
 

%>>>>---------------------------------------------------------------------------------------------------

%%Sahar's plots adds on

m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
    m1=m1+1;
    
    CATv_1{h1,h2}=catv(:,av{m1});
    
end
end

m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
    m1=m1+1;
    
    catm{h1,h2}=mean(annm(av{m1},:));
    
end
end



figure(10)
m1=0;
type=zeros(1,12);
% text(-10,10.2,'Wavelength($$\AA$$)','Interpreter', 'Latex')    
for h1=1:1:n_1
 for   h2=1:1:n_2
     check_s=CATv_1{h1,h2};
               size_ch=size(check_s);
               if (isnan(catm{h1,h2}) ~= 1) 
                if (size_ch(2) > 1)
    m1=m1+1;
    chi_s=zeros(1,12);
    a = catm{h1,h2};
    for n1=1:1:12
       chi = 0;
       for i = 1:1:630
           chi = chi+(a(i)-kinney(i,n1))^2/kinney(i,n1);
       end
       chi_s(n1)=chi;
    end   
    [dum,num_kinney]=min(chi_s(:));
   index=find(min(chi_s));
    %m1=m1+1;
    subplot(n_1,n_2,m1)
    plot(xx,catm{h1,h2},'K')%,xlim([1 630])
    type(m1)=num_kinney;
    %plot(xx,kinney(1:630,num_kinney),'r')%,xlim([1 630])
    %strmax = ['Type = ',num2str(num_kinney)];
    %text(3500,2,strmax,'HorizontalAlignment','right');
   ylabel('flux (arbitary unit)')
   xlabel('Wavelength($$\AA$$)','Interpreter', 'Latex' ) 
     end    
               end
 end
end
%text(-10,10.2,'Wavelength(\AA)')
SEDim = strcat(dir,'SED_total',n1st,'_by_self_org_res',n2st,'.png');
saveas(figure(10),SEDim,'png')
% figure(11)
% m1=0;
% 
% 
% for h1=1:1:n_1
%  for   h2=1:1:n_2
%   
%      m1=m1+1;
%     
% subplot(n_1,n_2,m1)
%    plot(kinney(1:630,m1),'r'),xlim([1 630])
%     
% end
% end


if n_1~=1
    m1=0;
    for h1=1:1:n_1
        for   h2=1:1:n_2
            m1=m1+1;
    
    
            MEDv_mass(h1,h2)=median(info_model(Tabv_1{h1,h2},31));
            MEDv_ssfr(h1,h2)=median(info_model(Tabv_1{h1,h2},33)-info_model(Tabv_1{h1,h2},31)-9);
            MEDv_age(h1,h2)=median(info_model(Tabv_1{h1,h2},37));
            MEDv_frac(h1,h2)=median(info_model(Tabv_1{h1,h2},49));

    
        end
    end

    for x=1:n_1
        for y=1:n_2
            
        figure(6)
        subplot(2,2,1)
        plot3(x,y,MEDv_mass(x,y),'-Ok')
        title('Stellar Mass')
        zlabel('log M_{star} [M_{Sun}]')
        hold on
        end
    end
     for x=1:n_1
        for y=1:n_2
    subplot(2,2,2)
    plot3(x,y,MEDv_age(x,y), '-Ok')
    title('Age')
    zlabel('log ^{t}D4000 [Gyr]')
    hold on
        end
     end
     for x=1:n_1
        for y=1:n_2
    subplot(2,2,3)
    plot3(x,y,MEDv_ssfr(x,y), '-Ok')
    title('SSFR')
    zlabel('log \phi [Gyr^{-1}]')
    hold on
        end
     end
     for x=1:n_1
        for y=1:n_2
    subplot(2,2,4)
    plot3(x,y,MEDv_frac(x,y), '-Ok')
    title('Fraction')
    hold on
        end
     end
    if (n_1 == n_2)
     x=1:n_1;
     y=1:n_2;
     figure(7)
     subplot(2,2,1)
        plot3(x,y,MEDv_mass(x,y),'-O')
        set(gca,'xtick',1:n_1);
        set(gca,'ytick',1:n_2);
        title('Stellar Mass')
        zlabel('log M_{star} [M_{Sun}]')
     subplot(2,2,2)
    plot3(x,y,MEDv_age(x,y), '-O')
    set(gca,'xtick',1:n_1);
    set(gca,'ytick',1:n_2);
    title('Age')
    zlabel('log ^{t}D4000 [Gyr]')
      subplot(2,2,3)
    plot3(x,y,MEDv_ssfr(x,y), '-O')
    set(gca,'xtick',1:n_1);
    set(gca,'ytick',1:n_2);
    title('SSFR')
    zlabel('log \phi [Gyr^{-1}]')
    subplot(2,2,4)
    plot3(x,y,MEDv_frac(x,y), '-O')
    set(gca,'xtick',1:n_1);
    set(gca,'ytick',1:n_2);
    title('Fraction')
    end
end
 
for h1=n_1:-1:1
 for   h2=1:1:n_2
    
   stellar_mass{h1,h2}=params(Tabv_1{h1,h2},31);
    D4000{h1,h2}=params(Tabv_1{h1,h2},37);
    ssfr{h1,h2}=params(Tabv_1{h1,h2},33)-params(Tabv_1{h1,h2},29)-9;
    
    
end
end

if (n_1*n_2 <13)
figure(12)
subplot(3,1,1)
m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
     
     m1=m1+1;
     if (type(m1) == 1) 
       errorbar(mean(ssfr{h1,h2}),mean(D4000{h1,h2}),std(D4000{h1,h2}),'sm')
       strmax = ['(',num2str(h1),',',num2str(h2),')'];
       text(mean(ssfr{h1,h2})+0.2,mean(D4000{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
       hold on
       herrorbar(mean(ssfr{h1,h2}),mean(D4000{h1,h2}),std(ssfr{h1,h2}),'sm')
     elseif (type(m1)==2) || (type(m1)==3)
        errorbar(mean(ssfr{h1,h2}),mean(D4000{h1,h2}),std(D4000{h1,h2}),'sr')
        strmax = ['(',num2str(h1),',',num2str(h2),')'];
        text(mean(ssfr{h1,h2})+0.2,mean(D4000{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
        hold on
       herrorbar(mean(ssfr{h1,h2}),mean(D4000{h1,h2}),std(ssfr{h1,h2}),'sr')
     elseif (type(m1)==4) || (type(m1)==5) || (type(m1)==6)
         errorbar(mean(ssfr{h1,h2}),mean(D4000{h1,h2}),std(D4000{h1,h2}),'sy')
         strmax = ['(',num2str(h1),',',num2str(h2),')'];
         text(mean(ssfr{h1,h2})+0.2,mean(D4000{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
         hold on
       herrorbar(mean(ssfr{h1,h2}),mean(D4000{h1,h2}),std(ssfr{h1,h2}),'sy')
     elseif (type(m1)==7) || (type(m1)==8) || (type(m1)==9)
         errorbar(mean(ssfr{h1,h2}),mean(D4000{h1,h2}),std(D4000{h1,h2}),'sb')
         strmax = ['(',num2str(h1),',',num2str(h2),')'];
         text(mean(ssfr{h1,h2})+0.2,mean(D4000{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
         hold on
       herrorbar(mean(ssfr{h1,h2}),mean(D4000{h1,h2}),std(ssfr{h1,h2}),'sb')
     else
         errorbar(mean(ssfr{h1,h2}),mean(D4000{h1,h2}),std(D4000{h1,h2}),'sk')
         strmax = ['(',num2str(h1),',',num2str(h2),')'];
         text(mean(ssfr{h1,h2})+0.2,mean(D4000{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
         hold on
       herrorbar(mean(ssfr{h1,h2}),mean(D4000{h1,h2}),std(ssfr{h1,h2}),'sk')
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
       errorbar(mean(stellar_mass{h1,h2}),mean(ssfr{h1,h2}),std(ssfr{h1,h2}),'sm')
       strmax = ['(',num2str(h1),',',num2str(h2),')'];
       text(mean(stellar_mass{h1,h2})+0.2,mean(ssfr{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
       hold on
       herrorbar(mean(stellar_mass{h1,h2}),mean(ssfr{h1,h2}),std(stellar_mass{h1,h2}),'sm')
     elseif (type(m1)==2) || (type(m1)==3)
        errorbar(mean(stellar_mass{h1,h2}),mean(ssfr{h1,h2}),std(ssfr{h1,h2}),'sr')
        strmax = ['(',num2str(h1),',',num2str(h2),')'];
        text(mean(stellar_mass{h1,h2})+0.2,mean(ssfr{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
        hold on
       herrorbar(mean(stellar_mass{h1,h2}),mean(ssfr{h1,h2}),std(stellar_mass{h1,h2}),'sr')
     elseif (type(m1)==4) || (type(m1)==5) || (type(m1)==6)
         errorbar(mean(stellar_mass{h1,h2}),mean(ssfr{h1,h2}),std(ssfr{h1,h2}),'sy')
         strmax = ['(',num2str(h1),',',num2str(h2),')'];
         text(mean(stellar_mass{h1,h2})+0.2,mean(ssfr{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
         hold on
       herrorbar(mean(stellar_mass{h1,h2}),mean(ssfr{h1,h2}),std(stellar_mass{h1,h2}),'sy')
     elseif (type(m1)==7) || (type(m1)==8) || (type(m1)==9)
         errorbar(mean(stellar_mass{h1,h2}),mean(ssfr{h1,h2}),std(ssfr{h1,h2}),'sb')
         strmax = ['(',num2str(h1),',',num2str(h2),')'];
         text(mean(stellar_mass{h1,h2})+0.2,mean(ssfr{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
         hold on
       herrorbar(mean(stellar_mass{h1,h2}),mean(ssfr{h1,h2}),std(stellar_mass{h1,h2}),'sb')
     else
         errorbar(mean(stellar_mass{h1,h2}),mean(ssfr{h1,h2}),std(ssfr{h1,h2}),'sk')
         strmax = ['(',num2str(h1),',',num2str(h2),')'];
         text(mean(stellar_mass{h1,h2})+0.2,mean(ssfr{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
         hold on
       herrorbar(mean(stellar_mass{h1,h2}),mean(ssfr{h1,h2}),std(stellar_mass{h1,h2}),'sk')
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
       errorbar(mean(stellar_mass{h1,h2}),mean(D4000{h1,h2}),std(D4000{h1,h2}),'sm')
       strmax = ['(',num2str(h1),',',num2str(h2),')'];
       text(mean(stellar_mass{h1,h2})+0.2,mean(D4000{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
       hold on 
       herrorbar(mean(stellar_mass{h1,h2}),mean(D4000{h1,h2}),std(stellar_mass{h1,h2}),'sm')
     elseif (type(m1)==2) || (type(m1)==3)
        errorbar(mean(stellar_mass{h1,h2}),mean(D4000{h1,h2}),std(D4000{h1,h2}),'sr')
        strmax = ['(',num2str(h1),',',num2str(h2),')'];
        text(mean(stellar_mass{h1,h2})+0.2,mean(D4000{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
        hold on 
       herrorbar(mean(stellar_mass{h1,h2}),mean(D4000{h1,h2}),std(stellar_mass{h1,h2}),'sr')
     elseif (type(m1)==4) || (type(m1)==5) || (type(m1)==6)
         errorbar(mean(stellar_mass{h1,h2}),mean(D4000{h1,h2}),std(D4000{h1,h2}),'sy')
         strmax = ['(',num2str(h1),',',num2str(h2),')'];
         text(mean(stellar_mass{h1,h2})+0.2,mean(D4000{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
         hold on 
       herrorbar(mean(stellar_mass{h1,h2}),mean(D4000{h1,h2}),std(stellar_mass{h1,h2}),'sy')
     elseif (type(m1)==7) || (type(m1)==8) || (type(m1)==9)
         errorbar(mean(stellar_mass{h1,h2}),mean(D4000{h1,h2}),std(D4000{h1,h2}),'sb')
         strmax = ['(',num2str(h1),',',num2str(h2),')'];
         text(mean(stellar_mass{h1,h2})+0.2,mean(D4000{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
         hold on 
       herrorbar(mean(stellar_mass{h1,h2}),mean(D4000{h1,h2}),std(stellar_mass{h1,h2}),'sb')
     else
         errorbar(mean(stellar_mass{h1,h2}),mean(D4000{h1,h2}),std(D4000{h1,h2}),'sk')
         strmax = ['(',num2str(h1),',',num2str(h2),')'];
         text(mean(stellar_mass{h1,h2})+0.2,mean(D4000{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
         hold on 
       herrorbar(mean(stellar_mass{h1,h2}),mean(D4000{h1,h2}),std(stellar_mass{h1,h2}),'sk')
     end
   hold on  
   xlabel('log M_{star} [M_{Sun}]')
   ylabel('log ^{t}D4000 [Gyr]')
 end
end
end   

rep = strcat(dir,'prop_vs_prop',n1st,'_by_self_org_res',n2st,'.png');
 saveas(figure(12),rep,'png')


%>>>>---------------------------------------------------------------------------------------------------
 net_name= strcat(dir,'net_',n1st,'_by_self_org_res_',n2st);
 save(net_name, 'net')
