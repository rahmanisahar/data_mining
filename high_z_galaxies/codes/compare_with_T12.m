clear, clc

dir='~/Desktop';
%%%Reading Networks
net1by22=load('~/Desktop/project/data_mining/high_z_galaxies/results_SED/1D/new_type_of_plots/net_1_by_22.mat');
net12=net1by22.net;



load ~/Desktop/project/data_mining/high_z_galaxies/data/info_model.txt
params = info_model;

load ~/Desktop/project/data_mining/high_z_galaxies/data/model_n.txt
catv=model_n(:,1:end);

catv=catv';   
catv_fix=(catv);
catv_fix_norm= (catv_fix);
annv=catv_fix_norm';


%%%% puting the data in the network
%%net12
n_112=1;  n_212=22;
sim_v12=sim(net12,annv);
for k1=1:n_112*n_212
 av12{k1}=find(sim_v12(k1,:)==1);
end


for k1=1:n_112*n_212
inpv12{k1}=annv(:,av12{k1});
end

 
m1=0;
for h1=n_112:-1:1
 for   h2=1:1:n_212
    m1=m1+1;
    Tabv12{h1,h2}=av12{m1};
 
end
end

 h1=1;
 for   h2=1:1:n_212
    if size(Tabv12{h1,h2}) > 0
   stellar_mass{h1,h2}=params(Tabv12{h1,h2},31);
   estellar_mass{h1,h2}=params(Tabv12{h1,h2},32); 
   D4000{h1,h2}=params(Tabv12{h1,h2},37);
    eD4000{h1,h2}=params(Tabv12{h1,h2},38);
    ssfr{h1,h2}=params(Tabv12{h1,h2},33)-params(Tabv12{h1,h2},29)+9;
    Afuv{h1,h2}=params(Tabv12{h1,h2},55);
    eAfuv{h1,h2}=params(Tabv12{h1,h2},56);
    else
      stellar_mass{h1,h2}=[];
      estellar_mass{h1,h2}=[];
      D4000{h1,h2}=[];
      eD4000{h1,h2}=[];
      ssfr{h1,h2}=[];
      Afuv{h1,h2}=[];
      eAfuv{h1,h2}=[];
    end
end

    %colour1=[236 231 242]./255;
    %colour2=[43 140 190]./255;
    %colour1=[222 235 247]./255;
    %colour2=[49 130 189]./255;
    %colour2=[197 27 138]./255;
    %colour1=[253 224 221]./255;
    %colour1=[255 237 160]./255; %number 2
    %colour2=[240 59 32]./255;   %number 2
    colour1=[255 0 0]./255; %number 1
    colour2=[0 0 255]./255;   %number 1 
%     colour1=[230 230 250]./255; %number 3
%     colour2=[75 0 130]./255;   %number 3
    delta_c=colour2-colour1;

figure(1)
set(gcf,'color','white')
 h1=1;
 for   h2=1:1:n_212
     
     if size(stellar_mass{h1,h2}) > 0  
         sz=size(stellar_mass{h1,h2});
       col=[(delta_c(1)*(h2-1)/(n_212-1))+colour1(1) (delta_c(2)*(h2-1)/(n_212-1))+colour1(2) (delta_c(3)*(h2-1)/(n_212-1))+colour1(3)];
        errorbar(median(ssfr{h1,h2}),median(D4000{h1,h2}),(std(D4000{h1,h2}))/sqrt(sz(2)),'s','Color',col,...
            'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerSize',10)
       %plot(median(ssfr{h1,h2}),median(D4000{h1,h2}),'s','MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerSize',10)
      % strmax = ['(',num2str(h1),',',num2str(h2),')'];
      % text(median(ssfr{h1,h2})+0.2,median(D4000{h1,h2})+0.2,strmax,'HorizontalAlignment','right');
        hold on
       h =herrorbar(median(ssfr{h1,h2}),median(D4000{h1,h2}),(std(ssfr{h1,h2}))/sqrt(sz(2)));
       set(h,'color',col)
       hold on
       ax = gca;
       ax.LineWidth = 1.5;
       ax.FontSize = 20;
       ax.LabelFontSizeMultiplier = 1.5;
       xlabel('log \phi [Gyr^{-1}]')
      ylabel('log ^{t}D4000 [Gyr]')
    end
 end
saveas(figure(1),'~/Desktop/plot3/f1.fig','fig')
saveas(figure(1),'~/Desktop/plot3/f1.png','png')

figure(2)
set(gcf,'color','white')
h1=1;
 for   h2=1:1:n_212
     if size(stellar_mass{h1,h2}) > 0  
         sz=size(stellar_mass{h1,h2});
         col=[(delta_c(1)*(h2-1)/(n_212-1))+colour1(1) (delta_c(2)*(h2-1)/(n_212-1))+colour1(2) (delta_c(3)*(h2-1)/(n_212-1))+colour1(3)];
       errorbar(median(stellar_mass{h1,h2}),median(ssfr{h1,h2}),(std(ssfr{h1,h2}))/sqrt(sz(2)),'s','Color',col,...
           'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerSize',10)
       %plot(median(stellar_mass{h1,h2}),median(ssfr{h1,h2}),'s','MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerSize',10)
       hold on
       h11=herrorbar(median(stellar_mass{h1,h2}),median(ssfr{h1,h2}),(std(stellar_mass{h1,h2}))/sqrt(sz(2)));
       set(h11,'color',col)
       hold on 
              ax = gca;
       ax.LineWidth = 1.5;
       ax.FontSize = 20;
       ax.LabelFontSizeMultiplier = 1.5;
       xlabel('log M_{star} [M_{Sun}]')
       ylabel('log \phi [Gyr^{-1}]')
    end
 end
 
 saveas(figure(2),'~/Desktop/plot3/f2.fig','fig')
saveas(figure(2),'~/Desktop/plot3/f2.png','png')
 figure(3)
set(gcf,'color','white')

h1=1;
 for   h2=1:1:n_212
     if size(stellar_mass{h1,h2}) > 0 
         sz=size(stellar_mass{h1,h2});
         col=[(delta_c(1)*(h2-1)/(n_212-1))+colour1(1) (delta_c(2)*(h2-1)/(n_212-1))+colour1(2) (delta_c(3)*(h2-1)/(n_212-1))+colour1(3)];
       errorbar(median(stellar_mass{h1,h2}),median(D4000{h1,h2}),(std(D4000{h1,h2}))/sqrt(sz(2)),'s','Color',col,'MarkerFaceColor',col,...
           'MarkerEdgeColor',col,'MarkerSize',10)
       hold on 
       %plot(median(stellar_mass{h1,h2}),median(D4000{h1,h2}),'s','MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerSize',10)
       
       h3=herrorbar(median(stellar_mass{h1,h2}),median(D4000{h1,h2}),(std(stellar_mass{h1,h2}))/sqrt(sz(2)));
       set(h3,'color',col)
   hold on
       ax = gca;
       ax.LineWidth = 1.5;
       ax.FontSize = 20;
       ax.LabelFontSizeMultiplier = 1.5;
   xlabel('log M_{star} [M_{Sun}]')
   ylabel('log ^{t}D4000 [Gyr]')
     end
 end
 
 saveas(figure(3),'~/Desktop/plot3/f3.fig','fig')
saveas(figure(3),'~/Desktop/plot3/f3.png','png')

% figure(4)
% h1=1;
%  for   h2=1:1:n_212
%      if size(stellar_mass{h1,h2}) > 0 
%          col=[(delta_c(1)*(h2-1)/(n_212-1))+colour1(1) (delta_c(2)*(h2-1)/(n_212-1))+colour1(2) (delta_c(3)*(h2-1)/(n_212-1))+colour1(3)];
%        errorbar(stellar_mass{h1,h2},D4000{h1,h2},eD4000{h1,h2},'rs','Color',col,'MarkerFaceColor',col,...
%            'MarkerEdgeColor',col,'MarkerSize',10)
%        hold on 
%        %plot(median(stellar_mass{h1,h2}),median(D4000{h1,h2}),'s','MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerSize',10)
%        
%        h3=herrorbar(stellar_mass{h1,h2},D4000{h1,h2},estellar_mass{h1,h2});
%        set(h3,'color',col)
%    hold on  
%    xlabel('log M_{star} [M_{Sun}]')
%    ylabel('log ^{t}D4000 [Gyr]')
%      end
%  end
   
 figure(4)
 set(gcf,'color','white')
 h1=1;
 for   h2=1:1:n_212
     
     if size(stellar_mass{h1,h2}) > 0  
         sz=size(stellar_mass{h1,h2});
       col=[(delta_c(1)*(h2-1)/(n_212-1))+colour1(1) (delta_c(2)*(h2-1)/(n_212-1))+colour1(2) (delta_c(3)*(h2-1)/(n_212-1))+colour1(3)];
        errorbar(median(ssfr{h1,h2}),median(Afuv{h1,h2}),(std(Afuv{h1,h2}))/sqrt(sz(2)),'s','Color',col,...
            'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerSize',10)
        hold on
       h =herrorbar(median(ssfr{h1,h2}),median(Afuv{h1,h2}),(std(ssfr{h1,h2}))/sqrt(sz(2)));
       set(h,'color',col)
       hold on
       ax = gca;
       ax.LineWidth = 1.5;
       ax.FontSize = 20;
       ax.LabelFontSizeMultiplier = 1.5;
       xlabel('log \phi [Gyr^{-1}]')
      ylabel('A_{FUV} [mag]')
    end
 end
   
 saveas(figure(4),'~/Desktop/plot3/f4.fig','fig')
saveas(figure(4),'~/Desktop/plot3/f4.png','png')
 
  figure(5)
  set(gcf,'color','white')
 h1=1;
 for   h2=1:1:n_212
     
     if size(stellar_mass{h1,h2}) > 0 
         sz=size(stellar_mass{h1,h2});
       col=[(delta_c(1)*(h2-1)/(n_212-1))+colour1(1) (delta_c(2)*(h2-1)/(n_212-1))+colour1(2) (delta_c(3)*(h2-1)/(n_212-1))+colour1(3)];
        errorbar(median(D4000{h1,h2}),median(Afuv{h1,h2}),(std(Afuv{h1,h2}))/sqrt(sz(2)),'s','Color',col,...
            'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerSize',10)
        hold on
       h =herrorbar(median(D4000{h1,h2}),median(Afuv{h1,h2}),(std(ssfr{h1,h2}))/sqrt(sz(2)));
       set(h,'color',col)
       hold on
       ax = gca;
       ax.LineWidth = 1.5;
       ax.FontSize = 20;
       ax.LabelFontSizeMultiplier = 1.5;
       xlabel('log ^{t}D4000 [Gyr]')
      ylabel('A_{FUV} [mag]')
    end
 end

 saveas(figure(5),'~/Desktop/plot3/f5.fig','fig')
saveas(figure(5),'~/Desktop/plot3/f5.png','png')
 
  figure(6)
  set(gcf,'color','white')
 h1=1;
 for   h2=1:1:n_212
     
     if size(stellar_mass{h1,h2}) > 0
         sz=size(stellar_mass{h1,h2});
       col=[(delta_c(1)*(h2-1)/(n_212-1))+colour1(1) (delta_c(2)*(h2-1)/(n_212-1))+colour1(2) (delta_c(3)*(h2-1)/(n_212-1))+colour1(3)];
        errorbar(median(stellar_mass{h1,h2}),median(Afuv{h1,h2}),(std(Afuv{h1,h2}))/sqrt(sz(2)),'s','Color',col,...
            'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerSize',10)
        hold on
       h =herrorbar(median(stellar_mass{h1,h2}),median(Afuv{h1,h2}),(std(stellar_mass{h1,h2}))/sqrt(sz(2)));
       set(h,'color',col)
   hold on
       ax = gca;
       ax.LineWidth = 1.5;
       ax.FontSize = 20;
       ax.LabelFontSizeMultiplier = 1.5;
       xlabel('log M_{star} [M_{Sun}]')
      ylabel('A_{FUV} [mag]')
    end
 end
 
 saveas(figure(6),'~/Desktop/plot3/f6.fig','fig')
saveas(figure(6),'~/Desktop/plot3/f6.png','png')

figure(10)
h1=1; h2=3;
 col=[(delta_c(1)*(h2-1)/(n_212-1))+colour1(1) (delta_c(2)*(h2-1)/(n_212-1))+colour1(2) (delta_c(3)*(h2-1)/(n_212-1))+colour1(3)];
plot(median(stellar_mass{h1,h2}),median(D4000{h1,h2}),'s','Color',col,...
            'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerSize',10)
        hold on
h2=22;
col=[(delta_c(1)*(h2-1)/(n_212-1))+colour1(1) (delta_c(2)*(h2-1)/(n_212-1))+colour1(2) (delta_c(3)*(h2-1)/(n_212-1))+colour1(3)];
plot(median(stellar_mass{h1,h2}),median(D4000{h1,h2}),'s','Color',col,...
            'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerSize',10)
        
