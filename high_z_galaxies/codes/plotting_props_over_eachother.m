dir='~/Desktop/project/data_mining/high_z_galaxies/results_SED/1D/new_type_of_plots/new_colours/';
%%%Reading Networks
net1by2=load('~/Desktop/project/data_mining/high_z_galaxies/results_SED/1D/new_type_of_plots/net_1_by_2.mat');
net12=net1by2.net;

net1by4=load('~/Desktop/project/data_mining/high_z_galaxies/results_SED/1D/new_type_of_plots/net_1_by_4.mat');
net14=net1by4.net;

net1by8=load('~/Desktop/project/data_mining/high_z_galaxies/results_SED/1D/new_type_of_plots/net_1_by_8.mat');
net18=net1by8.net;
%%%%Reading and preparing the data adn info models
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
n_112=1;  n_212=2;
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

  m1=0;
    for h1=1:1:n_112
        for   h2=1:1:n_212
            m1=m1+1;
    
    
            MEDv_mass12(h1,h2)=median(info_model(Tabv12{h1,h2},31));
            MEDv_ssfr12(h1,h2)=median(info_model(Tabv12{h1,h2},33)-info_model(Tabv12{h1,h2},31)+9);
            MEDv_age12(h1,h2)=median(info_model(Tabv12{h1,h2},37));
            MEDv_frac12(h1,h2)=median(info_model(Tabv12{h1,h2},49));

    
        end
    end
    
    %%%% puting the data in the network
%%net14
n_114=1;  n_214=4;
sim_v14=sim(net14,annv);
for k1=1:n_114*n_214
 av14{k1}=find(sim_v14(k1,:)==1);
end


for k1=1:n_114*n_214
inpv14{k1}=annv(:,av14{k1});
end

 
m1=0;
for h1=n_114:-1:1
 for   h2=1:1:n_214
    m1=m1+1;
    Tabv14{h1,h2}=av14{m1};
 
end
end

  m1=0;
    for h1=1:1:n_114
        for   h2=1:1:n_214
            m1=m1+1;
    
    
            MEDv_mass14(h1,h2)=median(info_model(Tabv14{h1,h2},31));
            MEDv_ssfr14(h1,h2)=median(info_model(Tabv14{h1,h2},33)-info_model(Tabv14{h1,h2},31)+9);
            MEDv_age14(h1,h2)=median(info_model(Tabv14{h1,h2},37));
            MEDv_frac14(h1,h2)=median(info_model(Tabv14{h1,h2},49));

    
        end
    end
    
    %%%% puting the data in the network
%%net18
n_118=1;  n_218=8;
sim_v18=sim(net18,annv);
for k1=1:n_118*n_218
 av18{k1}=find(sim_v18(k1,:)==1);
end


for k1=1:n_118*n_218
inpv18{k1}=annv(:,av18{k1});
end

 
m1=0;
for h1=n_118:-1:1
 for   h2=1:1:n_218
    m1=m1+1;
    Tabv18{h1,h2}=av18{m1};
 
end
end

  m1=0;
    for h1=1:1:n_118
        for   h2=1:1:n_218
            m1=m1+1;
    
    
            MEDv_mass18(h1,h2)=median(info_model(Tabv18{h1,h2},31));
            MEDv_ssfr18(h1,h2)=median(info_model(Tabv18{h1,h2},33)-info_model(Tabv18{h1,h2},31)+9);
            MEDv_age18(h1,h2)=median(info_model(Tabv18{h1,h2},37));
            MEDv_frac18(h1,h2)=median(info_model(Tabv18{h1,h2},49));

    
        end
    end
    x12= 1:n_212;
    x12=(x12 - 1)/(n_212-1);
    x14= 1:n_214;
    x14=(x14 - 1)/(n_214-1);
    x18= 1:n_218;
    x18=(x18 - 1)/(n_218-1);
    
    figure(1)
    subplot(2,2,1)
    plot(x12,MEDv_mass12, 'k-O', x14,MEDv_mass14, '-o', x18,MEDv_mass18,'r-o')
    set(gca,'xtick',[0 0.25 0.5 0.75 1]);
    title('Stellar Mass')
    ylabel('log(M_{star} [M_{Sun}])')
    subplot(2,2,2)
    plot(x12,MEDv_age12, 'k-O', x14,MEDv_age14, '-o', x18,MEDv_age18,'r-o')
    set(gca,'xtick',[0 0.25 0.5 0.75 1]);
    title('Age')
    ylabel('log(^{t}D4000 [Gyr])')
    subplot(2,2,3)
    plot(x12,MEDv_ssfr12,  'k-O', x14,MEDv_ssfr14, '-o', x18,MEDv_ssfr18,'r-o')
    set(gca,'xtick',[0 0.25 0.5 0.75 1]);
    title('SSFR')
    ylabel('log(\phi [Gyr^{-1}])')
    xlabel('Normalized neuron number')
    subplot(2,2,4)
    plot(x12,MEDv_frac12,  'k-O', x14,MEDv_frac14, '-o', x18,MEDv_frac18,'r-o')
    set(gca,'xtick',[0 0.25 0.5 0.75 1]);
    ylabel('log(f_{burst})')
    xlabel('Normalized neuron number')
    title('mass fraction of young SP model')
prop= strcat(dir,'prop2_over_plot32.fig');
saveas(figure(1),prop,'fig')