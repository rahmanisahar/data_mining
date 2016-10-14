% clear, clc
n_1=1; n_2=22;
net10by10=load('~/Desktop/project/data_mining/high_z_galaxies/results_SED/2D/net_1_by_22.mat');
net=net10by10.net;
dir='~/Desktop/project/data_mining/high_z_galaxies/results_SED/2D/new_type_of_plots/new_colours/';
load ~/Desktop/project/data_mining/high_z_galaxies/data/kinney_n.txt
cat=kinney_n(:,1:end);
cat=cat';   
cat_fix=(cat);
cat_fix_norm= (cat_fix);
annt=cat_fix_norm';
sim_t=sim(net,annt);

for k1=1:n_1*n_2
 at{k1}=find(sim_t(k1,:)==1);
end


for k1=1:n_1*n_2
inpt{k1}=annt(:,at{k1});
end
%  
 figure(1)
 plotsomnd_shar(net,annt)
 figure(2)
 plotsomhits_sahar(net,annt)
 
m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
    m1=m1+1;
    Tab_2{h1,h2}=at{m1};
 
end
end


%validate test
load ~/Desktop/project/data_mining/high_z_galaxies/data/model_n.txt
catv=model_n(:,1:end);
catv=catv';   
catv_fix=(catv);
catv_fix_norm= (catv_fix);
annv=catv_fix_norm';
sim_v=sim(net,annv);
for k1=1:n_1*n_2
 av{k1}=find(sim_v(k1,:)==1);
end


for k1=1:n_1*n_2
inpv{k1}=annv(:,av{k1});
end

 figure(3)
 plotsomhits_sahar(net,annv)
 
m1=0;
for h1=n_1:-1:1
 for   h2=1:1:n_2
    m1=m1+1;
    Tabv_2{h1,h2}=av{m1};
 
end
end
% 
 n1st=int2str(n_1);
 n2st=int2str(n_2);
colorim = strcat(dir,'dist_',n1st,'_by_self_org_res_',n2st,'.fig');
hitim_t= strcat(dir,'hit_t_',n1st,'_by_self_org_res_',n2st,'.fig');
hitim_v= strcat(dir,'hit_v_',n1st,'_by_self_org_res_',n2st,'.fig');
% saveas(figure(1),colorim,'fig')
% saveas(figure(2),hitim_t,'fig')
% saveas(figure(3),hitim_v,'fig')


