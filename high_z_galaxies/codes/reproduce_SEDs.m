dir='~/Desktop/project/data_mining/high_z_galaxies/results_SED/1D/new_type_of_plots/new_colours/';
%%%Reading Networks
net1by22=load('~/Desktop/project/data_mining/high_z_galaxies/results_SED/1D/new_type_of_plots/net_1_by_22.mat');
net=net1by22.net;


load ~/Desktop/project/data_mining/high_z_galaxies/data/model_n.txt
catv=model_n(:,1:end);
load ~/Desktop/project/data_mining/nearby_galaxies/H_paper_2012_som/xx.txt
x_ax=xx;
load ~/Desktop/project/data_mining/high_z_galaxies/data/kinney_n.txt
cat=kinney_n(:,1:end);


kinney=cat;
annm=catv';
catv=catv';
cat=cat';     

cat_fix=(cat);
catv_fix=(catv);

cat_fix_norm= (cat_fix);
catv_fix_norm= (catv_fix);

annt=cat_fix_norm';
annv=catv_fix_norm';

n_1=1;  n_2=22;
sim_t=sim(net,annt);
sim_v=sim(net,annv);

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
    check_s=CATv_1{h1,h2};
               size_ch=size(check_s);
    if (size_ch(2) > 1)
    catm{h1,h2}=median(annm(av{m1},:));
    elseif (size_ch(2) == 1)
        catm{h1,h2}=annm(av{m1},:);
    end
    
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
                if (size_ch(2) > 0)
    m1=m1+1;
%     chi_s=zeros(1,12);
%     a = catm{h1,h2};
%     for n1=1:1:12
%        chi = 0;
%        for i = 1:1:630
%            chi = chi+(a(i)-kinney(i,n1))^2/kinney(i,n1);
%        end
%        chi_s(n1)=chi;
%     end   
%     [dum,num_kinney]=min(chi_s(:));
%    index=find(min(chi_s));
    %m1=m1+1;
    subplot(n_1,n_2,m1)
    plot(xx,catm{h1,h2},'K')%,xlim([1 630])
%     type(m1)=num_kinney;
   ax = gca;
       ax.LineWidth = 1;
       ax.FontSize = 14;
       ax.LabelFontSizeMultiplier = 1.5;
  % ylabel('flux (arbitary unit)')
   xlabel('Wavelength($$\AA$$)','Interpreter', 'Latex' ) 
     end    
               end
 end
end
%text(-10,10.2,'Wavelength(\AA)')

