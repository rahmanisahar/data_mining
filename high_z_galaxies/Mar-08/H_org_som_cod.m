clear
clc

load info_model.txt
%>>>>  load your input dara in the Excel format
%>>>> new name 'cat':  M x N  (M=regions, N= parameters; this is the format of original file)
%>>>> by this command 'Null' is converted to NaN which is readable by Matlab 

%[cat] = xlsread('m31.csv');  
%>>>>---------------------------------------------------------------------------------------------------
%>>>> if all components are numbers then you can us 'txt' format then run this to load the data and inactive the above commend by % 

load kinney_n.txt
cat=kinney_n(:,1:end);

load model_n.txt
catv=model_n(:,1:end);

%>>>>---------------------------------------------------------------------------------------------------
%>>>>  Inverse  'cat' to N x M  for normalization and other processes
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
n_1=1
n_2=8


ost=1000
tnd=2

orl=.9
tlr=0.02


n_epoch =100*n_1*n_2;
%>>>>---------------------------------------------------------------------------------------------------
%>>>> MATLAB NETWORK  ; should not change
%net = newsom(PR,[D1,D2,...],TFCN,DFCN,OLR,OSTEPS,TLR,TND)
 net = newsom(annt,[n_2,n_1],'hextop','dist',orl,ost,tlr,tnd);
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

 figure(1)
 plotsomnd(net,annt)
figure(2)
plotsomhits(net,annt)
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

%  figure(3)
%  plotsomnd(net,annv)
figure(4)
plotsomhits(net,annv)


if n_1==1

m1=0;
for h1=1:1:n_1
 for   h2=1:1:n_2
    m1=m1+1;
    
    
    MEDv_mass(h1,h2)=median(info_model(Tabv_1{h1,h2},31));
    MEDv_sfr(h1,h2)=median(info_model(Tabv_1{h1,h2},33));
    MEDv_age(h1,h2)=median(info_model(Tabv_1{h1,h2},37));
    MEDv_frac(h1,h2)=median(info_model(Tabv_1{h1,h2},49));

    
end
end


figure(5)
subplot(2,2,1)
plot(MEDv_mass, '-O')
subplot(2,2,2)
plot(MEDv_age, '-O')
subplot(2,2,3)
plot(MEDv_sfr, '-O')
subplot(2,2,4)
plot(MEDv_frac, '-O')

end
 

%>>>>---------------------------------------------------------------------------------------------------



