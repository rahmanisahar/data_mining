
close all, clear, clc
net10by10=load('~/Desktop/project/data_mining/nearby_galaxies/SOM/grey_cycle/2d/subsets/without_reg1/net_without_reg1_without_stellar_mass.mat'); %Load network
net=net10by10.net;
catv = csvread('~/Desktop/project/data_mining/nearby_galaxies/m31/ascii_tables/subsets/subset1.csv',1,3);%col1,[1,col1,10,col2]); %load data
%catv = csvread('~/Desktop/project/data_mining/m101/ascii_table/m101_total_with_mean_per_arcsecsq_same_as_derived_ones.csv',1,1);
%load m101_va.txt
%catv=m101_va(:,1:end);
catv=catv';
catv_fix=fixunknowns(catv);
[catv1_fix_norm,~]= mapminmax(catv);
y_min = -1;
y_max = 1;
sz = size(catv);
catv_min = min(catv')' * ones(1,sz(2));
catv_max = max(catv')' * ones(1,sz(2));
catv_fix_norm = (y_max - y_min) * (catv - catv_min) ./ (catv_max - catv_min) + y_min;
annv=catv_fix_norm(:,1); %changing namme to introduce to network
annt=catv_fix_norm(:,2:10);
sz=size(annv); %finding size of original data
 n_1=10;
 n_2=10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Giving data to our network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_v=sim(net, annv);
sim_t=sim(net, annt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Saving informations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weights = net.IW{1,1};
numInputs = net.inputs{1}.processedSize;
numNeurons = net.layers{1}.size;
topologyFcn = net.layers{1}.topologyFcn;


neighbors = sparse(tril(net.layers{1}.distances <= 1.001) - eye(numNeurons));
neighbors2 = sparse(tril(net.layers{1}.distances >= 1.001 & net.layers{1}.distances <= 2.001 ));
neighbors3 = sparse(tril(net.layers{1}.distances >= 2.001 & net.layers{1}.distances <= 3.001 ));
neighbors4 = sparse(tril(net.layers{1}.distances >= 3.001 & net.layers{1}.distances <= 4.001 ));
neighbors5 = sparse(tril(net.layers{1}.distances >= 4.001 & net.layers{1}.distances <= 5.001 ));
neighbors6 = sparse(tril(net.layers{1}.distances >= 5.001 & net.layers{1}.distances <= 6.001 ));
neighbors7 = sparse(tril(net.layers{1}.distances >= 6.001 & net.layers{1}.distances <= 7.001 ));
neighbors8 = sparse(tril(net.layers{1}.distances >= 7.001 & net.layers{1}.distances <= 8.001 ));
neighbors9 = sparse(tril(net.layers{1}.distances >= 8.001 & net.layers{1}.distances <= 9.001 ));
neighbors10 = sparse(tril(net.layers{1}.distances >= 9.001 & net.layers{1}.distances <= 10.001 ));
neighbors11 = sparse(tril(net.layers{1}.distances >= 10.001 & net.layers{1}.distances <= 11.001 ));
neighbors12 = sparse(tril(net.layers{1}.distances >= 11.001 & net.layers{1}.distances <= 12.001 ));
neighbors13 = sparse(tril(net.layers{1}.distances >= 12.001 & net.layers{1}.distances <= 13.001 ));
neighbors14 = sparse(tril(net.layers{1}.distances >= 13.001 & net.layers{1}.distances <= 14.001 ));
neighbors15 = sparse(tril(net.layers{1}.distances >= 14.001 & net.layers{1}.distances <= 15.001 ));


k = 1;
  for i=1:numNeurons
    for j=find(neighbors(i,:))
      levels(k) = sqrt(sum((weights(i,:)-weights(j,:)).^2));
      levelss(i,j) = sqrt(sum((weights(i,:)-weights(j,:)).^2));
      k = k + 1;
    end
  end
  mm = minmax(levels);
  levels = (levels-mm(1)) ./ (mm(2)-mm(1)); %adjust weights between 0 to 1
  if mm(1) == mm(2), levels = zeros(size(levels)) + 0.5; end


for k1=1:n_1*n_2
 at{k1}=find(sim_t(k1,:)==1);
 sz=size(at{k1});
    if (sz(2) == 1)
 	at_2(k1)=at{k1}(1);
end

for k1=1:n_1*n_2
 av{k1}=find(sim_v(k1,:)==1);

 end
end

for k1=1:n_1*n_2
inpv{k1}=annv(:,av{k1});
end
 for k1=1:n_1*n_2
inpt{k1}=annt(:,at{k1});
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
    
    Tabt_1{h1,h2}=at{m1};
    
end
end

for h1=n_1:-1:1
 for   h2=1:1:n_2
    sz=size(Tabt_1{h1,h2});
    if (sz(2) == 1)
    
    Tabt(h1,h2)=Tabt_1{h1,h2}(1);
end
    
end
end


for h1=n_1:-1:1
 for   h2=1:1:n_2
    sz=size(Tabv_1{h1,h2});
    if (sz(2) == 1)
    
    Tabv(h1,h2)=Tabv_1{h1,h2}(1);

end
    
end
end

[nn1, mm1]=find(Tabv==1);
nonem_at=find(~cellfun(@isempty,at));





at_m=sparse(Tabt);
[d pred]=shortest_paths(at_m,[nn1,mm1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting and saving information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%csvwrite('~/Desktop/levels_without_stellar_mass_2.csv',(1-levels)*100); 
figure(2)
    plotsom_sahar(net,annv) %MATLAB som built-in SOM plots; shows density of each neurans

    %saveas(figure(2),'~/Desktop/project/data_mining/nearby_galaxies/SOM/grey_cycle/2d/subsets/without_reg1/where_is_reg_1_without_stellar_mass.png','png')
neighborss = sparse(net.layers{1}.distances <= 1.001 - eye(numNeurons));


  for i=1:numNeurons
    for j=find(neighborss(i,:))
    	if (i == j)
      	levelss(i,j)=0;
      else
      levelss(i,j) = sqrt(sum((weights(i,:)-weights(j,:)).^2));
      end

    end
  end


  levelss2 = (levelss-mm(1)) ./ (mm(2)-mm(1)); %adjust weights between 0 to 1
  if mm(1) == mm(2), levelss2 = zeros(size(levelss)) + 0.5; end
lev=(1-levelss2)*100;

lev(lev > 100)=100;
lev=lev/100.;



    n88=find(neighborss(88,:));
    sz=size(n88);

    %%PAth to reg 2.
t1=1; t2=1; t3=1; t4=1; t5=1; t6=1; t7=1; t8=1; t9=1; t10=1; t11=1; t12=1;
t13=1; t14=1; t15=1;
    	for i1=1:sz(2)
    		if (n88(i1) == 97)
    			display('The shortest path to reg 2 is 1');
    			weightp1_reg2(t1)=lev(88,n88(i1));
    			t1=t1+1;
    		else
    			n88_2=find(neighborss(n88(i1),:));
    			n88_2=n88_2(n88_2~=88);
    			sz2=size(n88_2);
    			for i2=1:sz2(2)
    				if (n88_2(i2) == 97)
    					display('The shortest path to reg 2 is 2');
    					weightp2_reg2(t2)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2));
    					t2=t2+1;
    				else
    					n88_3=find(neighborss(n88_2(i2),:));
    					n88_3=n88_3(n88_3~=n88_2(i2) & n88_3~=88);
    					sz3=size(n88_3);
    					for i3=1:sz3(2)
    						if (n88_3(i3) == 97)
    							display('The shortest path to reg 2 is 3');
    							weightp3_reg2(t3)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3));
    							t3=t3+1;
    						else
    							n88_4=find(neighborss(n88_3(i3),:));
    							n88_4=n88_4(n88_4~=n88_2(i2) & n88_4~=n88_3(i3) & n88_4~=88);
    							sz4=size(n88_4);
    							for i4=1:sz4(2)
    								if (n88_4(i4) == 97)
    								display('The shortest path to reg 2 is 4');
    								weightp4_reg2(t4)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4));
    								t4=t4+1;
        							else
    									n88_5=find(neighborss(n88_4(i4),:));
    									n88_5=n88_5(n88_5~=n88_2(i2) & n88_5~=n88_3(i3) & n88_5~=n88_4(i4) & n88_5~=88);
    									sz5=size(n88_5);
    									for i5=1:sz5(2)
    										if (n88_5(i5) == 97)
    										display('The shortest path to reg 2 is 5');
    										weightp5_reg2(t5)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5));
    										t5=t5+1;
    										else
    											n88_6=find(neighborss(n88_5(i5),:));
    											n88_6=n88_6(n88_6~=88 & n88_6~=n88_2(i2) & n88_6~=n88_3(i3) & n88_6~=n88_4(i4) & n88_6~=n88_5(i5));
    											sz6=size(n88_6);
    											for i6=1:sz6(2)
    												if (n88_6(i6) == 97)
    													display('The shortest path to reg 2 is 6');
    													weightp6_reg2(t6)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6));
    													t6=t6+1;
    													else
															n88_7=find(neighborss(n88_6(i6),:));
															n88_7=n88_7(n88_7~=88 & n88_7~=n88_2(i2) & n88_7~=n88_3(i3) & n88_7~=n88_4(i4) & n88_7~=n88_5(i5) & n88_7~=n88_6(i6));
															sz7=size(n88_7);
															for i7=1:sz7(2)
																if (n88_7(i7) == 97)
																	display('The shortest path to reg 2 is 7');
																	weightp7_reg2(t7)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7));
																	t7=t7+1;
																	else
																		n88_8=find(neighborss(n88_7(i7),:));
																		n88_8=n88_8(n88_8~=88 & n88_8~=n88_2(i2) & n88_8~=n88_3(i3) & n88_8~=n88_4(i4) & n88_8~=n88_5(i5) & n88_8~=n88_6(i6) & n88_8~=n88_7(i7));
																		sz8=size(n88_8);
																		for i8=1:sz8(2)
																			if (n88_8(i8) == 97)
																			display('The shortest path to reg 2 is 8');
																			weightp8_reg2(t8)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8));
																			t8=t8+1;
																			else
																				n88_9=find(neighborss(n88_8(i8),:));
																				n88_9=n88_9(n88_9~=88 & n88_9~=n88_2(i2) & n88_9~=n88_3(i3) & n88_9~=n88_4(i4) & n88_9~=n88_5(i5) & n88_9~=n88_6(i6) & n88_9~=n88_7(i7) & n88_9~=n88_8(i8));
																				sz9=size(n88_9);
																				for i9=1:sz9(2)
																					if (n88_9(i9) == 97)
																						display('The shortest path to reg 2 is 9');
																						weightp9_reg2(t9)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9));
																						t9=t9+1;
																						else
																							n88_10=find(neighborss(n88_9(i9),:));
																							n88_10=n88_10(n88_10~=88 & n88_10~=n88_2(i2) & n88_10~=n88_3(i3) & n88_10~=n88_4(i4) & n88_10~=n88_5(i5) & n88_10~=n88_6(i6) & n88_10~=n88_7(i7) & n88_10~=n88_8(i8) & n88_10~=n88_9(i9));
																							sz10=size(n88_10);
																							for i10=1:sz10(2)
																								if (n88_10(i10) == 97)
																									display('The shortest path to reg 2 is 10');
																									weightp10_reg2(t10)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10));
																									t10=t10+1;
																									else
																										n88_11=find(neighborss(n88_10(i10),:));
																										n88_11=n88_11(n88_11~=88 & n88_11~=n88_2(i2) & n88_11~=n88_3(i3) & n88_11~=n88_4(i4) & n88_11~=n88_5(i5) & n88_11~=n88_6(i6) & n88_11~=n88_7(i7) & n88_11~=n88_8(i8) & n88_11~=n88_9(i9) & n88_11~=n88_10(i10));
																										sz11=size(n88_11);
																										for i11=1:sz11(2)
																											if (n88_11(i11) == 97)
																												display('The shortest path to reg 2 is 11');
																												weightp11_reg2(t11)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11));
																												t11=t11+1;
																												else
																													n88_12=find(neighborss(n88_11(i11),:));
																													n88_12=n88_12(n88_12~=88 & n88_12~=n88_2(i2) & n88_12~=n88_3(i3) & n88_12~=n88_4(i4) & n88_12~=n88_5(i5) & n88_12~=n88_6(i6) & n88_12~=n88_7(i7) & n88_12~=n88_8(i8) & n88_12~=n88_9(i9) & n88_12~=n88_10(i10) & n88_12~=n88_11(i11));
																													sz12=size(n88_12);
																													for i12=1:sz12(2)
																														if (n88_12(i12) == 97)
																															display('The shortest path to reg 2 is 12');
																															weightp12_reg2(t12)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12));
																															t12=t12+1;
																															else
																																n88_13=find(neighborss(n88_12(i12),:));
																																n88_13=n88_13(n88_13~=88 & n88_13~=n88_2(i2) & n88_13~=n88_3(i3) & n88_13~=n88_4(i4) & n88_13~=n88_5(i5) & n88_13~=n88_6(i6) & n88_13~=n88_7(i7) & n88_13~=n88_8(i8) & n88_13~=n88_9(i9) & n88_13~=n88_10(i10) & n88_13~=n88_11(i11) & n88_13~=n88_12(i12));
																																sz13=size(n88_13);
																																for i13=1:sz13(2)
																																	if (n88_13(i13) == 97)
																																		display('The shortest path to reg 2 is 13');
																																		weightp13_reg2(t13)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n88_14=find(neighborss(n88_13(i13),:));
																																		% 	n88_14=n88_14(n88_14~=88 & n88_14~=n88_2(i2) & n88_14~=n88_3(i3) & n88_14~=n88_4(i4) & n88_14~=n88_5(i5) & n88_14~=n88_6(i6) & n88_14~=n88_7(i7) & n88_14~=n88_8(i8) & n88_14~=n88_9(i9) & n88_14~=n88_10(i10) & n88_14~=n88_11(i11) & n88_14~=n88_12(i12) & n88_14~=n88_13(i13));
																																		% 	sz14=size(n88_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n88_14(i14) == 97)
																																		% 			display('The shortest path to reg 2 is 14');
																																		% 			weightp14_reg2(t14)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n88_15=find(neighborss(n88_14(i14),:));
																												 					% 					n88_15=n88_15(n88_15~=88 & n88_15~=n88_2(i2) & n88_15~=n88_3(i3) & n88_15~=n88_4(i4) & n88_15~=n88_5(i5) & n88_15~=n88_6(i6) & n88_15~=n88_7(i7) & n88_15~=n88_8(i8) & n88_15~=n88_9(i9) & n88_15~=n88_10(i10) & n88_15~=n88_11(i11) & n88_15~=n88_12(i12) & n88_15~=n88_13(i13) & n88_15~=n88_14(i14));
																												 					% 					sz15=size(n88_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n88_15(i15) == 97)
																												 					% 							display('The shortest path to reg 2 is 15');
																												 					% 							weightp15_reg2(t15)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14))*lev(n88_14(i14),n88_15(i15));
																												 					% 							t15=t15+1;
																												 					% 						end
																																		% 				end
																																		% 			end
																																		% end
																																	end
																																end
																														end
																													end
																											end
																										end
																								end
																							end
																					end
																				end
																			end
																		end
																end
															end
                                    				end
                                				end

                                        	end
                                    	end

                                    end
                                end
                                
    						end
    					end
    				end
    			end
    		end
    	end

t1=1; t2=1; t3=1; t4=1; t5=1; t6=1; t7=1; t8=1; t9=1; t10=1; t11=1; t12=1;
t13=1; t14=1; t15=1;
loc_of_sec_reg = 91;

    	for i1=1:sz(2)
    		if (n88(i1) == loc_of_sec_reg)
    			display('The shortest path to reg 3 is 1');
    			weightp1_reg3(t1)=lev(88,n88(i1));
    			t1=t1+1;
    		else
    			n88_2=find(neighborss(n88(i1),:));
    			n88_2=n88_2(n88_2~=88);
    			sz2=size(n88_2);
    			for i2=1:sz2(2)
    				if (n88_2(i2) == loc_of_sec_reg)
    					display('The shortest path to reg 3 is 2');
    					weightp2_reg3(t2)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2));
    					t2=t2+1;
    				else
    					n88_3=find(neighborss(n88_2(i2),:));
    					n88_3=n88_3(n88_3~=n88_2(i2) & n88_3~=88);
    					sz3=size(n88_3);
    					for i3=1:sz3(2)
    						if (n88_3(i3) == loc_of_sec_reg)
    							display('The shortest path to reg 3 is 3');
    							weightp3_reg3(t3)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3));
    							t3=t3+1;
    						else
    							n88_4=find(neighborss(n88_3(i3),:));
    							n88_4=n88_4(n88_4~=n88_2(i2) & n88_4~=n88_3(i3) & n88_4~=88);
    							sz4=size(n88_4);
    							for i4=1:sz4(2)
    								if (n88_4(i4) == loc_of_sec_reg)
    								display('The shortest path to reg 3 is 4');
    								weightp4_reg3(t4)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4));
    								t4=t4+1;
        							else
    									n88_5=find(neighborss(n88_4(i4),:));
    									n88_5=n88_5(n88_5~=n88_2(i2) & n88_5~=n88_3(i3) & n88_5~=n88_4(i4) & n88_5~=88);
    									sz5=size(n88_5);
    									for i5=1:sz5(2)
    										if (n88_5(i5) == loc_of_sec_reg)
    										display('The shortest path to reg 3 is 5');
    										weightp5_reg3(t5)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5));
    										t5=t5+1;
    										else
    											n88_6=find(neighborss(n88_5(i5),:));
    											n88_6=n88_6(n88_6~=88 & n88_6~=n88_2(i2) & n88_6~=n88_3(i3) & n88_6~=n88_4(i4) & n88_6~=n88_5(i5));
    											sz6=size(n88_6);
    											for i6=1:sz6(2)
    												if (n88_6(i6) == loc_of_sec_reg)
    													display('The shortest path to reg 3 is 6');
    													weightp6_reg3(t6)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6));
    													t6=t6+1;
    													else
															n88_7=find(neighborss(n88_6(i6),:));
															n88_7=n88_7(n88_7~=88 & n88_7~=n88_2(i2) & n88_7~=n88_3(i3) & n88_7~=n88_4(i4) & n88_7~=n88_5(i5) & n88_7~=n88_6(i6));
															sz7=size(n88_7);
															for i7=1:sz7(2)
																if (n88_7(i7) == loc_of_sec_reg)
																	display('The shortest path to reg 3 is 7');
																	weightp7_reg3(t7)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7));
																	t7=t7+1;
																	else
																		n88_8=find(neighborss(n88_7(i7),:));
																		n88_8=n88_8(n88_8~=88 & n88_8~=n88_2(i2) & n88_8~=n88_3(i3) & n88_8~=n88_4(i4) & n88_8~=n88_5(i5) & n88_8~=n88_6(i6) & n88_8~=n88_7(i7));
																		sz8=size(n88_8);
																		for i8=1:sz8(2)
																			if (n88_8(i8) == loc_of_sec_reg)
																			display('The shortest path to reg 3 is 8');
																			weightp8_reg3(t8)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8));
																			t8=t8+1;
																			else
																				n88_9=find(neighborss(n88_8(i8),:));
																				n88_9=n88_9(n88_9~=88 & n88_9~=n88_2(i2) & n88_9~=n88_3(i3) & n88_9~=n88_4(i4) & n88_9~=n88_5(i5) & n88_9~=n88_6(i6) & n88_9~=n88_7(i7) & n88_9~=n88_8(i8));
																				sz9=size(n88_9);
																				for i9=1:sz9(2)
																					if (n88_9(i9) == loc_of_sec_reg)
																						display('The shortest path to reg 3 is 9');
																						weightp9_reg3(t9)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9));
																						t9=t9+1;
																						else
																							n88_10=find(neighborss(n88_9(i9),:));
																							n88_10=n88_10(n88_10~=88 & n88_10~=n88_2(i2) & n88_10~=n88_3(i3) & n88_10~=n88_4(i4) & n88_10~=n88_5(i5) & n88_10~=n88_6(i6) & n88_10~=n88_7(i7) & n88_10~=n88_8(i8) & n88_10~=n88_9(i9));
																							sz10=size(n88_10);
																							for i10=1:sz10(2)
																								if (n88_10(i10) == loc_of_sec_reg)
																									display('The shortest path to reg 3 is 10');
																									weightp10_reg3(t10)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10));
																									t10=t10+1;
																									else
																										n88_11=find(neighborss(n88_10(i10),:));
																										n88_11=n88_11(n88_11~=88 & n88_11~=n88_2(i2) & n88_11~=n88_3(i3) & n88_11~=n88_4(i4) & n88_11~=n88_5(i5) & n88_11~=n88_6(i6) & n88_11~=n88_7(i7) & n88_11~=n88_8(i8) & n88_11~=n88_9(i9) & n88_11~=n88_10(i10));
																										sz11=size(n88_11);
																										for i11=1:sz11(2)
																											if (n88_11(i11) == loc_of_sec_reg)
																												display('The shortest path to reg 3 is 11');
																												weightp11_reg3(t11)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11));
																												t11=t11+1;
																												else
																													n88_12=find(neighborss(n88_11(i11),:));
																													n88_12=n88_12(n88_12~=88 & n88_12~=n88_2(i2) & n88_12~=n88_3(i3) & n88_12~=n88_4(i4) & n88_12~=n88_5(i5) & n88_12~=n88_6(i6) & n88_12~=n88_7(i7) & n88_12~=n88_8(i8) & n88_12~=n88_9(i9) & n88_12~=n88_10(i10) & n88_12~=n88_11(i11));
																													sz12=size(n88_12);
																													for i12=1:sz12(2)
																														if (n88_12(i12) == loc_of_sec_reg)
																															display('The shortest path to reg 3 is 12');
																															weightp12_reg3(t12)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12));
																															t12=t12+1;
																															else
																																n88_13=find(neighborss(n88_12(i12),:));
																																n88_13=n88_13(n88_13~=88 & n88_13~=n88_2(i2) & n88_13~=n88_3(i3) & n88_13~=n88_4(i4) & n88_13~=n88_5(i5) & n88_13~=n88_6(i6) & n88_13~=n88_7(i7) & n88_13~=n88_8(i8) & n88_13~=n88_9(i9) & n88_13~=n88_10(i10) & n88_13~=n88_11(i11) & n88_13~=n88_12(i12));
																																sz13=size(n88_13);
																																for i13=1:sz13(2)
																																	if (n88_13(i13) == loc_of_sec_reg)
																																		display('The shortest path to reg 3 is 13');
																																		weightp13_reg3(t13)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n88_14=find(neighborss(n88_13(i13),:));
																																		% 	n88_14=n88_14(n88_14~=88 & n88_14~=n88_2(i2) & n88_14~=n88_3(i3) & n88_14~=n88_4(i4) & n88_14~=n88_5(i5) & n88_14~=n88_6(i6) & n88_14~=n88_7(i7) & n88_14~=n88_8(i8) & n88_14~=n88_9(i9) & n88_14~=n88_10(i10) & n88_14~=n88_11(i11) & n88_14~=n88_12(i12) & n88_14~=n88_13(i13));
																																		% 	sz14=size(n88_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n88_14(i14) == loc_of_sec_reg)
																																		% 			display('The shortest path to reg 3 is 14');
																																		% 			weightp14_reg3(t14)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n88_15=find(neighborss(n88_14(i14),:));
																												 					% 					n88_15=n88_15(n88_15~=88 & n88_15~=n88_2(i2) & n88_15~=n88_3(i3) & n88_15~=n88_4(i4) & n88_15~=n88_5(i5) & n88_15~=n88_6(i6) & n88_15~=n88_7(i7) & n88_15~=n88_8(i8) & n88_15~=n88_9(i9) & n88_15~=n88_10(i10) & n88_15~=n88_11(i11) & n88_15~=n88_12(i12) & n88_15~=n88_13(i13) & n88_15~=n88_14(i14));
																												 					% 					sz15=size(n88_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n88_15(i15) == loc_of_sec_reg)
																												 					% 							display('The shortest path to reg 3 is 15');
																												 					% 							weightp15_reg3(t15)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14))*lev(n88_14(i14),n88_15(i15));
																												 					% 							t15=t15+1;
																												 					% 						end
																																		% 				end
																																		% 			end
																																		% end
																																	end
																																end
																														end
																													end
																											end
																										end
																								end
																							end
																					end
																				end
																			end
																		end
																end
															end
                                    				end
                                				end

                                        	end
                                    	end

                                    end
                                end
                                
    						end
    					end
    				end
    			end
    		end
    	end
    
  
t1=1; t2=1; t3=1; t4=1; t5=1; t6=1; t7=1; t8=1; t9=1; t10=1; t11=1; t12=1;
t13=1; t14=1; t15=1;
loc_of_sec_reg = 46;

    	for i1=1:sz(2)
    		if (n88(i1) == loc_of_sec_reg)
    			display('The shortest path to reg 4 is 1');
    			weightp1_reg4(t1)=lev(88,n88(i1));
    			t1=t1+1;
    		else
    			n88_2=find(neighborss(n88(i1),:));
    			n88_2=n88_2(n88_2~=88);
    			sz2=size(n88_2);
    			for i2=1:sz2(2)
    				if (n88_2(i2) == loc_of_sec_reg)
    					display('The shortest path to reg 4 is 2');
    					weightp2_reg4(t2)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2));
    					t2=t2+1;
    				else
    					n88_3=find(neighborss(n88_2(i2),:));
    					n88_3=n88_3(n88_3~=n88_2(i2) & n88_3~=88);
    					sz3=size(n88_3);
    					for i3=1:sz3(2)
    						if (n88_3(i3) == loc_of_sec_reg)
    							display('The shortest path to reg 4 is 3');
    							weightp3_reg4(t3)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3));
    							t3=t3+1;
    						else
    							n88_4=find(neighborss(n88_3(i3),:));
    							n88_4=n88_4(n88_4~=n88_2(i2) & n88_4~=n88_3(i3) & n88_4~=88);
    							sz4=size(n88_4);
    							for i4=1:sz4(2)
    								if (n88_4(i4) == loc_of_sec_reg)
    								display('The shortest path to reg 4 is 4');
    								weightp4_reg4(t4)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4));
    								t4=t4+1;
        							else
    									n88_5=find(neighborss(n88_4(i4),:));
    									n88_5=n88_5(n88_5~=n88_2(i2) & n88_5~=n88_3(i3) & n88_5~=n88_4(i4) & n88_5~=88);
    									sz5=size(n88_5);
    									for i5=1:sz5(2)
    										if (n88_5(i5) == loc_of_sec_reg)
    										display('The shortest path to reg 4 is 5');
    										weightp5_reg4(t5)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5));
    										t5=t5+1;
    										else
    											n88_6=find(neighborss(n88_5(i5),:));
    											n88_6=n88_6(n88_6~=88 & n88_6~=n88_2(i2) & n88_6~=n88_3(i3) & n88_6~=n88_4(i4) & n88_6~=n88_5(i5));
    											sz6=size(n88_6);
    											for i6=1:sz6(2)
    												if (n88_6(i6) == loc_of_sec_reg)
    													display('The shortest path to reg 4 is 6');
    													weightp6_reg4(t6)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6));
    													t6=t6+1;
    													else
															n88_7=find(neighborss(n88_6(i6),:));
															n88_7=n88_7(n88_7~=88 & n88_7~=n88_2(i2) & n88_7~=n88_3(i3) & n88_7~=n88_4(i4) & n88_7~=n88_5(i5) & n88_7~=n88_6(i6));
															sz7=size(n88_7);
															for i7=1:sz7(2)
																if (n88_7(i7) == loc_of_sec_reg)
																	display('The shortest path to reg 4 is 7');
																	weightp7_reg4(t7)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7));
																	t7=t7+1;
																	else
																		n88_8=find(neighborss(n88_7(i7),:));
																		n88_8=n88_8(n88_8~=88 & n88_8~=n88_2(i2) & n88_8~=n88_3(i3) & n88_8~=n88_4(i4) & n88_8~=n88_5(i5) & n88_8~=n88_6(i6) & n88_8~=n88_7(i7));
																		sz8=size(n88_8);
																		for i8=1:sz8(2)
																			if (n88_8(i8) == loc_of_sec_reg)
																			display('The shortest path to reg 4 is 8');
																			weightp8_reg4(t8)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8));
																			t8=t8+1;
																			else
																				n88_9=find(neighborss(n88_8(i8),:));
																				n88_9=n88_9(n88_9~=88 & n88_9~=n88_2(i2) & n88_9~=n88_3(i3) & n88_9~=n88_4(i4) & n88_9~=n88_5(i5) & n88_9~=n88_6(i6) & n88_9~=n88_7(i7) & n88_9~=n88_8(i8));
																				sz9=size(n88_9);
																				for i9=1:sz9(2)
																					if (n88_9(i9) == loc_of_sec_reg)
																						display('The shortest path to reg 4 is 9');
																						weightp9_reg4(t9)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9));
																						t9=t9+1;
																						else
																							n88_10=find(neighborss(n88_9(i9),:));
																							n88_10=n88_10(n88_10~=88 & n88_10~=n88_2(i2) & n88_10~=n88_3(i3) & n88_10~=n88_4(i4) & n88_10~=n88_5(i5) & n88_10~=n88_6(i6) & n88_10~=n88_7(i7) & n88_10~=n88_8(i8) & n88_10~=n88_9(i9));
																							sz10=size(n88_10);
																							for i10=1:sz10(2)
																								if (n88_10(i10) == loc_of_sec_reg)
																									display('The shortest path to reg 4 is 10');
																									weightp10_reg4(t10)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10));
																									t10=t10+1;
																									else
																										n88_11=find(neighborss(n88_10(i10),:));
																										n88_11=n88_11(n88_11~=88 & n88_11~=n88_2(i2) & n88_11~=n88_3(i3) & n88_11~=n88_4(i4) & n88_11~=n88_5(i5) & n88_11~=n88_6(i6) & n88_11~=n88_7(i7) & n88_11~=n88_8(i8) & n88_11~=n88_9(i9) & n88_11~=n88_10(i10));
																										sz11=size(n88_11);
																										for i11=1:sz11(2)
																											if (n88_11(i11) == loc_of_sec_reg)
																												display('The shortest path to reg 4 is 11');
																												weightp11_reg4(t11)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11));
																												t11=t11+1;
																												else
																													n88_12=find(neighborss(n88_11(i11),:));
																													n88_12=n88_12(n88_12~=88 & n88_12~=n88_2(i2) & n88_12~=n88_3(i3) & n88_12~=n88_4(i4) & n88_12~=n88_5(i5) & n88_12~=n88_6(i6) & n88_12~=n88_7(i7) & n88_12~=n88_8(i8) & n88_12~=n88_9(i9) & n88_12~=n88_10(i10) & n88_12~=n88_11(i11));
																													sz12=size(n88_12);
																													for i12=1:sz12(2)
																														if (n88_12(i12) == loc_of_sec_reg)
																															display('The shortest path to reg 4 is 12');
																															weightp12_reg4(t12)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12));
																															t12=t12+1;
																															else
																																n88_13=find(neighborss(n88_12(i12),:));
																																n88_13=n88_13(n88_13~=88 & n88_13~=n88_2(i2) & n88_13~=n88_3(i3) & n88_13~=n88_4(i4) & n88_13~=n88_5(i5) & n88_13~=n88_6(i6) & n88_13~=n88_7(i7) & n88_13~=n88_8(i8) & n88_13~=n88_9(i9) & n88_13~=n88_10(i10) & n88_13~=n88_11(i11) & n88_13~=n88_12(i12));
																																sz13=size(n88_13);
																																for i13=1:sz13(2)
																																	if (n88_13(i13) == loc_of_sec_reg)
																																		display('The shortest path to reg 4 is 13');
																																		weightp13_reg4(t13)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n88_14=find(neighborss(n88_13(i13),:));
																																		% 	n88_14=n88_14(n88_14~=88 & n88_14~=n88_2(i2) & n88_14~=n88_3(i3) & n88_14~=n88_4(i4) & n88_14~=n88_5(i5) & n88_14~=n88_6(i6) & n88_14~=n88_7(i7) & n88_14~=n88_8(i8) & n88_14~=n88_9(i9) & n88_14~=n88_10(i10) & n88_14~=n88_11(i11) & n88_14~=n88_12(i12) & n88_14~=n88_13(i13));
																																		% 	sz14=size(n88_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n88_14(i14) == loc_of_sec_reg)
																																		% 			display('The shortest path to reg 4 is 14');
																																		% 			weightp14_reg4(t14)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n88_15=find(neighborss(n88_14(i14),:));
																												 					% 					n88_15=n88_15(n88_15~=88 & n88_15~=n88_2(i2) & n88_15~=n88_3(i3) & n88_15~=n88_4(i4) & n88_15~=n88_5(i5) & n88_15~=n88_6(i6) & n88_15~=n88_7(i7) & n88_15~=n88_8(i8) & n88_15~=n88_9(i9) & n88_15~=n88_10(i10) & n88_15~=n88_11(i11) & n88_15~=n88_12(i12) & n88_15~=n88_13(i13) & n88_15~=n88_14(i14));
																												 					% 					sz15=size(n88_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n88_15(i15) == loc_of_sec_reg)
																												 					% 							display('The shortest path to reg 4 is 15');
																												 					% 							weightp15_reg4(t15)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14))*lev(n88_14(i14),n88_15(i15));
																												 					% 							t15=t15+1;
																												 					% 						end
																																		% 				end
																																		% 			end
																																		% end
																																	end
																																end
																														end
																													end
																											end
																										end
																								end
																							end
																					end
																				end
																			end
																		end
																end
															end
                                    				end
                                				end

                                        	end
                                    	end

                                    end
                                end
                                
    						end
    					end
    				end
    			end
    		end
    	end
    
  
t1=1; t2=1; t3=1; t4=1; t5=1; t6=1; t7=1; t8=1; t9=1; t10=1; t11=1; t12=1;
t13=1; t14=1; t15=1;
loc_of_sec_reg = 5;

    	for i1=1:sz(2)
    		if (n88(i1) == loc_of_sec_reg)
    			display('The shortest path to reg 5 is 1');
    			weightp1_reg5(t1)=lev(88,n88(i1));
    			t1=t1+1;
    		else
    			n88_2=find(neighborss(n88(i1),:));
    			n88_2=n88_2(n88_2~=88);
    			sz2=size(n88_2);
    			for i2=1:sz2(2)
    				if (n88_2(i2) == loc_of_sec_reg)
    					display('The shortest path to reg 5 is 2');
    					weightp2_reg5(t2)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2));
    					t2=t2+1;
    				else
    					n88_3=find(neighborss(n88_2(i2),:));
    					n88_3=n88_3(n88_3~=n88_2(i2) & n88_3~=88);
    					sz3=size(n88_3);
    					for i3=1:sz3(2)
    						if (n88_3(i3) == loc_of_sec_reg)
    							display('The shortest path to reg 5 is 3');
    							weightp3_reg5(t3)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3));
    							t3=t3+1;
    						else
    							n88_4=find(neighborss(n88_3(i3),:));
    							n88_4=n88_4(n88_4~=n88_2(i2) & n88_4~=n88_3(i3) & n88_4~=88);
    							sz4=size(n88_4);
    							for i4=1:sz4(2)
    								if (n88_4(i4) == loc_of_sec_reg)
    								display('The shortest path to reg 5 is 4');
    								weightp4_reg5(t4)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4));
    								t4=t4+1;
        							else
    									n88_5=find(neighborss(n88_4(i4),:));
    									n88_5=n88_5(n88_5~=n88_2(i2) & n88_5~=n88_3(i3) & n88_5~=n88_4(i4) & n88_5~=88);
    									sz5=size(n88_5);
    									for i5=1:sz5(2)
    										if (n88_5(i5) == loc_of_sec_reg)
    										display('The shortest path to reg 5 is 5');
    										weightp5_reg5(t5)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5));
    										t5=t5+1;
    										else
    											n88_6=find(neighborss(n88_5(i5),:));
    											n88_6=n88_6(n88_6~=88 & n88_6~=n88_2(i2) & n88_6~=n88_3(i3) & n88_6~=n88_4(i4) & n88_6~=n88_5(i5));
    											sz6=size(n88_6);
    											for i6=1:sz6(2)
    												if (n88_6(i6) == loc_of_sec_reg)
    													display('The shortest path to reg 5 is 6');
    													weightp6_reg5(t6)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6));
    													t6=t6+1;
    													else
															n88_7=find(neighborss(n88_6(i6),:));
															n88_7=n88_7(n88_7~=88 & n88_7~=n88_2(i2) & n88_7~=n88_3(i3) & n88_7~=n88_4(i4) & n88_7~=n88_5(i5) & n88_7~=n88_6(i6));
															sz7=size(n88_7);
															for i7=1:sz7(2)
																if (n88_7(i7) == loc_of_sec_reg)
																	display('The shortest path to reg 5 is 7');
																	weightp7_reg5(t7)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7));
																	t7=t7+1;
																	else
																		n88_8=find(neighborss(n88_7(i7),:));
																		n88_8=n88_8(n88_8~=88 & n88_8~=n88_2(i2) & n88_8~=n88_3(i3) & n88_8~=n88_4(i4) & n88_8~=n88_5(i5) & n88_8~=n88_6(i6) & n88_8~=n88_7(i7));
																		sz8=size(n88_8);
																		for i8=1:sz8(2)
																			if (n88_8(i8) == loc_of_sec_reg)
																			display('The shortest path to reg 5 is 8');
																			weightp8_reg5(t8)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8));
																			t8=t8+1;
																			else
																				n88_9=find(neighborss(n88_8(i8),:));
																				n88_9=n88_9(n88_9~=88 & n88_9~=n88_2(i2) & n88_9~=n88_3(i3) & n88_9~=n88_4(i4) & n88_9~=n88_5(i5) & n88_9~=n88_6(i6) & n88_9~=n88_7(i7) & n88_9~=n88_8(i8));
																				sz9=size(n88_9);
																				for i9=1:sz9(2)
																					if (n88_9(i9) == loc_of_sec_reg)
																						display('The shortest path to reg 5 is 9');
																						weightp9_reg5(t9)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9));
																						t9=t9+1;
																						else
																							n88_10=find(neighborss(n88_9(i9),:));
																							n88_10=n88_10(n88_10~=88 & n88_10~=n88_2(i2) & n88_10~=n88_3(i3) & n88_10~=n88_4(i4) & n88_10~=n88_5(i5) & n88_10~=n88_6(i6) & n88_10~=n88_7(i7) & n88_10~=n88_8(i8) & n88_10~=n88_9(i9));
																							sz10=size(n88_10);
																							for i10=1:sz10(2)
																								if (n88_10(i10) == loc_of_sec_reg)
																									display('The shortest path to reg 5 is 10');
																									weightp10_reg5(t10)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10));
																									t10=t10+1;
																									else
																										n88_11=find(neighborss(n88_10(i10),:));
																										n88_11=n88_11(n88_11~=88 & n88_11~=n88_2(i2) & n88_11~=n88_3(i3) & n88_11~=n88_4(i4) & n88_11~=n88_5(i5) & n88_11~=n88_6(i6) & n88_11~=n88_7(i7) & n88_11~=n88_8(i8) & n88_11~=n88_9(i9) & n88_11~=n88_10(i10));
																										sz11=size(n88_11);
																										for i11=1:sz11(2)
																											if (n88_11(i11) == loc_of_sec_reg)
																												display('The shortest path to reg 5 is 11');
																												weightp11_reg5(t11)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11));
																												t11=t11+1;
																												else
																													n88_12=find(neighborss(n88_11(i11),:));
																													n88_12=n88_12(n88_12~=88 & n88_12~=n88_2(i2) & n88_12~=n88_3(i3) & n88_12~=n88_4(i4) & n88_12~=n88_5(i5) & n88_12~=n88_6(i6) & n88_12~=n88_7(i7) & n88_12~=n88_8(i8) & n88_12~=n88_9(i9) & n88_12~=n88_10(i10) & n88_12~=n88_11(i11));
																													sz12=size(n88_12);
																													for i12=1:sz12(2)
																														if (n88_12(i12) == loc_of_sec_reg)
																															display('The shortest path to reg 5 is 12');
																															weightp12_reg5(t12)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12));
																															t12=t12+1;
																															else
																																n88_13=find(neighborss(n88_12(i12),:));
																																n88_13=n88_13(n88_13~=88 & n88_13~=n88_2(i2) & n88_13~=n88_3(i3) & n88_13~=n88_4(i4) & n88_13~=n88_5(i5) & n88_13~=n88_6(i6) & n88_13~=n88_7(i7) & n88_13~=n88_8(i8) & n88_13~=n88_9(i9) & n88_13~=n88_10(i10) & n88_13~=n88_11(i11) & n88_13~=n88_12(i12));
																																sz13=size(n88_13);
																																for i13=1:sz13(2)
																																	if (n88_13(i13) == loc_of_sec_reg)
																																		display('The shortest path to reg 5 is 13');
																																		weightp13_reg5(t13)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n88_14=find(neighborss(n88_13(i13),:));
																																		% 	n88_14=n88_14(n88_14~=88 & n88_14~=n88_2(i2) & n88_14~=n88_3(i3) & n88_14~=n88_4(i4) & n88_14~=n88_5(i5) & n88_14~=n88_6(i6) & n88_14~=n88_7(i7) & n88_14~=n88_8(i8) & n88_14~=n88_9(i9) & n88_14~=n88_10(i10) & n88_14~=n88_11(i11) & n88_14~=n88_12(i12) & n88_14~=n88_13(i13));
																																		% 	sz14=size(n88_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n88_14(i14) == loc_of_sec_reg)
																																		% 			display('The shortest path to reg 5 is 14');
																																		% 			weightp14_reg5(t14)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n88_15=find(neighborss(n88_14(i14),:));
																												 					% 					n88_15=n88_15(n88_15~=88 & n88_15~=n88_2(i2) & n88_15~=n88_3(i3) & n88_15~=n88_4(i4) & n88_15~=n88_5(i5) & n88_15~=n88_6(i6) & n88_15~=n88_7(i7) & n88_15~=n88_8(i8) & n88_15~=n88_9(i9) & n88_15~=n88_10(i10) & n88_15~=n88_11(i11) & n88_15~=n88_12(i12) & n88_15~=n88_13(i13) & n88_15~=n88_14(i14));
																												 					% 					sz15=size(n88_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n88_15(i15) == loc_of_sec_reg)
																												 					% 							display('The shortest path to reg 5 is 15');
																												 					% 							weightp15_reg5(t15)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14))*lev(n88_14(i14),n88_15(i15));
																												 					% 							t15=t15+1;
																												 					% 						end
																																		% 				end
																																		% 			end
																																		% end
																																	end
																																end
																														end
																													end
																											end
																										end
																								end
																							end
																					end
																				end
																			end
																		end
																end
															end
                                    				end
                                				end

                                        	end
                                    	end

                                    end
                                end
                                
    						end
    					end
    				end
    			end
    		end
    	end
    
  
t1=1; t2=1; t3=1; t4=1; t5=1; t6=1; t7=1; t8=1; t9=1; t10=1; t11=1; t12=1;
t13=1; t14=1; t15=1;
loc_of_sec_reg = 1;

    	for i1=1:sz(2)
    		if (n88(i1) == loc_of_sec_reg)
    			display('The shortest path to reg 6 is 1');
    			weightp1_reg6(t1)=lev(88,n88(i1));
    			t1=t1+1;
    		else
    			n88_2=find(neighborss(n88(i1),:));
    			n88_2=n88_2(n88_2~=88);
    			sz2=size(n88_2);
    			for i2=1:sz2(2)
    				if (n88_2(i2) == loc_of_sec_reg)
    					display('The shortest path to reg 6 is 2');
    					weightp2_reg6(t2)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2));
    					t2=t2+1;
    				else
    					n88_3=find(neighborss(n88_2(i2),:));
    					n88_3=n88_3(n88_3~=n88_2(i2) & n88_3~=88);
    					sz3=size(n88_3);
    					for i3=1:sz3(2)
    						if (n88_3(i3) == loc_of_sec_reg)
    							display('The shortest path to reg 6 is 3');
    							weightp3_reg6(t3)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3));
    							t3=t3+1;
    						else
    							n88_4=find(neighborss(n88_3(i3),:));
    							n88_4=n88_4(n88_4~=n88_2(i2) & n88_4~=n88_3(i3) & n88_4~=88);
    							sz4=size(n88_4);
    							for i4=1:sz4(2)
    								if (n88_4(i4) == loc_of_sec_reg)
    								display('The shortest path to reg 6 is 4');
    								weightp4_reg6(t4)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4));
    								t4=t4+1;
        							else
    									n88_5=find(neighborss(n88_4(i4),:));
    									n88_5=n88_5(n88_5~=n88_2(i2) & n88_5~=n88_3(i3) & n88_5~=n88_4(i4) & n88_5~=88);
    									sz5=size(n88_5);
    									for i5=1:sz5(2)
    										if (n88_5(i5) == loc_of_sec_reg)
    										display('The shortest path to reg 6 is 5');
    										weightp5_reg6(t5)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5));
    										t5=t5+1;
    										else
    											n88_6=find(neighborss(n88_5(i5),:));
    											n88_6=n88_6(n88_6~=88 & n88_6~=n88_2(i2) & n88_6~=n88_3(i3) & n88_6~=n88_4(i4) & n88_6~=n88_5(i5));
    											sz6=size(n88_6);
    											for i6=1:sz6(2)
    												if (n88_6(i6) == loc_of_sec_reg)
    													display('The shortest path to reg 6 is 6');
    													weightp6_reg6(t6)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6));
    													t6=t6+1;
    													else
															n88_7=find(neighborss(n88_6(i6),:));
															n88_7=n88_7(n88_7~=88 & n88_7~=n88_2(i2) & n88_7~=n88_3(i3) & n88_7~=n88_4(i4) & n88_7~=n88_5(i5) & n88_7~=n88_6(i6));
															sz7=size(n88_7);
															for i7=1:sz7(2)
																if (n88_7(i7) == loc_of_sec_reg)
																	display('The shortest path to reg 6 is 7');
																	weightp7_reg6(t7)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7));
																	t7=t7+1;
																	else
																		n88_8=find(neighborss(n88_7(i7),:));
																		n88_8=n88_8(n88_8~=88 & n88_8~=n88_2(i2) & n88_8~=n88_3(i3) & n88_8~=n88_4(i4) & n88_8~=n88_5(i5) & n88_8~=n88_6(i6) & n88_8~=n88_7(i7));
																		sz8=size(n88_8);
																		for i8=1:sz8(2)
																			if (n88_8(i8) == loc_of_sec_reg)
																			display('The shortest path to reg 6 is 8');
																			weightp8_reg6(t8)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8));
																			t8=t8+1;
																			else
																				n88_9=find(neighborss(n88_8(i8),:));
																				n88_9=n88_9(n88_9~=88 & n88_9~=n88_2(i2) & n88_9~=n88_3(i3) & n88_9~=n88_4(i4) & n88_9~=n88_5(i5) & n88_9~=n88_6(i6) & n88_9~=n88_7(i7) & n88_9~=n88_8(i8));
																				sz9=size(n88_9);
																				for i9=1:sz9(2)
																					if (n88_9(i9) == loc_of_sec_reg)
																						display('The shortest path to reg 6 is 9');
																						weightp9_reg6(t9)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9));
																						t9=t9+1;
																						else
																							n88_10=find(neighborss(n88_9(i9),:));
																							n88_10=n88_10(n88_10~=88 & n88_10~=n88_2(i2) & n88_10~=n88_3(i3) & n88_10~=n88_4(i4) & n88_10~=n88_5(i5) & n88_10~=n88_6(i6) & n88_10~=n88_7(i7) & n88_10~=n88_8(i8) & n88_10~=n88_9(i9));
																							sz10=size(n88_10);
																							for i10=1:sz10(2)
																								if (n88_10(i10) == loc_of_sec_reg)
																									display('The shortest path to reg 6 is 10');
																									weightp10_reg6(t10)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10));
																									t10=t10+1;
																									else
																										n88_11=find(neighborss(n88_10(i10),:));
																										n88_11=n88_11(n88_11~=88 & n88_11~=n88_2(i2) & n88_11~=n88_3(i3) & n88_11~=n88_4(i4) & n88_11~=n88_5(i5) & n88_11~=n88_6(i6) & n88_11~=n88_7(i7) & n88_11~=n88_8(i8) & n88_11~=n88_9(i9) & n88_11~=n88_10(i10));
																										sz11=size(n88_11);
																										for i11=1:sz11(2)
																											if (n88_11(i11) == loc_of_sec_reg)
																												display('The shortest path to reg 6 is 11');
																												weightp11_reg6(t11)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11));
																												t11=t11+1;
																												else
																													n88_12=find(neighborss(n88_11(i11),:));
																													n88_12=n88_12(n88_12~=88 & n88_12~=n88_2(i2) & n88_12~=n88_3(i3) & n88_12~=n88_4(i4) & n88_12~=n88_5(i5) & n88_12~=n88_6(i6) & n88_12~=n88_7(i7) & n88_12~=n88_8(i8) & n88_12~=n88_9(i9) & n88_12~=n88_10(i10) & n88_12~=n88_11(i11));
																													sz12=size(n88_12);
																													for i12=1:sz12(2)
																														if (n88_12(i12) == loc_of_sec_reg)
																															display('The shortest path to reg 6 is 12');
																															weightp12_reg6(t12)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12));
																															t12=t12+1;
																															else
																																n88_13=find(neighborss(n88_12(i12),:));
																																n88_13=n88_13(n88_13~=88 & n88_13~=n88_2(i2) & n88_13~=n88_3(i3) & n88_13~=n88_4(i4) & n88_13~=n88_5(i5) & n88_13~=n88_6(i6) & n88_13~=n88_7(i7) & n88_13~=n88_8(i8) & n88_13~=n88_9(i9) & n88_13~=n88_10(i10) & n88_13~=n88_11(i11) & n88_13~=n88_12(i12));
																																sz13=size(n88_13);
																																for i13=1:sz13(2)
																																	if (n88_13(i13) == loc_of_sec_reg)
																																		display('The shortest path to reg 6 is 13');
																																		weightp13_reg6(t13)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n88_14=find(neighborss(n88_13(i13),:));
																																		% 	n88_14=n88_14(n88_14~=88 & n88_14~=n88_2(i2) & n88_14~=n88_3(i3) & n88_14~=n88_4(i4) & n88_14~=n88_5(i5) & n88_14~=n88_6(i6) & n88_14~=n88_7(i7) & n88_14~=n88_8(i8) & n88_14~=n88_9(i9) & n88_14~=n88_10(i10) & n88_14~=n88_11(i11) & n88_14~=n88_12(i12) & n88_14~=n88_13(i13));
																																		% 	sz14=size(n88_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n88_14(i14) == loc_of_sec_reg)
																																		% 			display('The shortest path to reg 6 is 14');
																																		% 			weightp14_reg6(t14)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n88_15=find(neighborss(n88_14(i14),:));
																												 					% 					n88_15=n88_15(n88_15~=88 & n88_15~=n88_2(i2) & n88_15~=n88_3(i3) & n88_15~=n88_4(i4) & n88_15~=n88_5(i5) & n88_15~=n88_6(i6) & n88_15~=n88_7(i7) & n88_15~=n88_8(i8) & n88_15~=n88_9(i9) & n88_15~=n88_10(i10) & n88_15~=n88_11(i11) & n88_15~=n88_12(i12) & n88_15~=n88_13(i13) & n88_15~=n88_14(i14));
																												 					% 					sz15=size(n88_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n88_15(i15) == loc_of_sec_reg)
																												 					% 							display('The shortest path to reg 6 is 15');
																												 					% 							weightp15_reg6(t15)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14))*lev(n88_14(i14),n88_15(i15));
																												 					% 							t15=t15+1;
																												 					% 						end
																																		% 				end
																																		% 			end
																																		% end
																																	end
																																end
																														end
																													end
																											end
																										end
																								end
																							end
																					end
																				end
																			end
																		end
																end
															end
                                    				end
                                				end

                                        	end
                                    	end

                                    end
                                end
                                
    						end
    					end
    				end
    			end
    		end
    	end
    
  
t1=1; t2=1; t3=1; t4=1; t5=1; t6=1; t7=1; t8=1; t9=1; t10=1; t11=1; t12=1;
t13=1; t14=1; t15=1;
loc_of_sec_reg = 64;

    	for i1=1:sz(2)
    		if (n88(i1) == loc_of_sec_reg)
    			display('The shortest path to reg 7 is 1');
    			weightp1_reg7(t1)=lev(88,n88(i1));
    			t1=t1+1;
    		else
    			n88_2=find(neighborss(n88(i1),:));
    			n88_2=n88_2(n88_2~=88);
    			sz2=size(n88_2);
    			for i2=1:sz2(2)
    				if (n88_2(i2) == loc_of_sec_reg)
    					display('The shortest path to reg 7 is 2');
    					weightp2_reg7(t2)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2));
    					t2=t2+1;
    				else
    					n88_3=find(neighborss(n88_2(i2),:));
    					n88_3=n88_3(n88_3~=n88_2(i2) & n88_3~=88);
    					sz3=size(n88_3);
    					for i3=1:sz3(2)
    						if (n88_3(i3) == loc_of_sec_reg)
    							display('The shortest path to reg 7 is 3');
    							weightp3_reg7(t3)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3));
    							t3=t3+1;
    						else
    							n88_4=find(neighborss(n88_3(i3),:));
    							n88_4=n88_4(n88_4~=n88_2(i2) & n88_4~=n88_3(i3) & n88_4~=88);
    							sz4=size(n88_4);
    							for i4=1:sz4(2)
    								if (n88_4(i4) == loc_of_sec_reg)
    								display('The shortest path to reg 7 is 4');
    								weightp4_reg7(t4)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4));
    								t4=t4+1;
        							else
    									n88_5=find(neighborss(n88_4(i4),:));
    									n88_5=n88_5(n88_5~=n88_2(i2) & n88_5~=n88_3(i3) & n88_5~=n88_4(i4) & n88_5~=88);
    									sz5=size(n88_5);
    									for i5=1:sz5(2)
    										if (n88_5(i5) == loc_of_sec_reg)
    										display('The shortest path to reg 7 is 5');
    										weightp5_reg7(t5)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5));
    										t5=t5+1;
    										else
    											n88_6=find(neighborss(n88_5(i5),:));
    											n88_6=n88_6(n88_6~=88 & n88_6~=n88_2(i2) & n88_6~=n88_3(i3) & n88_6~=n88_4(i4) & n88_6~=n88_5(i5));
    											sz6=size(n88_6);
    											for i6=1:sz6(2)
    												if (n88_6(i6) == loc_of_sec_reg)
    													display('The shortest path to reg 7 is 6');
    													weightp6_reg7(t6)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6));
    													t6=t6+1;
    													else
															n88_7=find(neighborss(n88_6(i6),:));
															n88_7=n88_7(n88_7~=88 & n88_7~=n88_2(i2) & n88_7~=n88_3(i3) & n88_7~=n88_4(i4) & n88_7~=n88_5(i5) & n88_7~=n88_6(i6));
															sz7=size(n88_7);
															for i7=1:sz7(2)
																if (n88_7(i7) == loc_of_sec_reg)
																	display('The shortest path to reg 7 is 7');
																	weightp7_reg7(t7)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7));
																	t7=t7+1;
																	else
																		n88_8=find(neighborss(n88_7(i7),:));
																		n88_8=n88_8(n88_8~=88 & n88_8~=n88_2(i2) & n88_8~=n88_3(i3) & n88_8~=n88_4(i4) & n88_8~=n88_5(i5) & n88_8~=n88_6(i6) & n88_8~=n88_7(i7));
																		sz8=size(n88_8);
																		for i8=1:sz8(2)
																			if (n88_8(i8) == loc_of_sec_reg)
																			display('The shortest path to reg 7 is 8');
																			weightp8_reg7(t8)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8));
																			t8=t8+1;
																			else
																				n88_9=find(neighborss(n88_8(i8),:));
																				n88_9=n88_9(n88_9~=88 & n88_9~=n88_2(i2) & n88_9~=n88_3(i3) & n88_9~=n88_4(i4) & n88_9~=n88_5(i5) & n88_9~=n88_6(i6) & n88_9~=n88_7(i7) & n88_9~=n88_8(i8));
																				sz9=size(n88_9);
																				for i9=1:sz9(2)
																					if (n88_9(i9) == loc_of_sec_reg)
																						display('The shortest path to reg 7 is 9');
																						weightp9_reg7(t9)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9));
																						t9=t9+1;
																						else
																							n88_10=find(neighborss(n88_9(i9),:));
																							n88_10=n88_10(n88_10~=88 & n88_10~=n88_2(i2) & n88_10~=n88_3(i3) & n88_10~=n88_4(i4) & n88_10~=n88_5(i5) & n88_10~=n88_6(i6) & n88_10~=n88_7(i7) & n88_10~=n88_8(i8) & n88_10~=n88_9(i9));
																							sz10=size(n88_10);
																							for i10=1:sz10(2)
																								if (n88_10(i10) == loc_of_sec_reg)
																									display('The shortest path to reg 7 is 10');
																									weightp10_reg7(t10)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10));
																									t10=t10+1;
																									else
																										n88_11=find(neighborss(n88_10(i10),:));
																										n88_11=n88_11(n88_11~=88 & n88_11~=n88_2(i2) & n88_11~=n88_3(i3) & n88_11~=n88_4(i4) & n88_11~=n88_5(i5) & n88_11~=n88_6(i6) & n88_11~=n88_7(i7) & n88_11~=n88_8(i8) & n88_11~=n88_9(i9) & n88_11~=n88_10(i10));
																										sz11=size(n88_11);
																										for i11=1:sz11(2)
																											if (n88_11(i11) == loc_of_sec_reg)
																												display('The shortest path to reg 7 is 11');
																												weightp11_reg7(t11)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11));
																												t11=t11+1;
																												else
																													n88_12=find(neighborss(n88_11(i11),:));
																													n88_12=n88_12(n88_12~=88 & n88_12~=n88_2(i2) & n88_12~=n88_3(i3) & n88_12~=n88_4(i4) & n88_12~=n88_5(i5) & n88_12~=n88_6(i6) & n88_12~=n88_7(i7) & n88_12~=n88_8(i8) & n88_12~=n88_9(i9) & n88_12~=n88_10(i10) & n88_12~=n88_11(i11));
																													sz12=size(n88_12);
																													for i12=1:sz12(2)
																														if (n88_12(i12) == loc_of_sec_reg)
																															display('The shortest path to reg 7 is 12');
																															weightp12_reg7(t12)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12));
																															t12=t12+1;
																															else
																																n88_13=find(neighborss(n88_12(i12),:));
																																n88_13=n88_13(n88_13~=88 & n88_13~=n88_2(i2) & n88_13~=n88_3(i3) & n88_13~=n88_4(i4) & n88_13~=n88_5(i5) & n88_13~=n88_6(i6) & n88_13~=n88_7(i7) & n88_13~=n88_8(i8) & n88_13~=n88_9(i9) & n88_13~=n88_10(i10) & n88_13~=n88_11(i11) & n88_13~=n88_12(i12));
																																sz13=size(n88_13);
																																for i13=1:sz13(2)
																																	if (n88_13(i13) == loc_of_sec_reg)
																																		display('The shortest path to reg 7 is 13');
																																		weightp13_reg7(t13)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n88_14=find(neighborss(n88_13(i13),:));
																																		% 	n88_14=n88_14(n88_14~=88 & n88_14~=n88_2(i2) & n88_14~=n88_3(i3) & n88_14~=n88_4(i4) & n88_14~=n88_5(i5) & n88_14~=n88_6(i6) & n88_14~=n88_7(i7) & n88_14~=n88_8(i8) & n88_14~=n88_9(i9) & n88_14~=n88_10(i10) & n88_14~=n88_11(i11) & n88_14~=n88_12(i12) & n88_14~=n88_13(i13));
																																		% 	sz14=size(n88_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n88_14(i14) == loc_of_sec_reg)
																																		% 			display('The shortest path to reg 7 is 14');
																																		% 			weightp14_reg7(t14)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n88_15=find(neighborss(n88_14(i14),:));
																												 					% 					n88_15=n88_15(n88_15~=88 & n88_15~=n88_2(i2) & n88_15~=n88_3(i3) & n88_15~=n88_4(i4) & n88_15~=n88_5(i5) & n88_15~=n88_6(i6) & n88_15~=n88_7(i7) & n88_15~=n88_8(i8) & n88_15~=n88_9(i9) & n88_15~=n88_10(i10) & n88_15~=n88_11(i11) & n88_15~=n88_12(i12) & n88_15~=n88_13(i13) & n88_15~=n88_14(i14));
																												 					% 					sz15=size(n88_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n88_15(i15) == loc_of_sec_reg)
																												 					% 							display('The shortest path to reg 7 is 15');
																												 					% 							weightp15_reg7(t15)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14))*lev(n88_14(i14),n88_15(i15));
																												 					% 							t15=t15+1;
																												 					% 						end
																																		% 				end
																																		% 			end
																																		% end
																																	end
																																end
																														end
																													end
																											end
																										end
																								end
																							end
																					end
																				end
																			end
																		end
																end
															end
                                    				end
                                				end

                                        	end
                                    	end

                                    end
                                end
                                
    						end
    					end
    				end
    			end
    		end
    	end
    
  
t1=1; t2=1; t3=1; t4=1; t5=1; t6=1; t7=1; t8=1; t9=1; t10=1; t11=1; t12=1;
t13=1; t14=1; t15=1;
loc_of_sec_reg = 31;

    	for i1=1:sz(2)
    		if (n88(i1) == loc_of_sec_reg)
    			display('The shortest path to reg 8 is 1');
    			weightp1_reg8(t1)=lev(88,n88(i1));
    			t1=t1+1;
    		else
    			n88_2=find(neighborss(n88(i1),:));
    			n88_2=n88_2(n88_2~=88);
    			sz2=size(n88_2);
    			for i2=1:sz2(2)
    				if (n88_2(i2) == loc_of_sec_reg)
    					display('The shortest path to reg 8 is 2');
    					weightp2_reg8(t2)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2));
    					t2=t2+1;
    				else
    					n88_3=find(neighborss(n88_2(i2),:));
    					n88_3=n88_3(n88_3~=n88_2(i2) & n88_3~=88);
    					sz3=size(n88_3);
    					for i3=1:sz3(2)
    						if (n88_3(i3) == loc_of_sec_reg)
    							display('The shortest path to reg 8 is 3');
    							weightp3_reg8(t3)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3));
    							t3=t3+1;
    						else
    							n88_4=find(neighborss(n88_3(i3),:));
    							n88_4=n88_4(n88_4~=n88_2(i2) & n88_4~=n88_3(i3) & n88_4~=88);
    							sz4=size(n88_4);
    							for i4=1:sz4(2)
    								if (n88_4(i4) == loc_of_sec_reg)
    								display('The shortest path to reg 8 is 4');
    								weightp4_reg8(t4)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4));
    								t4=t4+1;
        							else
    									n88_5=find(neighborss(n88_4(i4),:));
    									n88_5=n88_5(n88_5~=n88_2(i2) & n88_5~=n88_3(i3) & n88_5~=n88_4(i4) & n88_5~=88);
    									sz5=size(n88_5);
    									for i5=1:sz5(2)
    										if (n88_5(i5) == loc_of_sec_reg)
    										display('The shortest path to reg 8 is 5');
    										weightp5_reg8(t5)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5));
    										t5=t5+1;
    										else
    											n88_6=find(neighborss(n88_5(i5),:));
    											n88_6=n88_6(n88_6~=88 & n88_6~=n88_2(i2) & n88_6~=n88_3(i3) & n88_6~=n88_4(i4) & n88_6~=n88_5(i5));
    											sz6=size(n88_6);
    											for i6=1:sz6(2)
    												if (n88_6(i6) == loc_of_sec_reg)
    													display('The shortest path to reg 8 is 6');
    													weightp6_reg8(t6)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6));
    													t6=t6+1;
    													else
															n88_7=find(neighborss(n88_6(i6),:));
															n88_7=n88_7(n88_7~=88 & n88_7~=n88_2(i2) & n88_7~=n88_3(i3) & n88_7~=n88_4(i4) & n88_7~=n88_5(i5) & n88_7~=n88_6(i6));
															sz7=size(n88_7);
															for i7=1:sz7(2)
																if (n88_7(i7) == loc_of_sec_reg)
																	display('The shortest path to reg 8 is 7');
																	weightp7_reg8(t7)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7));
																	t7=t7+1;
																	else
																		n88_8=find(neighborss(n88_7(i7),:));
																		n88_8=n88_8(n88_8~=88 & n88_8~=n88_2(i2) & n88_8~=n88_3(i3) & n88_8~=n88_4(i4) & n88_8~=n88_5(i5) & n88_8~=n88_6(i6) & n88_8~=n88_7(i7));
																		sz8=size(n88_8);
																		for i8=1:sz8(2)
																			if (n88_8(i8) == loc_of_sec_reg)
																			display('The shortest path to reg 8 is 8');
																			weightp8_reg8(t8)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8));
																			t8=t8+1;
																			else
																				n88_9=find(neighborss(n88_8(i8),:));
																				n88_9=n88_9(n88_9~=88 & n88_9~=n88_2(i2) & n88_9~=n88_3(i3) & n88_9~=n88_4(i4) & n88_9~=n88_5(i5) & n88_9~=n88_6(i6) & n88_9~=n88_7(i7) & n88_9~=n88_8(i8));
																				sz9=size(n88_9);
																				for i9=1:sz9(2)
																					if (n88_9(i9) == loc_of_sec_reg)
																						display('The shortest path to reg 8 is 9');
																						weightp9_reg8(t9)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9));
																						t9=t9+1;
																						else
																							n88_10=find(neighborss(n88_9(i9),:));
																							n88_10=n88_10(n88_10~=88 & n88_10~=n88_2(i2) & n88_10~=n88_3(i3) & n88_10~=n88_4(i4) & n88_10~=n88_5(i5) & n88_10~=n88_6(i6) & n88_10~=n88_7(i7) & n88_10~=n88_8(i8) & n88_10~=n88_9(i9));
																							sz10=size(n88_10);
																							for i10=1:sz10(2)
																								if (n88_10(i10) == loc_of_sec_reg)
																									display('The shortest path to reg 8 is 10');
																									weightp10_reg8(t10)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10));
																									t10=t10+1;
																									else
																										n88_11=find(neighborss(n88_10(i10),:));
																										n88_11=n88_11(n88_11~=88 & n88_11~=n88_2(i2) & n88_11~=n88_3(i3) & n88_11~=n88_4(i4) & n88_11~=n88_5(i5) & n88_11~=n88_6(i6) & n88_11~=n88_7(i7) & n88_11~=n88_8(i8) & n88_11~=n88_9(i9) & n88_11~=n88_10(i10));
																										sz11=size(n88_11);
																										for i11=1:sz11(2)
																											if (n88_11(i11) == loc_of_sec_reg)
																												display('The shortest path to reg 8 is 11');
																												weightp11_reg8(t11)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11));
																												t11=t11+1;
																												else
																													n88_12=find(neighborss(n88_11(i11),:));
																													n88_12=n88_12(n88_12~=88 & n88_12~=n88_2(i2) & n88_12~=n88_3(i3) & n88_12~=n88_4(i4) & n88_12~=n88_5(i5) & n88_12~=n88_6(i6) & n88_12~=n88_7(i7) & n88_12~=n88_8(i8) & n88_12~=n88_9(i9) & n88_12~=n88_10(i10) & n88_12~=n88_11(i11));
																													sz12=size(n88_12);
																													for i12=1:sz12(2)
																														if (n88_12(i12) == loc_of_sec_reg)
																															display('The shortest path to reg 8 is 12');
																															weightp12_reg8(t12)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12));
																															t12=t12+1;
																															else
																																n88_13=find(neighborss(n88_12(i12),:));
																																n88_13=n88_13(n88_13~=88 & n88_13~=n88_2(i2) & n88_13~=n88_3(i3) & n88_13~=n88_4(i4) & n88_13~=n88_5(i5) & n88_13~=n88_6(i6) & n88_13~=n88_7(i7) & n88_13~=n88_8(i8) & n88_13~=n88_9(i9) & n88_13~=n88_10(i10) & n88_13~=n88_11(i11) & n88_13~=n88_12(i12));
																																sz13=size(n88_13);
																																for i13=1:sz13(2)
																																	if (n88_13(i13) == loc_of_sec_reg)
																																		display('The shortest path to reg 8 is 13');
																																		weightp13_reg8(t13)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n88_14=find(neighborss(n88_13(i13),:));
																																		% 	n88_14=n88_14(n88_14~=88 & n88_14~=n88_2(i2) & n88_14~=n88_3(i3) & n88_14~=n88_4(i4) & n88_14~=n88_5(i5) & n88_14~=n88_6(i6) & n88_14~=n88_7(i7) & n88_14~=n88_8(i8) & n88_14~=n88_9(i9) & n88_14~=n88_10(i10) & n88_14~=n88_11(i11) & n88_14~=n88_12(i12) & n88_14~=n88_13(i13));
																																		% 	sz14=size(n88_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n88_14(i14) == loc_of_sec_reg)
																																		% 			display('The shortest path to reg 8 is 14');
																																		% 			weightp14_reg8(t14)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n88_15=find(neighborss(n88_14(i14),:));
																												 					% 					n88_15=n88_15(n88_15~=88 & n88_15~=n88_2(i2) & n88_15~=n88_3(i3) & n88_15~=n88_4(i4) & n88_15~=n88_5(i5) & n88_15~=n88_6(i6) & n88_15~=n88_7(i7) & n88_15~=n88_8(i8) & n88_15~=n88_9(i9) & n88_15~=n88_10(i10) & n88_15~=n88_11(i11) & n88_15~=n88_12(i12) & n88_15~=n88_13(i13) & n88_15~=n88_14(i14));
																												 					% 					sz15=size(n88_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n88_15(i15) == loc_of_sec_reg)
																												 					% 							display('The shortest path to reg 8 is 15');
																												 					% 							weightp15_reg8(t15)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14))*lev(n88_14(i14),n88_15(i15));
																												 					% 							t15=t15+1;
																												 					% 						end
																																		% 				end
																																		% 			end
																																		% end
																																	end
																																end
																														end
																													end
																											end
																										end
																								end
																							end
																					end
																				end
																			end
																		end
																end
															end
                                    				end
                                				end

                                        	end
                                    	end

                                    end
                                end
                                
    						end
    					end
    				end
    			end
    		end
    	end
    
  
t1=1; t2=1; t3=1; t4=1; t5=1; t6=1; t7=1; t8=1; t9=1; t10=1; t11=1; t12=1;
t13=1; t14=1; t15=1;
loc_of_sec_reg = 80;

    	for i1=1:sz(2)
    		if (n88(i1) == loc_of_sec_reg)
    			display('The shortest path to reg 9 is 1');
    			weightp1_reg9(t1)=lev(88,n88(i1));
    			t1=t1+1;
    		else
    			n88_2=find(neighborss(n88(i1),:));
    			n88_2=n88_2(n88_2~=88);
    			sz2=size(n88_2);
    			for i2=1:sz2(2)
    				if (n88_2(i2) == loc_of_sec_reg)
    					display('The shortest path to reg 9 is 2');
    					weightp2_reg9(t2)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2));
    					t2=t2+1;
    				else
    					n88_3=find(neighborss(n88_2(i2),:));
    					n88_3=n88_3(n88_3~=n88_2(i2) & n88_3~=88);
    					sz3=size(n88_3);
    					for i3=1:sz3(2)
    						if (n88_3(i3) == loc_of_sec_reg)
    							display('The shortest path to reg 9 is 3');
    							weightp3_reg9(t3)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3));
    							t3=t3+1;
    						else
    							n88_4=find(neighborss(n88_3(i3),:));
    							n88_4=n88_4(n88_4~=n88_2(i2) & n88_4~=n88_3(i3) & n88_4~=88);
    							sz4=size(n88_4);
    							for i4=1:sz4(2)
    								if (n88_4(i4) == loc_of_sec_reg)
    								display('The shortest path to reg 9 is 4');
    								weightp4_reg9(t4)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4));
    								t4=t4+1;
        							else
    									n88_5=find(neighborss(n88_4(i4),:));
    									n88_5=n88_5(n88_5~=n88_2(i2) & n88_5~=n88_3(i3) & n88_5~=n88_4(i4) & n88_5~=88);
    									sz5=size(n88_5);
    									for i5=1:sz5(2)
    										if (n88_5(i5) == loc_of_sec_reg)
    										display('The shortest path to reg 9 is 5');
    										weightp5_reg9(t5)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5));
    										t5=t5+1;
    										else
    											n88_6=find(neighborss(n88_5(i5),:));
    											n88_6=n88_6(n88_6~=88 & n88_6~=n88_2(i2) & n88_6~=n88_3(i3) & n88_6~=n88_4(i4) & n88_6~=n88_5(i5));
    											sz6=size(n88_6);
    											for i6=1:sz6(2)
    												if (n88_6(i6) == loc_of_sec_reg)
    													display('The shortest path to reg 9 is 6');
    													weightp6_reg9(t6)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6));
    													t6=t6+1;
    													else
															n88_7=find(neighborss(n88_6(i6),:));
															n88_7=n88_7(n88_7~=88 & n88_7~=n88_2(i2) & n88_7~=n88_3(i3) & n88_7~=n88_4(i4) & n88_7~=n88_5(i5) & n88_7~=n88_6(i6));
															sz7=size(n88_7);
															for i7=1:sz7(2)
																if (n88_7(i7) == loc_of_sec_reg)
																	display('The shortest path to reg 9 is 7');
																	weightp7_reg9(t7)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7));
																	t7=t7+1;
																	else
																		n88_8=find(neighborss(n88_7(i7),:));
																		n88_8=n88_8(n88_8~=88 & n88_8~=n88_2(i2) & n88_8~=n88_3(i3) & n88_8~=n88_4(i4) & n88_8~=n88_5(i5) & n88_8~=n88_6(i6) & n88_8~=n88_7(i7));
																		sz8=size(n88_8);
																		for i8=1:sz8(2)
																			if (n88_8(i8) == loc_of_sec_reg)
																			display('The shortest path to reg 9 is 8');
																			weightp8_reg9(t8)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8));
																			t8=t8+1;
																			else
																				n88_9=find(neighborss(n88_8(i8),:));
																				n88_9=n88_9(n88_9~=88 & n88_9~=n88_2(i2) & n88_9~=n88_3(i3) & n88_9~=n88_4(i4) & n88_9~=n88_5(i5) & n88_9~=n88_6(i6) & n88_9~=n88_7(i7) & n88_9~=n88_8(i8));
																				sz9=size(n88_9);
																				for i9=1:sz9(2)
																					if (n88_9(i9) == loc_of_sec_reg)
																						display('The shortest path to reg 9 is 9');
																						weightp9_reg9(t9)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9));
																						t9=t9+1;
																						else
																							n88_10=find(neighborss(n88_9(i9),:));
																							n88_10=n88_10(n88_10~=88 & n88_10~=n88_2(i2) & n88_10~=n88_3(i3) & n88_10~=n88_4(i4) & n88_10~=n88_5(i5) & n88_10~=n88_6(i6) & n88_10~=n88_7(i7) & n88_10~=n88_8(i8) & n88_10~=n88_9(i9));
																							sz10=size(n88_10);
																							for i10=1:sz10(2)
																								if (n88_10(i10) == loc_of_sec_reg)
																									display('The shortest path to reg 9 is 10');
																									weightp10_reg9(t10)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10));
																									t10=t10+1;
																									else
																										n88_11=find(neighborss(n88_10(i10),:));
																										n88_11=n88_11(n88_11~=88 & n88_11~=n88_2(i2) & n88_11~=n88_3(i3) & n88_11~=n88_4(i4) & n88_11~=n88_5(i5) & n88_11~=n88_6(i6) & n88_11~=n88_7(i7) & n88_11~=n88_8(i8) & n88_11~=n88_9(i9) & n88_11~=n88_10(i10));
																										sz11=size(n88_11);
																										for i11=1:sz11(2)
																											if (n88_11(i11) == loc_of_sec_reg)
																												display('The shortest path to reg 9 is 11');
																												weightp11_reg9(t11)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11));
																												t11=t11+1;
																												else
																													n88_12=find(neighborss(n88_11(i11),:));
																													n88_12=n88_12(n88_12~=88 & n88_12~=n88_2(i2) & n88_12~=n88_3(i3) & n88_12~=n88_4(i4) & n88_12~=n88_5(i5) & n88_12~=n88_6(i6) & n88_12~=n88_7(i7) & n88_12~=n88_8(i8) & n88_12~=n88_9(i9) & n88_12~=n88_10(i10) & n88_12~=n88_11(i11));
																													sz12=size(n88_12);
																													for i12=1:sz12(2)
																														if (n88_12(i12) == loc_of_sec_reg)
																															display('The shortest path to reg 9 is 12');
																															weightp12_reg9(t12)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12));
																															t12=t12+1;
																															else
																																n88_13=find(neighborss(n88_12(i12),:));
																																n88_13=n88_13(n88_13~=88 & n88_13~=n88_2(i2) & n88_13~=n88_3(i3) & n88_13~=n88_4(i4) & n88_13~=n88_5(i5) & n88_13~=n88_6(i6) & n88_13~=n88_7(i7) & n88_13~=n88_8(i8) & n88_13~=n88_9(i9) & n88_13~=n88_10(i10) & n88_13~=n88_11(i11) & n88_13~=n88_12(i12));
																																sz13=size(n88_13);
																																for i13=1:sz13(2)
																																	if (n88_13(i13) == loc_of_sec_reg)
																																		display('The shortest path to reg 9 is 13');
																																		weightp13_reg9(t13)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n88_14=find(neighborss(n88_13(i13),:));
																																		% 	n88_14=n88_14(n88_14~=88 & n88_14~=n88_2(i2) & n88_14~=n88_3(i3) & n88_14~=n88_4(i4) & n88_14~=n88_5(i5) & n88_14~=n88_6(i6) & n88_14~=n88_7(i7) & n88_14~=n88_8(i8) & n88_14~=n88_9(i9) & n88_14~=n88_10(i10) & n88_14~=n88_11(i11) & n88_14~=n88_12(i12) & n88_14~=n88_13(i13));
																																		% 	sz14=size(n88_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n88_14(i14) == loc_of_sec_reg)
																																		% 			display('The shortest path to reg 9 is 14');
																																		% 			weightp14_reg9(t14)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n88_15=find(neighborss(n88_14(i14),:));
																												 					% 					n88_15=n88_15(n88_15~=88 & n88_15~=n88_2(i2) & n88_15~=n88_3(i3) & n88_15~=n88_4(i4) & n88_15~=n88_5(i5) & n88_15~=n88_6(i6) & n88_15~=n88_7(i7) & n88_15~=n88_8(i8) & n88_15~=n88_9(i9) & n88_15~=n88_10(i10) & n88_15~=n88_11(i11) & n88_15~=n88_12(i12) & n88_15~=n88_13(i13) & n88_15~=n88_14(i14));
																												 					% 					sz15=size(n88_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n88_15(i15) == loc_of_sec_reg)
																												 					% 							display('The shortest path to reg 9 is 15');
																												 					% 							weightp15_reg9(t15)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14))*lev(n88_14(i14),n88_15(i15));
																												 					% 							t15=t15+1;
																												 					% 						end
																																		% 				end
																																		% 			end
																																		% end
																																	end
																																end
																														end
																													end
																											end
																										end
																								end
																							end
																					end
																				end
																			end
																		end
																end
															end
                                    				end
                                				end

                                        	end
                                    	end

                                    end
                                end
                                
    						end
    					end
    				end
    			end
    		end
    	end
    
  
t1=1; t2=1; t3=1; t4=1; t5=1; t6=1; t7=1; t8=1; t9=1; t10=1; t11=1; t12=1;
t13=1; t14=1; t15=1;
loc_of_sec_reg = 10;

    	for i1=1:sz(2)
    		if (n88(i1) == loc_of_sec_reg)
    			display('The shortest path to reg 10 is 1');
    			weightp1_reg10(t1)=lev(88,n88(i1));
    			t1=t1+1;
    		else
    			n88_2=find(neighborss(n88(i1),:));
    			n88_2=n88_2(n88_2~=88);
    			sz2=size(n88_2);
    			for i2=1:sz2(2)
    				if (n88_2(i2) == loc_of_sec_reg)
    					display('The shortest path to reg 10 is 2');
    					weightp2_reg10(t2)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2));
    					t2=t2+1;
    				else
    					n88_3=find(neighborss(n88_2(i2),:));
    					n88_3=n88_3(n88_3~=n88_2(i2) & n88_3~=88);
    					sz3=size(n88_3);
    					for i3=1:sz3(2)
    						if (n88_3(i3) == loc_of_sec_reg)
    							display('The shortest path to reg 10 is 3');
    							weightp3_reg10(t3)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3));
    							t3=t3+1;
    						else
    							n88_4=find(neighborss(n88_3(i3),:));
    							n88_4=n88_4(n88_4~=n88_2(i2) & n88_4~=n88_3(i3) & n88_4~=88);
    							sz4=size(n88_4);
    							for i4=1:sz4(2)
    								if (n88_4(i4) == loc_of_sec_reg)
    								display('The shortest path to reg 10 is 4');
    								weightp4_reg10(t4)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4));
    								t4=t4+1;
        							else
    									n88_5=find(neighborss(n88_4(i4),:));
    									n88_5=n88_5(n88_5~=n88_2(i2) & n88_5~=n88_3(i3) & n88_5~=n88_4(i4) & n88_5~=88);
    									sz5=size(n88_5);
    									for i5=1:sz5(2)
    										if (n88_5(i5) == loc_of_sec_reg)
    										display('The shortest path to reg 10 is 5');
    										weightp5_reg10(t5)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5));
    										t5=t5+1;
    										else
    											n88_6=find(neighborss(n88_5(i5),:));
    											n88_6=n88_6(n88_6~=88 & n88_6~=n88_2(i2) & n88_6~=n88_3(i3) & n88_6~=n88_4(i4) & n88_6~=n88_5(i5));
    											sz6=size(n88_6);
    											for i6=1:sz6(2)
    												if (n88_6(i6) == loc_of_sec_reg)
    													display('The shortest path to reg 10 is 6');
    													weightp6_reg10(t6)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6));
    													t6=t6+1;
    													else
															n88_7=find(neighborss(n88_6(i6),:));
															n88_7=n88_7(n88_7~=88 & n88_7~=n88_2(i2) & n88_7~=n88_3(i3) & n88_7~=n88_4(i4) & n88_7~=n88_5(i5) & n88_7~=n88_6(i6));
															sz7=size(n88_7);
															for i7=1:sz7(2)
																if (n88_7(i7) == loc_of_sec_reg)
																	display('The shortest path to reg 10 is 7');
																	weightp7_reg10(t7)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7));
																	t7=t7+1;
																	else
																		n88_8=find(neighborss(n88_7(i7),:));
																		n88_8=n88_8(n88_8~=88 & n88_8~=n88_2(i2) & n88_8~=n88_3(i3) & n88_8~=n88_4(i4) & n88_8~=n88_5(i5) & n88_8~=n88_6(i6) & n88_8~=n88_7(i7));
																		sz8=size(n88_8);
																		for i8=1:sz8(2)
																			if (n88_8(i8) == loc_of_sec_reg)
																			display('The shortest path to reg 10 is 8');
																			weightp8_reg10(t8)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8));
																			t8=t8+1;
																			else
																				n88_9=find(neighborss(n88_8(i8),:));
																				n88_9=n88_9(n88_9~=88 & n88_9~=n88_2(i2) & n88_9~=n88_3(i3) & n88_9~=n88_4(i4) & n88_9~=n88_5(i5) & n88_9~=n88_6(i6) & n88_9~=n88_7(i7) & n88_9~=n88_8(i8));
																				sz9=size(n88_9);
																				for i9=1:sz9(2)
																					if (n88_9(i9) == loc_of_sec_reg)
																						display('The shortest path to reg 10 is 9');
																						weightp9_reg10(t9)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9));
																						t9=t9+1;
																						else
																							n88_10=find(neighborss(n88_9(i9),:));
																							n88_10=n88_10(n88_10~=88 & n88_10~=n88_2(i2) & n88_10~=n88_3(i3) & n88_10~=n88_4(i4) & n88_10~=n88_5(i5) & n88_10~=n88_6(i6) & n88_10~=n88_7(i7) & n88_10~=n88_8(i8) & n88_10~=n88_9(i9));
																							sz10=size(n88_10);
																							for i10=1:sz10(2)
																								if (n88_10(i10) == loc_of_sec_reg)
																									display('The shortest path to reg 10 is 10');
																									weightp10_reg10(t10)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10));
																									t10=t10+1;
																									else
																										n88_11=find(neighborss(n88_10(i10),:));
																										n88_11=n88_11(n88_11~=88 & n88_11~=n88_2(i2) & n88_11~=n88_3(i3) & n88_11~=n88_4(i4) & n88_11~=n88_5(i5) & n88_11~=n88_6(i6) & n88_11~=n88_7(i7) & n88_11~=n88_8(i8) & n88_11~=n88_9(i9) & n88_11~=n88_10(i10));
																										sz11=size(n88_11);
																										for i11=1:sz11(2)
																											if (n88_11(i11) == loc_of_sec_reg)
																												display('The shortest path to reg 10 is 11');
																												weightp11_reg10(t11)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11));
																												t11=t11+1;
																												else
																													n88_12=find(neighborss(n88_11(i11),:));
																													n88_12=n88_12(n88_12~=88 & n88_12~=n88_2(i2) & n88_12~=n88_3(i3) & n88_12~=n88_4(i4) & n88_12~=n88_5(i5) & n88_12~=n88_6(i6) & n88_12~=n88_7(i7) & n88_12~=n88_8(i8) & n88_12~=n88_9(i9) & n88_12~=n88_10(i10) & n88_12~=n88_11(i11));
																													sz12=size(n88_12);
																													for i12=1:sz12(2)
																														if (n88_12(i12) == loc_of_sec_reg)
																															display('The shortest path to reg 10 is 12');
																															weightp12_reg10(t12)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12));
																															t12=t12+1;
																															else
																																n88_13=find(neighborss(n88_12(i12),:));
																																n88_13=n88_13(n88_13~=88 & n88_13~=n88_2(i2) & n88_13~=n88_3(i3) & n88_13~=n88_4(i4) & n88_13~=n88_5(i5) & n88_13~=n88_6(i6) & n88_13~=n88_7(i7) & n88_13~=n88_8(i8) & n88_13~=n88_9(i9) & n88_13~=n88_10(i10) & n88_13~=n88_11(i11) & n88_13~=n88_12(i12));
																																sz13=size(n88_13);
																																for i13=1:sz13(2)
																																	if (n88_13(i13) == loc_of_sec_reg)
																																		display('The shortest path to reg 10 is 13');
																																		weightp13_reg10(t13)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n88_14=find(neighborss(n88_13(i13),:));
																																		% 	n88_14=n88_14(n88_14~=88 & n88_14~=n88_2(i2) & n88_14~=n88_3(i3) & n88_14~=n88_4(i4) & n88_14~=n88_5(i5) & n88_14~=n88_6(i6) & n88_14~=n88_7(i7) & n88_14~=n88_8(i8) & n88_14~=n88_9(i9) & n88_14~=n88_10(i10) & n88_14~=n88_11(i11) & n88_14~=n88_12(i12) & n88_14~=n88_13(i13));
																																		% 	sz14=size(n88_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n88_14(i14) == loc_of_sec_reg)
																																		% 			display('The shortest path to reg 10 is 14');
																																		% 			weightp14_reg10(t14)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n88_15=find(neighborss(n88_14(i14),:));
																												 					% 					n88_15=n88_15(n88_15~=88 & n88_15~=n88_2(i2) & n88_15~=n88_3(i3) & n88_15~=n88_4(i4) & n88_15~=n88_5(i5) & n88_15~=n88_6(i6) & n88_15~=n88_7(i7) & n88_15~=n88_8(i8) & n88_15~=n88_9(i9) & n88_15~=n88_10(i10) & n88_15~=n88_11(i11) & n88_15~=n88_12(i12) & n88_15~=n88_13(i13) & n88_15~=n88_14(i14));
																												 					% 					sz15=size(n88_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n88_15(i15) == loc_of_sec_reg)
																												 					% 							display('The shortest path to reg 10 is 15');
																												 					% 							weightp15_reg10(t15)=lev(88,n88(i1))*lev(n88(i1),n88_2(i2))*lev(n88_2(i2),n88_3(i3))*lev(n88_3(i3),n88_4(i4))*lev(n88_4(i4),n88_5(i5))*lev(n88_5(i5),n88_6(i6))*lev(n88_6(i6),n88_7(i7))*lev(n88_7(i7),n88_8(i8))*lev(n88_8(i8),n88_9(i9))*lev(n88_9(i9),n88_10(i10))*lev(n88_10(i10),n88_11(i11))*lev(n88_11(i11),n88_12(i12))*lev(n88_12(i12),n88_13(i13))*lev(n88_13(i13),n88_14(i14))*lev(n88_14(i14),n88_15(i15));
																												 					% 							t15=t15+1;
																												 					% 						end
																																		% 				end
																																		% 			end
																																		% end
																																	end
																																end
																														end
																													end
																											end
																										end
																								end
																							end
																					end
																				end
																			end
																		end
																end
															end
                                    				end
                                				end

                                        	end
                                    	end

                                    end
                                end
                                
    						end
    					end
    				end
    			end
    		end
    	end
    
  

   
    
  

   