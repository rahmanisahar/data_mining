%%% Finding probabilities

net10by10=load('~/Desktop/net_subset10.mat');
%net10by10=load('~/Desktop/project/data_mining/nearby_galaxies/SOM/grey_cycle/2d/subsets/net_subset10.mat'); %Load network
net=net10by10.net;


weights = net.IW{1,1};
numInputs = net.inputs{1}.processedSize;
numNeurons = net.layers{1}.size;
topologyFcn = net.layers{1}.topologyFcn;


neighbors = sparse(tril(net.layers{1}.distances <= 1.001) - eye(numNeurons));
% neighbors2 = sparse(tril(net.layers{1}.distances >= 1.001 & net.layers{1}.distances <= 2.001 ));
% neighbors3 = sparse(tril(net.layers{1}.distances >= 2.001 & net.layers{1}.distances <= 3.001 ));
% neighbors4 = sparse(tril(net.layers{1}.distances >= 3.001 & net.layers{1}.distances <= 4.001 ));
% neighbors5 = sparse(tril(net.layers{1}.distances >= 4.001 & net.layers{1}.distances <= 5.001 ));
% neighbors6 = sparse(tril(net.layers{1}.distances >= 5.001 & net.layers{1}.distances <= 6.001 ));
% neighbors7 = sparse(tril(net.layers{1}.distances >= 6.001 & net.layers{1}.distances <= 7.001 ));
% neighbors8 = sparse(tril(net.layers{1}.distances >= 7.001 & net.layers{1}.distances <= 8.001 ));
% neighbors9 = sparse(tril(net.layers{1}.distances >= 8.001 & net.layers{1}.distances <= 9.001 ));
% neighbors10 = sparse(tril(net.layers{1}.distances >= 9.001 & net.layers{1}.distances <= 10.001 ));
% neighbors11 = sparse(tril(net.layers{1}.distances >= 10.001 & net.layers{1}.distances <= 11.001 ));
% neighbors12 = sparse(tril(net.layers{1}.distances >= 11.001 & net.layers{1}.distances <= 12.001 ));
% neighbors13 = sparse(tril(net.layers{1}.distances >= 12.001 & net.layers{1}.distances <= 13.001 ));
% neighbors14 = sparse(tril(net.layers{1}.distances >= 13.001 & net.layers{1}.distances <= 14.001 ));
% neighbors15 = sparse(tril(net.layers{1}.distances >= 14.001 & net.layers{1}.distances <= 15.001 ));

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
%lev=lev/100.;

loc_of_sec_reg_arr=[42 85 91 98 59 82 90 36 10];
loc_of_reg1=1;
    n_reg1=find(neighborss(loc_of_reg1,:));
    sz=size(n_reg1);

    %%PAth to reg 2.
t1=1; t2=1; t3=1; t4=1; t5=1; t6=1; t7=1; t8=1; t9=1; t10=1; t11=1; t12=1;
t13=1; t14=1; t15=1;
loc_of_sec_reg = loc_of_sec_reg_arr(1);
    	for i1=1:sz(2)
    		if (n_reg1(i1) == loc_of_sec_reg)
    			display('A path to reg 2 is 1');
    			weightp1_reg2(t1)=lev(loc_of_reg1,n_reg1(i1));
    			t1=t1+1;
    		else
    			n_reg1_2=find(neighborss(n_reg1(i1),:));
    			n_reg1_2=n_reg1_2(n_reg1_2~=loc_of_reg1);
    			sz2=size(n_reg1_2);
    			for i2=1:sz2(2)
    				if (n_reg1_2(i2) == loc_of_sec_reg)
    					display('A path to reg 2 is 2');
    					weightp2_reg2(t2)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2));
    					t2=t2+1;
    				else
    					n_reg1_3=find(neighborss(n_reg1_2(i2),:));
    					n_reg1_3=n_reg1_3(n_reg1_3~=n_reg1_2(i2) & n_reg1_3~=loc_of_reg1);
    					sz3=size(n_reg1_3);
    					for i3=1:sz3(2)
    						if (n_reg1_3(i3) == loc_of_sec_reg)
    							display('A path to reg 2 is 3');
    							weightp3_reg2(t3)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3));
    							t3=t3+1;
    						else
    							n_reg1_4=find(neighborss(n_reg1_3(i3),:));
    							n_reg1_4=n_reg1_4(n_reg1_4~=n_reg1_2(i2) & n_reg1_4~=n_reg1_3(i3) & n_reg1_4~=loc_of_reg1);
    							sz4=size(n_reg1_4);
    							for i4=1:sz4(2)
    								if (n_reg1_4(i4) == loc_of_sec_reg)
    								display('A path to reg 2 is 4');
    								weightp4_reg2(t4)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4));
    								t4=t4+1;
        							else
    									n_reg1_5=find(neighborss(n_reg1_4(i4),:));
    									n_reg1_5=n_reg1_5(n_reg1_5~=n_reg1_2(i2) & n_reg1_5~=n_reg1_3(i3) & n_reg1_5~=n_reg1_4(i4) & n_reg1_5~=loc_of_reg1);
    									sz5=size(n_reg1_5);
    									for i5=1:sz5(2)
    										if (n_reg1_5(i5) == loc_of_sec_reg)
    										display('A path to reg 2 is 5');
    										weightp5_reg2(t5)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5));
    										t5=t5+1;
    										else
    											n_reg1_6=find(neighborss(n_reg1_5(i5),:));
    											n_reg1_6=n_reg1_6(n_reg1_6~=loc_of_reg1 & n_reg1_6~=n_reg1_2(i2) & n_reg1_6~=n_reg1_3(i3) & n_reg1_6~=n_reg1_4(i4) & n_reg1_6~=n_reg1_5(i5));
    											sz6=size(n_reg1_6);
    											for i6=1:sz6(2)
    												if (n_reg1_6(i6) == loc_of_sec_reg)
    													display('A path to reg 2 is 6');
    													weightp6_reg2(t6)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6));
    													t6=t6+1;
    													else
															n_reg1_7=find(neighborss(n_reg1_6(i6),:));
															n_reg1_7=n_reg1_7(n_reg1_7~=loc_of_reg1 & n_reg1_7~=n_reg1_2(i2) & n_reg1_7~=n_reg1_3(i3) & n_reg1_7~=n_reg1_4(i4) & n_reg1_7~=n_reg1_5(i5) & n_reg1_7~=n_reg1_6(i6));
															sz7=size(n_reg1_7);
															for i7=1:sz7(2)
																if (n_reg1_7(i7) == loc_of_sec_reg)
																	display('A path to reg 2 is 7');
																	weightp7_reg2(t7)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7));
																	t7=t7+1;
																	else
																		n_reg1_8=find(neighborss(n_reg1_7(i7),:));
																		n_reg1_8=n_reg1_8(n_reg1_8~=loc_of_reg1 & n_reg1_8~=n_reg1_2(i2) & n_reg1_8~=n_reg1_3(i3) & n_reg1_8~=n_reg1_4(i4) & n_reg1_8~=n_reg1_5(i5) & n_reg1_8~=n_reg1_6(i6) & n_reg1_8~=n_reg1_7(i7));
																		sz8=size(n_reg1_8);
																		for i8=1:sz8(2)
																			if (n_reg1_8(i8) == loc_of_sec_reg)
																			display('A path to reg 2 is 8');
																			weightp8_reg2(t8)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8));
																			t8=t8+1;
																			else
																				n_reg1_9=find(neighborss(n_reg1_8(i8),:));
																				n_reg1_9=n_reg1_9(n_reg1_9~=loc_of_reg1 & n_reg1_9~=n_reg1_2(i2) & n_reg1_9~=n_reg1_3(i3) & n_reg1_9~=n_reg1_4(i4) & n_reg1_9~=n_reg1_5(i5) & n_reg1_9~=n_reg1_6(i6) & n_reg1_9~=n_reg1_7(i7) & n_reg1_9~=n_reg1_8(i8));
																				sz9=size(n_reg1_9);
																				for i9=1:sz9(2)
																					if (n_reg1_9(i9) == loc_of_sec_reg)
																						display('A path to reg 2 is 9');
																						weightp9_reg2(t9)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9));
																						t9=t9+1;
																						else
																							n_reg1_10=find(neighborss(n_reg1_9(i9),:));
																							n_reg1_10=n_reg1_10(n_reg1_10~=loc_of_reg1 & n_reg1_10~=n_reg1_2(i2) & n_reg1_10~=n_reg1_3(i3) & n_reg1_10~=n_reg1_4(i4) & n_reg1_10~=n_reg1_5(i5) & n_reg1_10~=n_reg1_6(i6) & n_reg1_10~=n_reg1_7(i7) & n_reg1_10~=n_reg1_8(i8) & n_reg1_10~=n_reg1_9(i9));
																							sz10=size(n_reg1_10);
																							for i10=1:sz10(2)
																								if (n_reg1_10(i10) == loc_of_sec_reg)
																									display('A path to reg 2 is 10');
																									weightp10_reg2(t10)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10));
																									t10=t10+1;
																									else
																										n_reg1_11=find(neighborss(n_reg1_10(i10),:));
																										n_reg1_11=n_reg1_11(n_reg1_11~=loc_of_reg1 & n_reg1_11~=n_reg1_2(i2) & n_reg1_11~=n_reg1_3(i3) & n_reg1_11~=n_reg1_4(i4) & n_reg1_11~=n_reg1_5(i5) & n_reg1_11~=n_reg1_6(i6) & n_reg1_11~=n_reg1_7(i7) & n_reg1_11~=n_reg1_8(i8) & n_reg1_11~=n_reg1_9(i9) & n_reg1_11~=n_reg1_10(i10));
																										sz11=size(n_reg1_11);
																										for i11=1:sz11(2)
																											if (n_reg1_11(i11) == loc_of_sec_reg)
																												display('A path to reg 2 is 11');
																												weightp11_reg2(t11)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11));
																												t11=t11+1;
																												else
																													n_reg1_12=find(neighborss(n_reg1_11(i11),:));
																													n_reg1_12=n_reg1_12(n_reg1_12~=loc_of_reg1 & n_reg1_12~=n_reg1_2(i2) & n_reg1_12~=n_reg1_3(i3) & n_reg1_12~=n_reg1_4(i4) & n_reg1_12~=n_reg1_5(i5) & n_reg1_12~=n_reg1_6(i6) & n_reg1_12~=n_reg1_7(i7) & n_reg1_12~=n_reg1_8(i8) & n_reg1_12~=n_reg1_9(i9) & n_reg1_12~=n_reg1_10(i10) & n_reg1_12~=n_reg1_11(i11));
																													sz12=size(n_reg1_12);
																													for i12=1:sz12(2)
																														if (n_reg1_12(i12) == loc_of_sec_reg)
																															display('A path to reg 2 is 12');
																															weightp12_reg2(t12)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12));
																															t12=t12+1;
																															else
																																n_reg1_13=find(neighborss(n_reg1_12(i12),:));
																																n_reg1_13=n_reg1_13(n_reg1_13~=loc_of_reg1 & n_reg1_13~=n_reg1_2(i2) & n_reg1_13~=n_reg1_3(i3) & n_reg1_13~=n_reg1_4(i4) & n_reg1_13~=n_reg1_5(i5) & n_reg1_13~=n_reg1_6(i6) & n_reg1_13~=n_reg1_7(i7) & n_reg1_13~=n_reg1_8(i8) & n_reg1_13~=n_reg1_9(i9) & n_reg1_13~=n_reg1_10(i10) & n_reg1_13~=n_reg1_11(i11) & n_reg1_13~=n_reg1_12(i12));
																																sz13=size(n_reg1_13);
																																for i13=1:sz13(2)
																																	if (n_reg1_13(i13) == loc_of_sec_reg)
																																		display('A path to reg 2 is 13');
																																		weightp13_reg2(t13)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n_reg1_14=find(neighborss(n_reg1_13(i13),:));
																																		% 	n_reg1_14=n_reg1_14(n_reg1_14~=loc_of_reg1 & n_reg1_14~=n_reg1_2(i2) & n_reg1_14~=n_reg1_3(i3) & n_reg1_14~=n_reg1_4(i4) & n_reg1_14~=n_reg1_5(i5) & n_reg1_14~=n_reg1_6(i6) & n_reg1_14~=n_reg1_7(i7) & n_reg1_14~=n_reg1_8(i8) & n_reg1_14~=n_reg1_9(i9) & n_reg1_14~=n_reg1_10(i10) & n_reg1_14~=n_reg1_11(i11) & n_reg1_14~=n_reg1_12(i12) & n_reg1_14~=n_reg1_13(i13));
																																		% 	sz14=size(n_reg1_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n_reg1_14(i14) == loc_of_sec_reg)
																																		% 			display('A path to reg 2 is 14');
																																		% 			weightp14_reg2(t14)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n_reg1_15=find(neighborss(n_reg1_14(i14),:));
																												 					% 					n_reg1_15=n_reg1_15(n_reg1_15~=loc_of_reg1 & n_reg1_15~=n_reg1_2(i2) & n_reg1_15~=n_reg1_3(i3) & n_reg1_15~=n_reg1_4(i4) & n_reg1_15~=n_reg1_5(i5) & n_reg1_15~=n_reg1_6(i6) & n_reg1_15~=n_reg1_7(i7) & n_reg1_15~=n_reg1_8(i8) & n_reg1_15~=n_reg1_9(i9) & n_reg1_15~=n_reg1_10(i10) & n_reg1_15~=n_reg1_11(i11) & n_reg1_15~=n_reg1_12(i12) & n_reg1_15~=n_reg1_13(i13) & n_reg1_15~=n_reg1_14(i14));
																												 					% 					sz15=size(n_reg1_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n_reg1_15(i15) == loc_of_sec_reg)
																												 					% 							display('A path to reg 2 is 15');
																												 					% 							weightp15_reg2(t15)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14))*lev(n_reg1_14(i14),n_reg1_15(i15));
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
loc_of_sec_reg = loc_of_sec_reg_arr(2);

    	for i1=1:sz(2)
    		if (n_reg1(i1) == loc_of_sec_reg)
    			display('A path to reg 3 is 1');
    			weightp1_reg3(t1)=lev(loc_of_reg1,n_reg1(i1));
    			t1=t1+1;
    		else
    			n_reg1_2=find(neighborss(n_reg1(i1),:));
    			n_reg1_2=n_reg1_2(n_reg1_2~=loc_of_reg1);
    			sz2=size(n_reg1_2);
    			for i2=1:sz2(2)
    				if (n_reg1_2(i2) == loc_of_sec_reg)
    					display('A path to reg 3 is 2');
    					weightp2_reg3(t2)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2));
    					t2=t2+1;
    				else
    					n_reg1_3=find(neighborss(n_reg1_2(i2),:));
    					n_reg1_3=n_reg1_3(n_reg1_3~=n_reg1_2(i2) & n_reg1_3~=loc_of_reg1);
    					sz3=size(n_reg1_3);
    					for i3=1:sz3(2)
    						if (n_reg1_3(i3) == loc_of_sec_reg)
    							display('A path to reg 3 is 3');
    							weightp3_reg3(t3)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3));
    							t3=t3+1;
    						else
    							n_reg1_4=find(neighborss(n_reg1_3(i3),:));
    							n_reg1_4=n_reg1_4(n_reg1_4~=n_reg1_2(i2) & n_reg1_4~=n_reg1_3(i3) & n_reg1_4~=loc_of_reg1);
    							sz4=size(n_reg1_4);
    							for i4=1:sz4(2)
    								if (n_reg1_4(i4) == loc_of_sec_reg)
    								display('A path to reg 3 is 4');
    								weightp4_reg3(t4)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4));
    								t4=t4+1;
        							else
    									n_reg1_5=find(neighborss(n_reg1_4(i4),:));
    									n_reg1_5=n_reg1_5(n_reg1_5~=n_reg1_2(i2) & n_reg1_5~=n_reg1_3(i3) & n_reg1_5~=n_reg1_4(i4) & n_reg1_5~=loc_of_reg1);
    									sz5=size(n_reg1_5);
    									for i5=1:sz5(2)
    										if (n_reg1_5(i5) == loc_of_sec_reg)
    										display('A path to reg 3 is 5');
    										weightp5_reg3(t5)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5));
    										t5=t5+1;
    										else
    											n_reg1_6=find(neighborss(n_reg1_5(i5),:));
    											n_reg1_6=n_reg1_6(n_reg1_6~=loc_of_reg1 & n_reg1_6~=n_reg1_2(i2) & n_reg1_6~=n_reg1_3(i3) & n_reg1_6~=n_reg1_4(i4) & n_reg1_6~=n_reg1_5(i5));
    											sz6=size(n_reg1_6);
    											for i6=1:sz6(2)
    												if (n_reg1_6(i6) == loc_of_sec_reg)
    													display('A path to reg 3 is 6');
    													weightp6_reg3(t6)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6));
    													t6=t6+1;
    													else
															n_reg1_7=find(neighborss(n_reg1_6(i6),:));
															n_reg1_7=n_reg1_7(n_reg1_7~=loc_of_reg1 & n_reg1_7~=n_reg1_2(i2) & n_reg1_7~=n_reg1_3(i3) & n_reg1_7~=n_reg1_4(i4) & n_reg1_7~=n_reg1_5(i5) & n_reg1_7~=n_reg1_6(i6));
															sz7=size(n_reg1_7);
															for i7=1:sz7(2)
																if (n_reg1_7(i7) == loc_of_sec_reg)
																	display('A path to reg 3 is 7');
																	weightp7_reg3(t7)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7));
																	t7=t7+1;
																	else
																		n_reg1_8=find(neighborss(n_reg1_7(i7),:));
																		n_reg1_8=n_reg1_8(n_reg1_8~=loc_of_reg1 & n_reg1_8~=n_reg1_2(i2) & n_reg1_8~=n_reg1_3(i3) & n_reg1_8~=n_reg1_4(i4) & n_reg1_8~=n_reg1_5(i5) & n_reg1_8~=n_reg1_6(i6) & n_reg1_8~=n_reg1_7(i7));
																		sz8=size(n_reg1_8);
																		for i8=1:sz8(2)
																			if (n_reg1_8(i8) == loc_of_sec_reg)
																			display('A path to reg 3 is 8');
																			weightp8_reg3(t8)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8));
																			t8=t8+1;
																			else
																				n_reg1_9=find(neighborss(n_reg1_8(i8),:));
																				n_reg1_9=n_reg1_9(n_reg1_9~=loc_of_reg1 & n_reg1_9~=n_reg1_2(i2) & n_reg1_9~=n_reg1_3(i3) & n_reg1_9~=n_reg1_4(i4) & n_reg1_9~=n_reg1_5(i5) & n_reg1_9~=n_reg1_6(i6) & n_reg1_9~=n_reg1_7(i7) & n_reg1_9~=n_reg1_8(i8));
																				sz9=size(n_reg1_9);
																				for i9=1:sz9(2)
																					if (n_reg1_9(i9) == loc_of_sec_reg)
																						display('A path to reg 3 is 9');
																						weightp9_reg3(t9)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9));
																						t9=t9+1;
																						else
																							n_reg1_10=find(neighborss(n_reg1_9(i9),:));
																							n_reg1_10=n_reg1_10(n_reg1_10~=loc_of_reg1 & n_reg1_10~=n_reg1_2(i2) & n_reg1_10~=n_reg1_3(i3) & n_reg1_10~=n_reg1_4(i4) & n_reg1_10~=n_reg1_5(i5) & n_reg1_10~=n_reg1_6(i6) & n_reg1_10~=n_reg1_7(i7) & n_reg1_10~=n_reg1_8(i8) & n_reg1_10~=n_reg1_9(i9));
																							sz10=size(n_reg1_10);
																							for i10=1:sz10(2)
																								if (n_reg1_10(i10) == loc_of_sec_reg)
																									display('A path to reg 3 is 10');
																									weightp10_reg3(t10)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10));
																									t10=t10+1;
																									else
																										n_reg1_11=find(neighborss(n_reg1_10(i10),:));
																										n_reg1_11=n_reg1_11(n_reg1_11~=loc_of_reg1 & n_reg1_11~=n_reg1_2(i2) & n_reg1_11~=n_reg1_3(i3) & n_reg1_11~=n_reg1_4(i4) & n_reg1_11~=n_reg1_5(i5) & n_reg1_11~=n_reg1_6(i6) & n_reg1_11~=n_reg1_7(i7) & n_reg1_11~=n_reg1_8(i8) & n_reg1_11~=n_reg1_9(i9) & n_reg1_11~=n_reg1_10(i10));
																										sz11=size(n_reg1_11);
																										for i11=1:sz11(2)
																											if (n_reg1_11(i11) == loc_of_sec_reg)
																												display('A path to reg 3 is 11');
																												weightp11_reg3(t11)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11));
																												t11=t11+1;
																												else
																													n_reg1_12=find(neighborss(n_reg1_11(i11),:));
																													n_reg1_12=n_reg1_12(n_reg1_12~=loc_of_reg1 & n_reg1_12~=n_reg1_2(i2) & n_reg1_12~=n_reg1_3(i3) & n_reg1_12~=n_reg1_4(i4) & n_reg1_12~=n_reg1_5(i5) & n_reg1_12~=n_reg1_6(i6) & n_reg1_12~=n_reg1_7(i7) & n_reg1_12~=n_reg1_8(i8) & n_reg1_12~=n_reg1_9(i9) & n_reg1_12~=n_reg1_10(i10) & n_reg1_12~=n_reg1_11(i11));
																													sz12=size(n_reg1_12);
																													for i12=1:sz12(2)
																														if (n_reg1_12(i12) == loc_of_sec_reg)
																															display('A path to reg 3 is 12');
																															weightp12_reg3(t12)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12));
																															t12=t12+1;
																															else
																																n_reg1_13=find(neighborss(n_reg1_12(i12),:));
																																n_reg1_13=n_reg1_13(n_reg1_13~=loc_of_reg1 & n_reg1_13~=n_reg1_2(i2) & n_reg1_13~=n_reg1_3(i3) & n_reg1_13~=n_reg1_4(i4) & n_reg1_13~=n_reg1_5(i5) & n_reg1_13~=n_reg1_6(i6) & n_reg1_13~=n_reg1_7(i7) & n_reg1_13~=n_reg1_8(i8) & n_reg1_13~=n_reg1_9(i9) & n_reg1_13~=n_reg1_10(i10) & n_reg1_13~=n_reg1_11(i11) & n_reg1_13~=n_reg1_12(i12));
																																sz13=size(n_reg1_13);
																																for i13=1:sz13(2)
																																	if (n_reg1_13(i13) == loc_of_sec_reg)
																																		display('A path to reg 3 is 13');
																																		weightp13_reg3(t13)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n_reg1_14=find(neighborss(n_reg1_13(i13),:));
																																		% 	n_reg1_14=n_reg1_14(n_reg1_14~=loc_of_reg1 & n_reg1_14~=n_reg1_2(i2) & n_reg1_14~=n_reg1_3(i3) & n_reg1_14~=n_reg1_4(i4) & n_reg1_14~=n_reg1_5(i5) & n_reg1_14~=n_reg1_6(i6) & n_reg1_14~=n_reg1_7(i7) & n_reg1_14~=n_reg1_8(i8) & n_reg1_14~=n_reg1_9(i9) & n_reg1_14~=n_reg1_10(i10) & n_reg1_14~=n_reg1_11(i11) & n_reg1_14~=n_reg1_12(i12) & n_reg1_14~=n_reg1_13(i13));
																																		% 	sz14=size(n_reg1_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n_reg1_14(i14) == loc_of_sec_reg)
																																		% 			display('A path to reg 3 is 14');
																																		% 			weightp14_reg3(t14)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n_reg1_15=find(neighborss(n_reg1_14(i14),:));
																												 					% 					n_reg1_15=n_reg1_15(n_reg1_15~=loc_of_reg1 & n_reg1_15~=n_reg1_2(i2) & n_reg1_15~=n_reg1_3(i3) & n_reg1_15~=n_reg1_4(i4) & n_reg1_15~=n_reg1_5(i5) & n_reg1_15~=n_reg1_6(i6) & n_reg1_15~=n_reg1_7(i7) & n_reg1_15~=n_reg1_8(i8) & n_reg1_15~=n_reg1_9(i9) & n_reg1_15~=n_reg1_10(i10) & n_reg1_15~=n_reg1_11(i11) & n_reg1_15~=n_reg1_12(i12) & n_reg1_15~=n_reg1_13(i13) & n_reg1_15~=n_reg1_14(i14));
																												 					% 					sz15=size(n_reg1_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n_reg1_15(i15) == loc_of_sec_reg)
																												 					% 							display('A path to reg 3 is 15');
																												 					% 							weightp15_reg3(t15)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14))*lev(n_reg1_14(i14),n_reg1_15(i15));
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
loc_of_sec_reg = loc_of_sec_reg_arr(3);

    	for i1=1:sz(2)
    		if (n_reg1(i1) == loc_of_sec_reg)
    			display('A path to reg 4 is 1');
    			weightp1_reg4(t1)=lev(loc_of_reg1,n_reg1(i1));
    			t1=t1+1;
    		else
    			n_reg1_2=find(neighborss(n_reg1(i1),:));
    			n_reg1_2=n_reg1_2(n_reg1_2~=loc_of_reg1);
    			sz2=size(n_reg1_2);
    			for i2=1:sz2(2)
    				if (n_reg1_2(i2) == loc_of_sec_reg)
    					display('A path to reg 4 is 2');
    					weightp2_reg4(t2)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2));
    					t2=t2+1;
    				else
    					n_reg1_3=find(neighborss(n_reg1_2(i2),:));
    					n_reg1_3=n_reg1_3(n_reg1_3~=n_reg1_2(i2) & n_reg1_3~=loc_of_reg1);
    					sz3=size(n_reg1_3);
    					for i3=1:sz3(2)
    						if (n_reg1_3(i3) == loc_of_sec_reg)
    							display('A path to reg 4 is 3');
    							weightp3_reg4(t3)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3));
    							t3=t3+1;
    						else
    							n_reg1_4=find(neighborss(n_reg1_3(i3),:));
    							n_reg1_4=n_reg1_4(n_reg1_4~=n_reg1_2(i2) & n_reg1_4~=n_reg1_3(i3) & n_reg1_4~=loc_of_reg1);
    							sz4=size(n_reg1_4);
    							for i4=1:sz4(2)
    								if (n_reg1_4(i4) == loc_of_sec_reg)
    								display('A path to reg 4 is 4');
    								weightp4_reg4(t4)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4));
    								t4=t4+1;
        							else
    									n_reg1_5=find(neighborss(n_reg1_4(i4),:));
    									n_reg1_5=n_reg1_5(n_reg1_5~=n_reg1_2(i2) & n_reg1_5~=n_reg1_3(i3) & n_reg1_5~=n_reg1_4(i4) & n_reg1_5~=loc_of_reg1);
    									sz5=size(n_reg1_5);
    									for i5=1:sz5(2)
    										if (n_reg1_5(i5) == loc_of_sec_reg)
    										display('A path to reg 4 is 5');
    										weightp5_reg4(t5)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5));
    										t5=t5+1;
    										else
    											n_reg1_6=find(neighborss(n_reg1_5(i5),:));
    											n_reg1_6=n_reg1_6(n_reg1_6~=loc_of_reg1 & n_reg1_6~=n_reg1_2(i2) & n_reg1_6~=n_reg1_3(i3) & n_reg1_6~=n_reg1_4(i4) & n_reg1_6~=n_reg1_5(i5));
    											sz6=size(n_reg1_6);
    											for i6=1:sz6(2)
    												if (n_reg1_6(i6) == loc_of_sec_reg)
    													display('A path to reg 4 is 6');
    													weightp6_reg4(t6)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6));
    													t6=t6+1;
    													else
															n_reg1_7=find(neighborss(n_reg1_6(i6),:));
															n_reg1_7=n_reg1_7(n_reg1_7~=loc_of_reg1 & n_reg1_7~=n_reg1_2(i2) & n_reg1_7~=n_reg1_3(i3) & n_reg1_7~=n_reg1_4(i4) & n_reg1_7~=n_reg1_5(i5) & n_reg1_7~=n_reg1_6(i6));
															sz7=size(n_reg1_7);
															for i7=1:sz7(2)
																if (n_reg1_7(i7) == loc_of_sec_reg)
																	display('A path to reg 4 is 7');
																	weightp7_reg4(t7)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7));
																	t7=t7+1;
																	else
																		n_reg1_8=find(neighborss(n_reg1_7(i7),:));
																		n_reg1_8=n_reg1_8(n_reg1_8~=loc_of_reg1 & n_reg1_8~=n_reg1_2(i2) & n_reg1_8~=n_reg1_3(i3) & n_reg1_8~=n_reg1_4(i4) & n_reg1_8~=n_reg1_5(i5) & n_reg1_8~=n_reg1_6(i6) & n_reg1_8~=n_reg1_7(i7));
																		sz8=size(n_reg1_8);
																		for i8=1:sz8(2)
																			if (n_reg1_8(i8) == loc_of_sec_reg)
																			display('A path to reg 4 is 8');
																			weightp8_reg4(t8)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8));
																			t8=t8+1;
																			else
																				n_reg1_9=find(neighborss(n_reg1_8(i8),:));
																				n_reg1_9=n_reg1_9(n_reg1_9~=loc_of_reg1 & n_reg1_9~=n_reg1_2(i2) & n_reg1_9~=n_reg1_3(i3) & n_reg1_9~=n_reg1_4(i4) & n_reg1_9~=n_reg1_5(i5) & n_reg1_9~=n_reg1_6(i6) & n_reg1_9~=n_reg1_7(i7) & n_reg1_9~=n_reg1_8(i8));
																				sz9=size(n_reg1_9);
																				for i9=1:sz9(2)
																					if (n_reg1_9(i9) == loc_of_sec_reg)
																						display('A path to reg 4 is 9');
																						weightp9_reg4(t9)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9));
																						t9=t9+1;
																						else
																							n_reg1_10=find(neighborss(n_reg1_9(i9),:));
																							n_reg1_10=n_reg1_10(n_reg1_10~=loc_of_reg1 & n_reg1_10~=n_reg1_2(i2) & n_reg1_10~=n_reg1_3(i3) & n_reg1_10~=n_reg1_4(i4) & n_reg1_10~=n_reg1_5(i5) & n_reg1_10~=n_reg1_6(i6) & n_reg1_10~=n_reg1_7(i7) & n_reg1_10~=n_reg1_8(i8) & n_reg1_10~=n_reg1_9(i9));
																							sz10=size(n_reg1_10);
																							for i10=1:sz10(2)
																								if (n_reg1_10(i10) == loc_of_sec_reg)
																									display('A path to reg 4 is 10');
																									weightp10_reg4(t10)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10));
																									t10=t10+1;
																									else
																										n_reg1_11=find(neighborss(n_reg1_10(i10),:));
																										n_reg1_11=n_reg1_11(n_reg1_11~=loc_of_reg1 & n_reg1_11~=n_reg1_2(i2) & n_reg1_11~=n_reg1_3(i3) & n_reg1_11~=n_reg1_4(i4) & n_reg1_11~=n_reg1_5(i5) & n_reg1_11~=n_reg1_6(i6) & n_reg1_11~=n_reg1_7(i7) & n_reg1_11~=n_reg1_8(i8) & n_reg1_11~=n_reg1_9(i9) & n_reg1_11~=n_reg1_10(i10));
																										sz11=size(n_reg1_11);
																										for i11=1:sz11(2)
																											if (n_reg1_11(i11) == loc_of_sec_reg)
																												display('A path to reg 4 is 11');
																												weightp11_reg4(t11)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11));
																												t11=t11+1;
																												else
																													n_reg1_12=find(neighborss(n_reg1_11(i11),:));
																													n_reg1_12=n_reg1_12(n_reg1_12~=loc_of_reg1 & n_reg1_12~=n_reg1_2(i2) & n_reg1_12~=n_reg1_3(i3) & n_reg1_12~=n_reg1_4(i4) & n_reg1_12~=n_reg1_5(i5) & n_reg1_12~=n_reg1_6(i6) & n_reg1_12~=n_reg1_7(i7) & n_reg1_12~=n_reg1_8(i8) & n_reg1_12~=n_reg1_9(i9) & n_reg1_12~=n_reg1_10(i10) & n_reg1_12~=n_reg1_11(i11));
																													sz12=size(n_reg1_12);
																													for i12=1:sz12(2)
																														if (n_reg1_12(i12) == loc_of_sec_reg)
																															display('A path to reg 4 is 12');
																															weightp12_reg4(t12)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12));
																															t12=t12+1;
																															else
																																n_reg1_13=find(neighborss(n_reg1_12(i12),:));
																																n_reg1_13=n_reg1_13(n_reg1_13~=loc_of_reg1 & n_reg1_13~=n_reg1_2(i2) & n_reg1_13~=n_reg1_3(i3) & n_reg1_13~=n_reg1_4(i4) & n_reg1_13~=n_reg1_5(i5) & n_reg1_13~=n_reg1_6(i6) & n_reg1_13~=n_reg1_7(i7) & n_reg1_13~=n_reg1_8(i8) & n_reg1_13~=n_reg1_9(i9) & n_reg1_13~=n_reg1_10(i10) & n_reg1_13~=n_reg1_11(i11) & n_reg1_13~=n_reg1_12(i12));
																																sz13=size(n_reg1_13);
																																for i13=1:sz13(2)
																																	if (n_reg1_13(i13) == loc_of_sec_reg)
																																		display('A path to reg 4 is 13');
																																		weightp13_reg4(t13)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n_reg1_14=find(neighborss(n_reg1_13(i13),:));
																																		% 	n_reg1_14=n_reg1_14(n_reg1_14~=loc_of_reg1 & n_reg1_14~=n_reg1_2(i2) & n_reg1_14~=n_reg1_3(i3) & n_reg1_14~=n_reg1_4(i4) & n_reg1_14~=n_reg1_5(i5) & n_reg1_14~=n_reg1_6(i6) & n_reg1_14~=n_reg1_7(i7) & n_reg1_14~=n_reg1_8(i8) & n_reg1_14~=n_reg1_9(i9) & n_reg1_14~=n_reg1_10(i10) & n_reg1_14~=n_reg1_11(i11) & n_reg1_14~=n_reg1_12(i12) & n_reg1_14~=n_reg1_13(i13));
																																		% 	sz14=size(n_reg1_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n_reg1_14(i14) == loc_of_sec_reg)
																																		% 			display('A path to reg 4 is 14');
																																		% 			weightp14_reg4(t14)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n_reg1_15=find(neighborss(n_reg1_14(i14),:));
																												 					% 					n_reg1_15=n_reg1_15(n_reg1_15~=loc_of_reg1 & n_reg1_15~=n_reg1_2(i2) & n_reg1_15~=n_reg1_3(i3) & n_reg1_15~=n_reg1_4(i4) & n_reg1_15~=n_reg1_5(i5) & n_reg1_15~=n_reg1_6(i6) & n_reg1_15~=n_reg1_7(i7) & n_reg1_15~=n_reg1_8(i8) & n_reg1_15~=n_reg1_9(i9) & n_reg1_15~=n_reg1_10(i10) & n_reg1_15~=n_reg1_11(i11) & n_reg1_15~=n_reg1_12(i12) & n_reg1_15~=n_reg1_13(i13) & n_reg1_15~=n_reg1_14(i14));
																												 					% 					sz15=size(n_reg1_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n_reg1_15(i15) == loc_of_sec_reg)
																												 					% 							display('A path to reg 4 is 15');
																												 					% 							weightp15_reg4(t15)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14))*lev(n_reg1_14(i14),n_reg1_15(i15));
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
loc_of_sec_reg = loc_of_sec_reg_arr(4);

    	for i1=1:sz(2)
    		if (n_reg1(i1) == loc_of_sec_reg)
    			display('A path to reg 5 is 1');
    			weightp1_reg5(t1)=lev(loc_of_reg1,n_reg1(i1));
    			t1=t1+1;
    		else
    			n_reg1_2=find(neighborss(n_reg1(i1),:));
    			n_reg1_2=n_reg1_2(n_reg1_2~=loc_of_reg1);
    			sz2=size(n_reg1_2);
    			for i2=1:sz2(2)
    				if (n_reg1_2(i2) == loc_of_sec_reg)
    					display('A path to reg 5 is 2');
    					weightp2_reg5(t2)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2));
    					t2=t2+1;
    				else
    					n_reg1_3=find(neighborss(n_reg1_2(i2),:));
    					n_reg1_3=n_reg1_3(n_reg1_3~=n_reg1_2(i2) & n_reg1_3~=loc_of_reg1);
    					sz3=size(n_reg1_3);
    					for i3=1:sz3(2)
    						if (n_reg1_3(i3) == loc_of_sec_reg)
    							display('A path to reg 5 is 3');
    							weightp3_reg5(t3)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3));
    							t3=t3+1;
    						else
    							n_reg1_4=find(neighborss(n_reg1_3(i3),:));
    							n_reg1_4=n_reg1_4(n_reg1_4~=n_reg1_2(i2) & n_reg1_4~=n_reg1_3(i3) & n_reg1_4~=loc_of_reg1);
    							sz4=size(n_reg1_4);
    							for i4=1:sz4(2)
    								if (n_reg1_4(i4) == loc_of_sec_reg)
    								display('A path to reg 5 is 4');
    								weightp4_reg5(t4)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4));
    								t4=t4+1;
        							else
    									n_reg1_5=find(neighborss(n_reg1_4(i4),:));
    									n_reg1_5=n_reg1_5(n_reg1_5~=n_reg1_2(i2) & n_reg1_5~=n_reg1_3(i3) & n_reg1_5~=n_reg1_4(i4) & n_reg1_5~=loc_of_reg1);
    									sz5=size(n_reg1_5);
    									for i5=1:sz5(2)
    										if (n_reg1_5(i5) == loc_of_sec_reg)
    										display('A path to reg 5 is 5');
    										weightp5_reg5(t5)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5));
    										t5=t5+1;
    										else
    											n_reg1_6=find(neighborss(n_reg1_5(i5),:));
    											n_reg1_6=n_reg1_6(n_reg1_6~=loc_of_reg1 & n_reg1_6~=n_reg1_2(i2) & n_reg1_6~=n_reg1_3(i3) & n_reg1_6~=n_reg1_4(i4) & n_reg1_6~=n_reg1_5(i5));
    											sz6=size(n_reg1_6);
    											for i6=1:sz6(2)
    												if (n_reg1_6(i6) == loc_of_sec_reg)
    													display('A path to reg 5 is 6');
    													weightp6_reg5(t6)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6));
    													t6=t6+1;
    													else
															n_reg1_7=find(neighborss(n_reg1_6(i6),:));
															n_reg1_7=n_reg1_7(n_reg1_7~=loc_of_reg1 & n_reg1_7~=n_reg1_2(i2) & n_reg1_7~=n_reg1_3(i3) & n_reg1_7~=n_reg1_4(i4) & n_reg1_7~=n_reg1_5(i5) & n_reg1_7~=n_reg1_6(i6));
															sz7=size(n_reg1_7);
															for i7=1:sz7(2)
																if (n_reg1_7(i7) == loc_of_sec_reg)
																	display('A path to reg 5 is 7');
																	weightp7_reg5(t7)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7));
																	t7=t7+1;
																	else
																		n_reg1_8=find(neighborss(n_reg1_7(i7),:));
																		n_reg1_8=n_reg1_8(n_reg1_8~=loc_of_reg1 & n_reg1_8~=n_reg1_2(i2) & n_reg1_8~=n_reg1_3(i3) & n_reg1_8~=n_reg1_4(i4) & n_reg1_8~=n_reg1_5(i5) & n_reg1_8~=n_reg1_6(i6) & n_reg1_8~=n_reg1_7(i7));
																		sz8=size(n_reg1_8);
																		for i8=1:sz8(2)
																			if (n_reg1_8(i8) == loc_of_sec_reg)
																			display('A path to reg 5 is 8');
																			weightp8_reg5(t8)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8));
																			t8=t8+1;
																			else
																				n_reg1_9=find(neighborss(n_reg1_8(i8),:));
																				n_reg1_9=n_reg1_9(n_reg1_9~=loc_of_reg1 & n_reg1_9~=n_reg1_2(i2) & n_reg1_9~=n_reg1_3(i3) & n_reg1_9~=n_reg1_4(i4) & n_reg1_9~=n_reg1_5(i5) & n_reg1_9~=n_reg1_6(i6) & n_reg1_9~=n_reg1_7(i7) & n_reg1_9~=n_reg1_8(i8));
																				sz9=size(n_reg1_9);
																				for i9=1:sz9(2)
																					if (n_reg1_9(i9) == loc_of_sec_reg)
																						display('A path to reg 5 is 9');
																						weightp9_reg5(t9)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9));
																						t9=t9+1;
																						else
																							n_reg1_10=find(neighborss(n_reg1_9(i9),:));
																							n_reg1_10=n_reg1_10(n_reg1_10~=loc_of_reg1 & n_reg1_10~=n_reg1_2(i2) & n_reg1_10~=n_reg1_3(i3) & n_reg1_10~=n_reg1_4(i4) & n_reg1_10~=n_reg1_5(i5) & n_reg1_10~=n_reg1_6(i6) & n_reg1_10~=n_reg1_7(i7) & n_reg1_10~=n_reg1_8(i8) & n_reg1_10~=n_reg1_9(i9));
																							sz10=size(n_reg1_10);
																							for i10=1:sz10(2)
																								if (n_reg1_10(i10) == loc_of_sec_reg)
																									display('A path to reg 5 is 10');
																									weightp10_reg5(t10)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10));
																									t10=t10+1;
																									else
																										n_reg1_11=find(neighborss(n_reg1_10(i10),:));
																										n_reg1_11=n_reg1_11(n_reg1_11~=loc_of_reg1 & n_reg1_11~=n_reg1_2(i2) & n_reg1_11~=n_reg1_3(i3) & n_reg1_11~=n_reg1_4(i4) & n_reg1_11~=n_reg1_5(i5) & n_reg1_11~=n_reg1_6(i6) & n_reg1_11~=n_reg1_7(i7) & n_reg1_11~=n_reg1_8(i8) & n_reg1_11~=n_reg1_9(i9) & n_reg1_11~=n_reg1_10(i10));
																										sz11=size(n_reg1_11);
																										for i11=1:sz11(2)
																											if (n_reg1_11(i11) == loc_of_sec_reg)
																												display('A path to reg 5 is 11');
																												weightp11_reg5(t11)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11));
																												t11=t11+1;
																												else
																													n_reg1_12=find(neighborss(n_reg1_11(i11),:));
																													n_reg1_12=n_reg1_12(n_reg1_12~=loc_of_reg1 & n_reg1_12~=n_reg1_2(i2) & n_reg1_12~=n_reg1_3(i3) & n_reg1_12~=n_reg1_4(i4) & n_reg1_12~=n_reg1_5(i5) & n_reg1_12~=n_reg1_6(i6) & n_reg1_12~=n_reg1_7(i7) & n_reg1_12~=n_reg1_8(i8) & n_reg1_12~=n_reg1_9(i9) & n_reg1_12~=n_reg1_10(i10) & n_reg1_12~=n_reg1_11(i11));
																													sz12=size(n_reg1_12);
																													for i12=1:sz12(2)
																														if (n_reg1_12(i12) == loc_of_sec_reg)
																															display('A path to reg 5 is 12');
																															weightp12_reg5(t12)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12));
																															t12=t12+1;
																															else
																																n_reg1_13=find(neighborss(n_reg1_12(i12),:));
																																n_reg1_13=n_reg1_13(n_reg1_13~=loc_of_reg1 & n_reg1_13~=n_reg1_2(i2) & n_reg1_13~=n_reg1_3(i3) & n_reg1_13~=n_reg1_4(i4) & n_reg1_13~=n_reg1_5(i5) & n_reg1_13~=n_reg1_6(i6) & n_reg1_13~=n_reg1_7(i7) & n_reg1_13~=n_reg1_8(i8) & n_reg1_13~=n_reg1_9(i9) & n_reg1_13~=n_reg1_10(i10) & n_reg1_13~=n_reg1_11(i11) & n_reg1_13~=n_reg1_12(i12));
																																sz13=size(n_reg1_13);
																																for i13=1:sz13(2)
																																	if (n_reg1_13(i13) == loc_of_sec_reg)
																																		display('A path to reg 5 is 13');
																																		weightp13_reg5(t13)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n_reg1_14=find(neighborss(n_reg1_13(i13),:));
																																		% 	n_reg1_14=n_reg1_14(n_reg1_14~=loc_of_reg1 & n_reg1_14~=n_reg1_2(i2) & n_reg1_14~=n_reg1_3(i3) & n_reg1_14~=n_reg1_4(i4) & n_reg1_14~=n_reg1_5(i5) & n_reg1_14~=n_reg1_6(i6) & n_reg1_14~=n_reg1_7(i7) & n_reg1_14~=n_reg1_8(i8) & n_reg1_14~=n_reg1_9(i9) & n_reg1_14~=n_reg1_10(i10) & n_reg1_14~=n_reg1_11(i11) & n_reg1_14~=n_reg1_12(i12) & n_reg1_14~=n_reg1_13(i13));
																																		% 	sz14=size(n_reg1_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n_reg1_14(i14) == loc_of_sec_reg)
																																		% 			display('A path to reg 5 is 14');
																																		% 			weightp14_reg5(t14)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n_reg1_15=find(neighborss(n_reg1_14(i14),:));
																												 					% 					n_reg1_15=n_reg1_15(n_reg1_15~=loc_of_reg1 & n_reg1_15~=n_reg1_2(i2) & n_reg1_15~=n_reg1_3(i3) & n_reg1_15~=n_reg1_4(i4) & n_reg1_15~=n_reg1_5(i5) & n_reg1_15~=n_reg1_6(i6) & n_reg1_15~=n_reg1_7(i7) & n_reg1_15~=n_reg1_8(i8) & n_reg1_15~=n_reg1_9(i9) & n_reg1_15~=n_reg1_10(i10) & n_reg1_15~=n_reg1_11(i11) & n_reg1_15~=n_reg1_12(i12) & n_reg1_15~=n_reg1_13(i13) & n_reg1_15~=n_reg1_14(i14));
																												 					% 					sz15=size(n_reg1_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n_reg1_15(i15) == loc_of_sec_reg)
																												 					% 							display('A path to reg 5 is 15');
																												 					% 							weightp15_reg5(t15)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14))*lev(n_reg1_14(i14),n_reg1_15(i15));
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
loc_of_sec_reg = loc_of_sec_reg_arr(5);

    	for i1=1:sz(2)
    		if (n_reg1(i1) == loc_of_sec_reg)
    			display('A path to reg 6 is 1');
    			weightp1_reg6(t1)=lev(loc_of_reg1,n_reg1(i1));
    			t1=t1+1;
    		else
    			n_reg1_2=find(neighborss(n_reg1(i1),:));
    			n_reg1_2=n_reg1_2(n_reg1_2~=loc_of_reg1);
    			sz2=size(n_reg1_2);
    			for i2=1:sz2(2)
    				if (n_reg1_2(i2) == loc_of_sec_reg)
    					display('A path to reg 6 is 2');
    					weightp2_reg6(t2)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2));
    					t2=t2+1;
    				else
    					n_reg1_3=find(neighborss(n_reg1_2(i2),:));
    					n_reg1_3=n_reg1_3(n_reg1_3~=n_reg1_2(i2) & n_reg1_3~=loc_of_reg1);
    					sz3=size(n_reg1_3);
    					for i3=1:sz3(2)
    						if (n_reg1_3(i3) == loc_of_sec_reg)
    							display('A path to reg 6 is 3');
    							weightp3_reg6(t3)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3));
    							t3=t3+1;
    						else
    							n_reg1_4=find(neighborss(n_reg1_3(i3),:));
    							n_reg1_4=n_reg1_4(n_reg1_4~=n_reg1_2(i2) & n_reg1_4~=n_reg1_3(i3) & n_reg1_4~=loc_of_reg1);
    							sz4=size(n_reg1_4);
    							for i4=1:sz4(2)
    								if (n_reg1_4(i4) == loc_of_sec_reg)
    								display('A path to reg 6 is 4');
    								weightp4_reg6(t4)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4));
    								t4=t4+1;
        							else
    									n_reg1_5=find(neighborss(n_reg1_4(i4),:));
    									n_reg1_5=n_reg1_5(n_reg1_5~=n_reg1_2(i2) & n_reg1_5~=n_reg1_3(i3) & n_reg1_5~=n_reg1_4(i4) & n_reg1_5~=loc_of_reg1);
    									sz5=size(n_reg1_5);
    									for i5=1:sz5(2)
    										if (n_reg1_5(i5) == loc_of_sec_reg)
    										display('A path to reg 6 is 5');
    										weightp5_reg6(t5)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5));
    										t5=t5+1;
    										else
    											n_reg1_6=find(neighborss(n_reg1_5(i5),:));
    											n_reg1_6=n_reg1_6(n_reg1_6~=loc_of_reg1 & n_reg1_6~=n_reg1_2(i2) & n_reg1_6~=n_reg1_3(i3) & n_reg1_6~=n_reg1_4(i4) & n_reg1_6~=n_reg1_5(i5));
    											sz6=size(n_reg1_6);
    											for i6=1:sz6(2)
    												if (n_reg1_6(i6) == loc_of_sec_reg)
    													display('A path to reg 6 is 6');
    													weightp6_reg6(t6)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6));
    													t6=t6+1;
    													else
															n_reg1_7=find(neighborss(n_reg1_6(i6),:));
															n_reg1_7=n_reg1_7(n_reg1_7~=loc_of_reg1 & n_reg1_7~=n_reg1_2(i2) & n_reg1_7~=n_reg1_3(i3) & n_reg1_7~=n_reg1_4(i4) & n_reg1_7~=n_reg1_5(i5) & n_reg1_7~=n_reg1_6(i6));
															sz7=size(n_reg1_7);
															for i7=1:sz7(2)
																if (n_reg1_7(i7) == loc_of_sec_reg)
																	display('A path to reg 6 is 7');
																	weightp7_reg6(t7)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7));
																	t7=t7+1;
																	else
																		n_reg1_8=find(neighborss(n_reg1_7(i7),:));
																		n_reg1_8=n_reg1_8(n_reg1_8~=loc_of_reg1 & n_reg1_8~=n_reg1_2(i2) & n_reg1_8~=n_reg1_3(i3) & n_reg1_8~=n_reg1_4(i4) & n_reg1_8~=n_reg1_5(i5) & n_reg1_8~=n_reg1_6(i6) & n_reg1_8~=n_reg1_7(i7));
																		sz8=size(n_reg1_8);
																		for i8=1:sz8(2)
																			if (n_reg1_8(i8) == loc_of_sec_reg)
																			display('A path to reg 6 is 8');
																			weightp8_reg6(t8)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8));
																			t8=t8+1;
																			else
																				n_reg1_9=find(neighborss(n_reg1_8(i8),:));
																				n_reg1_9=n_reg1_9(n_reg1_9~=loc_of_reg1 & n_reg1_9~=n_reg1_2(i2) & n_reg1_9~=n_reg1_3(i3) & n_reg1_9~=n_reg1_4(i4) & n_reg1_9~=n_reg1_5(i5) & n_reg1_9~=n_reg1_6(i6) & n_reg1_9~=n_reg1_7(i7) & n_reg1_9~=n_reg1_8(i8));
																				sz9=size(n_reg1_9);
																				for i9=1:sz9(2)
																					if (n_reg1_9(i9) == loc_of_sec_reg)
																						display('A path to reg 6 is 9');
																						weightp9_reg6(t9)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9));
																						t9=t9+1;
																						else
																							n_reg1_10=find(neighborss(n_reg1_9(i9),:));
																							n_reg1_10=n_reg1_10(n_reg1_10~=loc_of_reg1 & n_reg1_10~=n_reg1_2(i2) & n_reg1_10~=n_reg1_3(i3) & n_reg1_10~=n_reg1_4(i4) & n_reg1_10~=n_reg1_5(i5) & n_reg1_10~=n_reg1_6(i6) & n_reg1_10~=n_reg1_7(i7) & n_reg1_10~=n_reg1_8(i8) & n_reg1_10~=n_reg1_9(i9));
																							sz10=size(n_reg1_10);
																							for i10=1:sz10(2)
																								if (n_reg1_10(i10) == loc_of_sec_reg)
																									display('A path to reg 6 is 10');
																									weightp10_reg6(t10)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10));
																									t10=t10+1;
																									else
																										n_reg1_11=find(neighborss(n_reg1_10(i10),:));
																										n_reg1_11=n_reg1_11(n_reg1_11~=loc_of_reg1 & n_reg1_11~=n_reg1_2(i2) & n_reg1_11~=n_reg1_3(i3) & n_reg1_11~=n_reg1_4(i4) & n_reg1_11~=n_reg1_5(i5) & n_reg1_11~=n_reg1_6(i6) & n_reg1_11~=n_reg1_7(i7) & n_reg1_11~=n_reg1_8(i8) & n_reg1_11~=n_reg1_9(i9) & n_reg1_11~=n_reg1_10(i10));
																										sz11=size(n_reg1_11);
																										for i11=1:sz11(2)
																											if (n_reg1_11(i11) == loc_of_sec_reg)
																												display('A path to reg 6 is 11');
																												weightp11_reg6(t11)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11));
																												t11=t11+1;
																												else
																													n_reg1_12=find(neighborss(n_reg1_11(i11),:));
																													n_reg1_12=n_reg1_12(n_reg1_12~=loc_of_reg1 & n_reg1_12~=n_reg1_2(i2) & n_reg1_12~=n_reg1_3(i3) & n_reg1_12~=n_reg1_4(i4) & n_reg1_12~=n_reg1_5(i5) & n_reg1_12~=n_reg1_6(i6) & n_reg1_12~=n_reg1_7(i7) & n_reg1_12~=n_reg1_8(i8) & n_reg1_12~=n_reg1_9(i9) & n_reg1_12~=n_reg1_10(i10) & n_reg1_12~=n_reg1_11(i11));
																													sz12=size(n_reg1_12);
																													for i12=1:sz12(2)
																														if (n_reg1_12(i12) == loc_of_sec_reg)
																															display('A path to reg 6 is 12');
																															weightp12_reg6(t12)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12));
																															t12=t12+1;
																															else
																																n_reg1_13=find(neighborss(n_reg1_12(i12),:));
																																n_reg1_13=n_reg1_13(n_reg1_13~=loc_of_reg1 & n_reg1_13~=n_reg1_2(i2) & n_reg1_13~=n_reg1_3(i3) & n_reg1_13~=n_reg1_4(i4) & n_reg1_13~=n_reg1_5(i5) & n_reg1_13~=n_reg1_6(i6) & n_reg1_13~=n_reg1_7(i7) & n_reg1_13~=n_reg1_8(i8) & n_reg1_13~=n_reg1_9(i9) & n_reg1_13~=n_reg1_10(i10) & n_reg1_13~=n_reg1_11(i11) & n_reg1_13~=n_reg1_12(i12));
																																sz13=size(n_reg1_13);
																																for i13=1:sz13(2)
																																	if (n_reg1_13(i13) == loc_of_sec_reg)
																																		display('A path to reg 6 is 13');
																																		weightp13_reg6(t13)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n_reg1_14=find(neighborss(n_reg1_13(i13),:));
																																		% 	n_reg1_14=n_reg1_14(n_reg1_14~=loc_of_reg1 & n_reg1_14~=n_reg1_2(i2) & n_reg1_14~=n_reg1_3(i3) & n_reg1_14~=n_reg1_4(i4) & n_reg1_14~=n_reg1_5(i5) & n_reg1_14~=n_reg1_6(i6) & n_reg1_14~=n_reg1_7(i7) & n_reg1_14~=n_reg1_8(i8) & n_reg1_14~=n_reg1_9(i9) & n_reg1_14~=n_reg1_10(i10) & n_reg1_14~=n_reg1_11(i11) & n_reg1_14~=n_reg1_12(i12) & n_reg1_14~=n_reg1_13(i13));
																																		% 	sz14=size(n_reg1_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n_reg1_14(i14) == loc_of_sec_reg)
																																		% 			display('A path to reg 6 is 14');
																																		% 			weightp14_reg6(t14)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n_reg1_15=find(neighborss(n_reg1_14(i14),:));
																												 					% 					n_reg1_15=n_reg1_15(n_reg1_15~=loc_of_reg1 & n_reg1_15~=n_reg1_2(i2) & n_reg1_15~=n_reg1_3(i3) & n_reg1_15~=n_reg1_4(i4) & n_reg1_15~=n_reg1_5(i5) & n_reg1_15~=n_reg1_6(i6) & n_reg1_15~=n_reg1_7(i7) & n_reg1_15~=n_reg1_8(i8) & n_reg1_15~=n_reg1_9(i9) & n_reg1_15~=n_reg1_10(i10) & n_reg1_15~=n_reg1_11(i11) & n_reg1_15~=n_reg1_12(i12) & n_reg1_15~=n_reg1_13(i13) & n_reg1_15~=n_reg1_14(i14));
																												 					% 					sz15=size(n_reg1_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n_reg1_15(i15) == loc_of_sec_reg)
																												 					% 							display('A path to reg 6 is 15');
																												 					% 							weightp15_reg6(t15)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14))*lev(n_reg1_14(i14),n_reg1_15(i15));
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
loc_of_sec_reg = loc_of_sec_reg_arr(6);

    	for i1=1:sz(2)
    		if (n_reg1(i1) == loc_of_sec_reg)
    			display('A path to reg 7 is 1');
    			weightp1_reg7(t1)=lev(loc_of_reg1,n_reg1(i1));
    			t1=t1+1;
    		else
    			n_reg1_2=find(neighborss(n_reg1(i1),:));
    			n_reg1_2=n_reg1_2(n_reg1_2~=loc_of_reg1);
    			sz2=size(n_reg1_2);
    			for i2=1:sz2(2)
    				if (n_reg1_2(i2) == loc_of_sec_reg)
    					display('A path to reg 7 is 2');
    					weightp2_reg7(t2)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2));
    					t2=t2+1;
    				else
    					n_reg1_3=find(neighborss(n_reg1_2(i2),:));
    					n_reg1_3=n_reg1_3(n_reg1_3~=n_reg1_2(i2) & n_reg1_3~=loc_of_reg1);
    					sz3=size(n_reg1_3);
    					for i3=1:sz3(2)
    						if (n_reg1_3(i3) == loc_of_sec_reg)
    							display('A path to reg 7 is 3');
    							weightp3_reg7(t3)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3));
    							t3=t3+1;
    						else
    							n_reg1_4=find(neighborss(n_reg1_3(i3),:));
    							n_reg1_4=n_reg1_4(n_reg1_4~=n_reg1_2(i2) & n_reg1_4~=n_reg1_3(i3) & n_reg1_4~=loc_of_reg1);
    							sz4=size(n_reg1_4);
    							for i4=1:sz4(2)
    								if (n_reg1_4(i4) == loc_of_sec_reg)
    								display('A path to reg 7 is 4');
    								weightp4_reg7(t4)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4));
    								t4=t4+1;
        							else
    									n_reg1_5=find(neighborss(n_reg1_4(i4),:));
    									n_reg1_5=n_reg1_5(n_reg1_5~=n_reg1_2(i2) & n_reg1_5~=n_reg1_3(i3) & n_reg1_5~=n_reg1_4(i4) & n_reg1_5~=loc_of_reg1);
    									sz5=size(n_reg1_5);
    									for i5=1:sz5(2)
    										if (n_reg1_5(i5) == loc_of_sec_reg)
    										display('A path to reg 7 is 5');
    										weightp5_reg7(t5)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5));
    										t5=t5+1;
    										else
    											n_reg1_6=find(neighborss(n_reg1_5(i5),:));
    											n_reg1_6=n_reg1_6(n_reg1_6~=loc_of_reg1 & n_reg1_6~=n_reg1_2(i2) & n_reg1_6~=n_reg1_3(i3) & n_reg1_6~=n_reg1_4(i4) & n_reg1_6~=n_reg1_5(i5));
    											sz6=size(n_reg1_6);
    											for i6=1:sz6(2)
    												if (n_reg1_6(i6) == loc_of_sec_reg)
    													display('A path to reg 7 is 6');
    													weightp6_reg7(t6)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6));
    													t6=t6+1;
    													else
															n_reg1_7=find(neighborss(n_reg1_6(i6),:));
															n_reg1_7=n_reg1_7(n_reg1_7~=loc_of_reg1 & n_reg1_7~=n_reg1_2(i2) & n_reg1_7~=n_reg1_3(i3) & n_reg1_7~=n_reg1_4(i4) & n_reg1_7~=n_reg1_5(i5) & n_reg1_7~=n_reg1_6(i6));
															sz7=size(n_reg1_7);
															for i7=1:sz7(2)
																if (n_reg1_7(i7) == loc_of_sec_reg)
																	display('A path to reg 7 is 7');
																	weightp7_reg7(t7)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7));
																	t7=t7+1;
																	else
																		n_reg1_8=find(neighborss(n_reg1_7(i7),:));
																		n_reg1_8=n_reg1_8(n_reg1_8~=loc_of_reg1 & n_reg1_8~=n_reg1_2(i2) & n_reg1_8~=n_reg1_3(i3) & n_reg1_8~=n_reg1_4(i4) & n_reg1_8~=n_reg1_5(i5) & n_reg1_8~=n_reg1_6(i6) & n_reg1_8~=n_reg1_7(i7));
																		sz8=size(n_reg1_8);
																		for i8=1:sz8(2)
																			if (n_reg1_8(i8) == loc_of_sec_reg)
																			display('A path to reg 7 is 8');
																			weightp8_reg7(t8)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8));
																			t8=t8+1;
																			else
																				n_reg1_9=find(neighborss(n_reg1_8(i8),:));
																				n_reg1_9=n_reg1_9(n_reg1_9~=loc_of_reg1 & n_reg1_9~=n_reg1_2(i2) & n_reg1_9~=n_reg1_3(i3) & n_reg1_9~=n_reg1_4(i4) & n_reg1_9~=n_reg1_5(i5) & n_reg1_9~=n_reg1_6(i6) & n_reg1_9~=n_reg1_7(i7) & n_reg1_9~=n_reg1_8(i8));
																				sz9=size(n_reg1_9);
																				for i9=1:sz9(2)
																					if (n_reg1_9(i9) == loc_of_sec_reg)
																						display('A path to reg 7 is 9');
																						weightp9_reg7(t9)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9));
																						t9=t9+1;
																						else
																							n_reg1_10=find(neighborss(n_reg1_9(i9),:));
																							n_reg1_10=n_reg1_10(n_reg1_10~=loc_of_reg1 & n_reg1_10~=n_reg1_2(i2) & n_reg1_10~=n_reg1_3(i3) & n_reg1_10~=n_reg1_4(i4) & n_reg1_10~=n_reg1_5(i5) & n_reg1_10~=n_reg1_6(i6) & n_reg1_10~=n_reg1_7(i7) & n_reg1_10~=n_reg1_8(i8) & n_reg1_10~=n_reg1_9(i9));
																							sz10=size(n_reg1_10);
																							for i10=1:sz10(2)
																								if (n_reg1_10(i10) == loc_of_sec_reg)
																									display('A path to reg 7 is 10');
																									weightp10_reg7(t10)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10));
																									t10=t10+1;
																									else
																										n_reg1_11=find(neighborss(n_reg1_10(i10),:));
																										n_reg1_11=n_reg1_11(n_reg1_11~=loc_of_reg1 & n_reg1_11~=n_reg1_2(i2) & n_reg1_11~=n_reg1_3(i3) & n_reg1_11~=n_reg1_4(i4) & n_reg1_11~=n_reg1_5(i5) & n_reg1_11~=n_reg1_6(i6) & n_reg1_11~=n_reg1_7(i7) & n_reg1_11~=n_reg1_8(i8) & n_reg1_11~=n_reg1_9(i9) & n_reg1_11~=n_reg1_10(i10));
																										sz11=size(n_reg1_11);
																										for i11=1:sz11(2)
																											if (n_reg1_11(i11) == loc_of_sec_reg)
																												display('A path to reg 7 is 11');
																												weightp11_reg7(t11)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11));
																												t11=t11+1;
																												else
																													n_reg1_12=find(neighborss(n_reg1_11(i11),:));
																													n_reg1_12=n_reg1_12(n_reg1_12~=loc_of_reg1 & n_reg1_12~=n_reg1_2(i2) & n_reg1_12~=n_reg1_3(i3) & n_reg1_12~=n_reg1_4(i4) & n_reg1_12~=n_reg1_5(i5) & n_reg1_12~=n_reg1_6(i6) & n_reg1_12~=n_reg1_7(i7) & n_reg1_12~=n_reg1_8(i8) & n_reg1_12~=n_reg1_9(i9) & n_reg1_12~=n_reg1_10(i10) & n_reg1_12~=n_reg1_11(i11));
																													sz12=size(n_reg1_12);
																													for i12=1:sz12(2)
																														if (n_reg1_12(i12) == loc_of_sec_reg)
																															display('A path to reg 7 is 12');
																															weightp12_reg7(t12)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12));
																															t12=t12+1;
																															else
																																n_reg1_13=find(neighborss(n_reg1_12(i12),:));
																																n_reg1_13=n_reg1_13(n_reg1_13~=loc_of_reg1 & n_reg1_13~=n_reg1_2(i2) & n_reg1_13~=n_reg1_3(i3) & n_reg1_13~=n_reg1_4(i4) & n_reg1_13~=n_reg1_5(i5) & n_reg1_13~=n_reg1_6(i6) & n_reg1_13~=n_reg1_7(i7) & n_reg1_13~=n_reg1_8(i8) & n_reg1_13~=n_reg1_9(i9) & n_reg1_13~=n_reg1_10(i10) & n_reg1_13~=n_reg1_11(i11) & n_reg1_13~=n_reg1_12(i12));
																																sz13=size(n_reg1_13);
																																for i13=1:sz13(2)
																																	if (n_reg1_13(i13) == loc_of_sec_reg)
																																		display('A path to reg 7 is 13');
																																		weightp13_reg7(t13)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n_reg1_14=find(neighborss(n_reg1_13(i13),:));
																																		% 	n_reg1_14=n_reg1_14(n_reg1_14~=loc_of_reg1 & n_reg1_14~=n_reg1_2(i2) & n_reg1_14~=n_reg1_3(i3) & n_reg1_14~=n_reg1_4(i4) & n_reg1_14~=n_reg1_5(i5) & n_reg1_14~=n_reg1_6(i6) & n_reg1_14~=n_reg1_7(i7) & n_reg1_14~=n_reg1_8(i8) & n_reg1_14~=n_reg1_9(i9) & n_reg1_14~=n_reg1_10(i10) & n_reg1_14~=n_reg1_11(i11) & n_reg1_14~=n_reg1_12(i12) & n_reg1_14~=n_reg1_13(i13));
																																		% 	sz14=size(n_reg1_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n_reg1_14(i14) == loc_of_sec_reg)
																																		% 			display('A path to reg 7 is 14');
																																		% 			weightp14_reg7(t14)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n_reg1_15=find(neighborss(n_reg1_14(i14),:));
																												 					% 					n_reg1_15=n_reg1_15(n_reg1_15~=loc_of_reg1 & n_reg1_15~=n_reg1_2(i2) & n_reg1_15~=n_reg1_3(i3) & n_reg1_15~=n_reg1_4(i4) & n_reg1_15~=n_reg1_5(i5) & n_reg1_15~=n_reg1_6(i6) & n_reg1_15~=n_reg1_7(i7) & n_reg1_15~=n_reg1_8(i8) & n_reg1_15~=n_reg1_9(i9) & n_reg1_15~=n_reg1_10(i10) & n_reg1_15~=n_reg1_11(i11) & n_reg1_15~=n_reg1_12(i12) & n_reg1_15~=n_reg1_13(i13) & n_reg1_15~=n_reg1_14(i14));
																												 					% 					sz15=size(n_reg1_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n_reg1_15(i15) == loc_of_sec_reg)
																												 					% 							display('A path to reg 7 is 15');
																												 					% 							weightp15_reg7(t15)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14))*lev(n_reg1_14(i14),n_reg1_15(i15));
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
loc_of_sec_reg = loc_of_sec_reg_arr(7);

    	for i1=1:sz(2)
    		if (n_reg1(i1) == loc_of_sec_reg)
    			display('A path to reg 8 is 1');
    			weightp1_reg8(t1)=lev(loc_of_reg1,n_reg1(i1));
    			t1=t1+1;
    		else
    			n_reg1_2=find(neighborss(n_reg1(i1),:));
    			n_reg1_2=n_reg1_2(n_reg1_2~=loc_of_reg1);
    			sz2=size(n_reg1_2);
    			for i2=1:sz2(2)
    				if (n_reg1_2(i2) == loc_of_sec_reg)
    					display('A path to reg 8 is 2');
    					weightp2_reg8(t2)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2));
    					t2=t2+1;
    				else
    					n_reg1_3=find(neighborss(n_reg1_2(i2),:));
    					n_reg1_3=n_reg1_3(n_reg1_3~=n_reg1_2(i2) & n_reg1_3~=loc_of_reg1);
    					sz3=size(n_reg1_3);
    					for i3=1:sz3(2)
    						if (n_reg1_3(i3) == loc_of_sec_reg)
    							display('A path to reg 8 is 3');
    							weightp3_reg8(t3)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3));
    							t3=t3+1;
    						else
    							n_reg1_4=find(neighborss(n_reg1_3(i3),:));
    							n_reg1_4=n_reg1_4(n_reg1_4~=n_reg1_2(i2) & n_reg1_4~=n_reg1_3(i3) & n_reg1_4~=loc_of_reg1);
    							sz4=size(n_reg1_4);
    							for i4=1:sz4(2)
    								if (n_reg1_4(i4) == loc_of_sec_reg)
    								display('A path to reg 8 is 4');
    								weightp4_reg8(t4)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4));
    								t4=t4+1;
        							else
    									n_reg1_5=find(neighborss(n_reg1_4(i4),:));
    									n_reg1_5=n_reg1_5(n_reg1_5~=n_reg1_2(i2) & n_reg1_5~=n_reg1_3(i3) & n_reg1_5~=n_reg1_4(i4) & n_reg1_5~=loc_of_reg1);
    									sz5=size(n_reg1_5);
    									for i5=1:sz5(2)
    										if (n_reg1_5(i5) == loc_of_sec_reg)
    										display('A path to reg 8 is 5');
    										weightp5_reg8(t5)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5));
    										t5=t5+1;
    										else
    											n_reg1_6=find(neighborss(n_reg1_5(i5),:));
    											n_reg1_6=n_reg1_6(n_reg1_6~=loc_of_reg1 & n_reg1_6~=n_reg1_2(i2) & n_reg1_6~=n_reg1_3(i3) & n_reg1_6~=n_reg1_4(i4) & n_reg1_6~=n_reg1_5(i5));
    											sz6=size(n_reg1_6);
    											for i6=1:sz6(2)
    												if (n_reg1_6(i6) == loc_of_sec_reg)
    													display('A path to reg 8 is 6');
    													weightp6_reg8(t6)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6));
    													t6=t6+1;
    													else
															n_reg1_7=find(neighborss(n_reg1_6(i6),:));
															n_reg1_7=n_reg1_7(n_reg1_7~=loc_of_reg1 & n_reg1_7~=n_reg1_2(i2) & n_reg1_7~=n_reg1_3(i3) & n_reg1_7~=n_reg1_4(i4) & n_reg1_7~=n_reg1_5(i5) & n_reg1_7~=n_reg1_6(i6));
															sz7=size(n_reg1_7);
															for i7=1:sz7(2)
																if (n_reg1_7(i7) == loc_of_sec_reg)
																	display('A path to reg 8 is 7');
																	weightp7_reg8(t7)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7));
																	t7=t7+1;
																	else
																		n_reg1_8=find(neighborss(n_reg1_7(i7),:));
																		n_reg1_8=n_reg1_8(n_reg1_8~=loc_of_reg1 & n_reg1_8~=n_reg1_2(i2) & n_reg1_8~=n_reg1_3(i3) & n_reg1_8~=n_reg1_4(i4) & n_reg1_8~=n_reg1_5(i5) & n_reg1_8~=n_reg1_6(i6) & n_reg1_8~=n_reg1_7(i7));
																		sz8=size(n_reg1_8);
																		for i8=1:sz8(2)
																			if (n_reg1_8(i8) == loc_of_sec_reg)
																			display('A path to reg 8 is 8');
																			weightp8_reg8(t8)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8));
																			t8=t8+1;
																			else
																				n_reg1_9=find(neighborss(n_reg1_8(i8),:));
																				n_reg1_9=n_reg1_9(n_reg1_9~=loc_of_reg1 & n_reg1_9~=n_reg1_2(i2) & n_reg1_9~=n_reg1_3(i3) & n_reg1_9~=n_reg1_4(i4) & n_reg1_9~=n_reg1_5(i5) & n_reg1_9~=n_reg1_6(i6) & n_reg1_9~=n_reg1_7(i7) & n_reg1_9~=n_reg1_8(i8));
																				sz9=size(n_reg1_9);
																				for i9=1:sz9(2)
																					if (n_reg1_9(i9) == loc_of_sec_reg)
																						display('A path to reg 8 is 9');
																						weightp9_reg8(t9)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9));
																						t9=t9+1;
																						else
																							n_reg1_10=find(neighborss(n_reg1_9(i9),:));
																							n_reg1_10=n_reg1_10(n_reg1_10~=loc_of_reg1 & n_reg1_10~=n_reg1_2(i2) & n_reg1_10~=n_reg1_3(i3) & n_reg1_10~=n_reg1_4(i4) & n_reg1_10~=n_reg1_5(i5) & n_reg1_10~=n_reg1_6(i6) & n_reg1_10~=n_reg1_7(i7) & n_reg1_10~=n_reg1_8(i8) & n_reg1_10~=n_reg1_9(i9));
																							sz10=size(n_reg1_10);
																							for i10=1:sz10(2)
																								if (n_reg1_10(i10) == loc_of_sec_reg)
																									display('A path to reg 8 is 10');
																									weightp10_reg8(t10)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10));
																									t10=t10+1;
																									else
																										n_reg1_11=find(neighborss(n_reg1_10(i10),:));
																										n_reg1_11=n_reg1_11(n_reg1_11~=loc_of_reg1 & n_reg1_11~=n_reg1_2(i2) & n_reg1_11~=n_reg1_3(i3) & n_reg1_11~=n_reg1_4(i4) & n_reg1_11~=n_reg1_5(i5) & n_reg1_11~=n_reg1_6(i6) & n_reg1_11~=n_reg1_7(i7) & n_reg1_11~=n_reg1_8(i8) & n_reg1_11~=n_reg1_9(i9) & n_reg1_11~=n_reg1_10(i10));
																										sz11=size(n_reg1_11);
																										for i11=1:sz11(2)
																											if (n_reg1_11(i11) == loc_of_sec_reg)
																												display('A path to reg 8 is 11');
																												weightp11_reg8(t11)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11));
																												t11=t11+1;
																												else
																													n_reg1_12=find(neighborss(n_reg1_11(i11),:));
																													n_reg1_12=n_reg1_12(n_reg1_12~=loc_of_reg1 & n_reg1_12~=n_reg1_2(i2) & n_reg1_12~=n_reg1_3(i3) & n_reg1_12~=n_reg1_4(i4) & n_reg1_12~=n_reg1_5(i5) & n_reg1_12~=n_reg1_6(i6) & n_reg1_12~=n_reg1_7(i7) & n_reg1_12~=n_reg1_8(i8) & n_reg1_12~=n_reg1_9(i9) & n_reg1_12~=n_reg1_10(i10) & n_reg1_12~=n_reg1_11(i11));
																													sz12=size(n_reg1_12);
																													for i12=1:sz12(2)
																														if (n_reg1_12(i12) == loc_of_sec_reg)
																															display('A path to reg 8 is 12');
																															weightp12_reg8(t12)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12));
																															t12=t12+1;
																															else
																																n_reg1_13=find(neighborss(n_reg1_12(i12),:));
																																n_reg1_13=n_reg1_13(n_reg1_13~=loc_of_reg1 & n_reg1_13~=n_reg1_2(i2) & n_reg1_13~=n_reg1_3(i3) & n_reg1_13~=n_reg1_4(i4) & n_reg1_13~=n_reg1_5(i5) & n_reg1_13~=n_reg1_6(i6) & n_reg1_13~=n_reg1_7(i7) & n_reg1_13~=n_reg1_8(i8) & n_reg1_13~=n_reg1_9(i9) & n_reg1_13~=n_reg1_10(i10) & n_reg1_13~=n_reg1_11(i11) & n_reg1_13~=n_reg1_12(i12));
																																sz13=size(n_reg1_13);
																																for i13=1:sz13(2)
																																	if (n_reg1_13(i13) == loc_of_sec_reg)
																																		display('A path to reg 8 is 13');
																																		weightp13_reg8(t13)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n_reg1_14=find(neighborss(n_reg1_13(i13),:));
																																		% 	n_reg1_14=n_reg1_14(n_reg1_14~=loc_of_reg1 & n_reg1_14~=n_reg1_2(i2) & n_reg1_14~=n_reg1_3(i3) & n_reg1_14~=n_reg1_4(i4) & n_reg1_14~=n_reg1_5(i5) & n_reg1_14~=n_reg1_6(i6) & n_reg1_14~=n_reg1_7(i7) & n_reg1_14~=n_reg1_8(i8) & n_reg1_14~=n_reg1_9(i9) & n_reg1_14~=n_reg1_10(i10) & n_reg1_14~=n_reg1_11(i11) & n_reg1_14~=n_reg1_12(i12) & n_reg1_14~=n_reg1_13(i13));
																																		% 	sz14=size(n_reg1_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n_reg1_14(i14) == loc_of_sec_reg)
																																		% 			display('A path to reg 8 is 14');
																																		% 			weightp14_reg8(t14)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n_reg1_15=find(neighborss(n_reg1_14(i14),:));
																												 					% 					n_reg1_15=n_reg1_15(n_reg1_15~=loc_of_reg1 & n_reg1_15~=n_reg1_2(i2) & n_reg1_15~=n_reg1_3(i3) & n_reg1_15~=n_reg1_4(i4) & n_reg1_15~=n_reg1_5(i5) & n_reg1_15~=n_reg1_6(i6) & n_reg1_15~=n_reg1_7(i7) & n_reg1_15~=n_reg1_8(i8) & n_reg1_15~=n_reg1_9(i9) & n_reg1_15~=n_reg1_10(i10) & n_reg1_15~=n_reg1_11(i11) & n_reg1_15~=n_reg1_12(i12) & n_reg1_15~=n_reg1_13(i13) & n_reg1_15~=n_reg1_14(i14));
																												 					% 					sz15=size(n_reg1_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n_reg1_15(i15) == loc_of_sec_reg)
																												 					% 							display('A path to reg 8 is 15');
																												 					% 							weightp15_reg8(t15)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14))*lev(n_reg1_14(i14),n_reg1_15(i15));
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
loc_of_sec_reg = loc_of_sec_reg_arr(8);

    	for i1=1:sz(2)
    		if (n_reg1(i1) == loc_of_sec_reg)
    			display('A path to reg 9 is 1');
    			weightp1_reg9(t1)=lev(loc_of_reg1,n_reg1(i1));
    			t1=t1+1;
    		else
    			n_reg1_2=find(neighborss(n_reg1(i1),:));
    			n_reg1_2=n_reg1_2(n_reg1_2~=loc_of_reg1);
    			sz2=size(n_reg1_2);
    			for i2=1:sz2(2)
    				if (n_reg1_2(i2) == loc_of_sec_reg)
    					display('A path to reg 9 is 2');
    					weightp2_reg9(t2)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2));
    					t2=t2+1;
    				else
    					n_reg1_3=find(neighborss(n_reg1_2(i2),:));
    					n_reg1_3=n_reg1_3(n_reg1_3~=n_reg1_2(i2) & n_reg1_3~=loc_of_reg1);
    					sz3=size(n_reg1_3);
    					for i3=1:sz3(2)
    						if (n_reg1_3(i3) == loc_of_sec_reg)
    							display('A path to reg 9 is 3');
    							weightp3_reg9(t3)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3));
    							t3=t3+1;
    						else
    							n_reg1_4=find(neighborss(n_reg1_3(i3),:));
    							n_reg1_4=n_reg1_4(n_reg1_4~=n_reg1_2(i2) & n_reg1_4~=n_reg1_3(i3) & n_reg1_4~=loc_of_reg1);
    							sz4=size(n_reg1_4);
    							for i4=1:sz4(2)
    								if (n_reg1_4(i4) == loc_of_sec_reg)
    								display('A path to reg 9 is 4');
    								weightp4_reg9(t4)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4));
    								t4=t4+1;
        							else
    									n_reg1_5=find(neighborss(n_reg1_4(i4),:));
    									n_reg1_5=n_reg1_5(n_reg1_5~=n_reg1_2(i2) & n_reg1_5~=n_reg1_3(i3) & n_reg1_5~=n_reg1_4(i4) & n_reg1_5~=loc_of_reg1);
    									sz5=size(n_reg1_5);
    									for i5=1:sz5(2)
    										if (n_reg1_5(i5) == loc_of_sec_reg)
    										display('A path to reg 9 is 5');
    										weightp5_reg9(t5)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5));
    										t5=t5+1;
    										else
    											n_reg1_6=find(neighborss(n_reg1_5(i5),:));
    											n_reg1_6=n_reg1_6(n_reg1_6~=loc_of_reg1 & n_reg1_6~=n_reg1_2(i2) & n_reg1_6~=n_reg1_3(i3) & n_reg1_6~=n_reg1_4(i4) & n_reg1_6~=n_reg1_5(i5));
    											sz6=size(n_reg1_6);
    											for i6=1:sz6(2)
    												if (n_reg1_6(i6) == loc_of_sec_reg)
    													display('A path to reg 9 is 6');
    													weightp6_reg9(t6)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6));
    													t6=t6+1;
    													else
															n_reg1_7=find(neighborss(n_reg1_6(i6),:));
															n_reg1_7=n_reg1_7(n_reg1_7~=loc_of_reg1 & n_reg1_7~=n_reg1_2(i2) & n_reg1_7~=n_reg1_3(i3) & n_reg1_7~=n_reg1_4(i4) & n_reg1_7~=n_reg1_5(i5) & n_reg1_7~=n_reg1_6(i6));
															sz7=size(n_reg1_7);
															for i7=1:sz7(2)
																if (n_reg1_7(i7) == loc_of_sec_reg)
																	display('A path to reg 9 is 7');
																	weightp7_reg9(t7)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7));
																	t7=t7+1;
																	else
																		n_reg1_8=find(neighborss(n_reg1_7(i7),:));
																		n_reg1_8=n_reg1_8(n_reg1_8~=loc_of_reg1 & n_reg1_8~=n_reg1_2(i2) & n_reg1_8~=n_reg1_3(i3) & n_reg1_8~=n_reg1_4(i4) & n_reg1_8~=n_reg1_5(i5) & n_reg1_8~=n_reg1_6(i6) & n_reg1_8~=n_reg1_7(i7));
																		sz8=size(n_reg1_8);
																		for i8=1:sz8(2)
																			if (n_reg1_8(i8) == loc_of_sec_reg)
																			display('A path to reg 9 is 8');
																			weightp8_reg9(t8)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8));
																			t8=t8+1;
																			else
																				n_reg1_9=find(neighborss(n_reg1_8(i8),:));
																				n_reg1_9=n_reg1_9(n_reg1_9~=loc_of_reg1 & n_reg1_9~=n_reg1_2(i2) & n_reg1_9~=n_reg1_3(i3) & n_reg1_9~=n_reg1_4(i4) & n_reg1_9~=n_reg1_5(i5) & n_reg1_9~=n_reg1_6(i6) & n_reg1_9~=n_reg1_7(i7) & n_reg1_9~=n_reg1_8(i8));
																				sz9=size(n_reg1_9);
																				for i9=1:sz9(2)
																					if (n_reg1_9(i9) == loc_of_sec_reg)
																						display('A path to reg 9 is 9');
																						weightp9_reg9(t9)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9));
																						t9=t9+1;
																						else
																							n_reg1_10=find(neighborss(n_reg1_9(i9),:));
																							n_reg1_10=n_reg1_10(n_reg1_10~=loc_of_reg1 & n_reg1_10~=n_reg1_2(i2) & n_reg1_10~=n_reg1_3(i3) & n_reg1_10~=n_reg1_4(i4) & n_reg1_10~=n_reg1_5(i5) & n_reg1_10~=n_reg1_6(i6) & n_reg1_10~=n_reg1_7(i7) & n_reg1_10~=n_reg1_8(i8) & n_reg1_10~=n_reg1_9(i9));
																							sz10=size(n_reg1_10);
																							for i10=1:sz10(2)
																								if (n_reg1_10(i10) == loc_of_sec_reg)
																									display('A path to reg 9 is 10');
																									weightp10_reg9(t10)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10));
																									t10=t10+1;
																									else
																										n_reg1_11=find(neighborss(n_reg1_10(i10),:));
																										n_reg1_11=n_reg1_11(n_reg1_11~=loc_of_reg1 & n_reg1_11~=n_reg1_2(i2) & n_reg1_11~=n_reg1_3(i3) & n_reg1_11~=n_reg1_4(i4) & n_reg1_11~=n_reg1_5(i5) & n_reg1_11~=n_reg1_6(i6) & n_reg1_11~=n_reg1_7(i7) & n_reg1_11~=n_reg1_8(i8) & n_reg1_11~=n_reg1_9(i9) & n_reg1_11~=n_reg1_10(i10));
																										sz11=size(n_reg1_11);
																										for i11=1:sz11(2)
																											if (n_reg1_11(i11) == loc_of_sec_reg)
																												display('A path to reg 9 is 11');
																												weightp11_reg9(t11)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11));
																												t11=t11+1;
																												else
																													n_reg1_12=find(neighborss(n_reg1_11(i11),:));
																													n_reg1_12=n_reg1_12(n_reg1_12~=loc_of_reg1 & n_reg1_12~=n_reg1_2(i2) & n_reg1_12~=n_reg1_3(i3) & n_reg1_12~=n_reg1_4(i4) & n_reg1_12~=n_reg1_5(i5) & n_reg1_12~=n_reg1_6(i6) & n_reg1_12~=n_reg1_7(i7) & n_reg1_12~=n_reg1_8(i8) & n_reg1_12~=n_reg1_9(i9) & n_reg1_12~=n_reg1_10(i10) & n_reg1_12~=n_reg1_11(i11));
																													sz12=size(n_reg1_12);
																													for i12=1:sz12(2)
																														if (n_reg1_12(i12) == loc_of_sec_reg)
																															display('A path to reg 9 is 12');
																															weightp12_reg9(t12)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12));
																															t12=t12+1;
																															else
																																n_reg1_13=find(neighborss(n_reg1_12(i12),:));
																																n_reg1_13=n_reg1_13(n_reg1_13~=loc_of_reg1 & n_reg1_13~=n_reg1_2(i2) & n_reg1_13~=n_reg1_3(i3) & n_reg1_13~=n_reg1_4(i4) & n_reg1_13~=n_reg1_5(i5) & n_reg1_13~=n_reg1_6(i6) & n_reg1_13~=n_reg1_7(i7) & n_reg1_13~=n_reg1_8(i8) & n_reg1_13~=n_reg1_9(i9) & n_reg1_13~=n_reg1_10(i10) & n_reg1_13~=n_reg1_11(i11) & n_reg1_13~=n_reg1_12(i12));
																																sz13=size(n_reg1_13);
																																for i13=1:sz13(2)
																																	if (n_reg1_13(i13) == loc_of_sec_reg)
																																		display('A path to reg 9 is 13');
																																		weightp13_reg9(t13)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n_reg1_14=find(neighborss(n_reg1_13(i13),:));
																																		% 	n_reg1_14=n_reg1_14(n_reg1_14~=loc_of_reg1 & n_reg1_14~=n_reg1_2(i2) & n_reg1_14~=n_reg1_3(i3) & n_reg1_14~=n_reg1_4(i4) & n_reg1_14~=n_reg1_5(i5) & n_reg1_14~=n_reg1_6(i6) & n_reg1_14~=n_reg1_7(i7) & n_reg1_14~=n_reg1_8(i8) & n_reg1_14~=n_reg1_9(i9) & n_reg1_14~=n_reg1_10(i10) & n_reg1_14~=n_reg1_11(i11) & n_reg1_14~=n_reg1_12(i12) & n_reg1_14~=n_reg1_13(i13));
																																		% 	sz14=size(n_reg1_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n_reg1_14(i14) == loc_of_sec_reg)
																																		% 			display('A path to reg 9 is 14');
																																		% 			weightp14_reg9(t14)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n_reg1_15=find(neighborss(n_reg1_14(i14),:));
																												 					% 					n_reg1_15=n_reg1_15(n_reg1_15~=loc_of_reg1 & n_reg1_15~=n_reg1_2(i2) & n_reg1_15~=n_reg1_3(i3) & n_reg1_15~=n_reg1_4(i4) & n_reg1_15~=n_reg1_5(i5) & n_reg1_15~=n_reg1_6(i6) & n_reg1_15~=n_reg1_7(i7) & n_reg1_15~=n_reg1_8(i8) & n_reg1_15~=n_reg1_9(i9) & n_reg1_15~=n_reg1_10(i10) & n_reg1_15~=n_reg1_11(i11) & n_reg1_15~=n_reg1_12(i12) & n_reg1_15~=n_reg1_13(i13) & n_reg1_15~=n_reg1_14(i14));
																												 					% 					sz15=size(n_reg1_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n_reg1_15(i15) == loc_of_sec_reg)
																												 					% 							display('A path to reg 9 is 15');
																												 					% 							weightp15_reg9(t15)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14))*lev(n_reg1_14(i14),n_reg1_15(i15));
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
loc_of_sec_reg = loc_of_sec_reg_arr(9);

    	for i1=1:sz(2)
    		if (n_reg1(i1) == loc_of_sec_reg)
    			display('A path to reg 10 is 1');
    			weightp1_reg10(t1)=lev(loc_of_reg1,n_reg1(i1));
    			t1=t1+1;
    		else
    			n_reg1_2=find(neighborss(n_reg1(i1),:));
    			n_reg1_2=n_reg1_2(n_reg1_2~=loc_of_reg1);
    			sz2=size(n_reg1_2);
    			for i2=1:sz2(2)
    				if (n_reg1_2(i2) == loc_of_sec_reg)
    					display('A path to reg 10 is 2');
    					weightp2_reg10(t2)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2));
    					t2=t2+1;
    				else
    					n_reg1_3=find(neighborss(n_reg1_2(i2),:));
    					n_reg1_3=n_reg1_3(n_reg1_3~=n_reg1_2(i2) & n_reg1_3~=loc_of_reg1);
    					sz3=size(n_reg1_3);
    					for i3=1:sz3(2)
    						if (n_reg1_3(i3) == loc_of_sec_reg)
    							display('A path to reg 10 is 3');
    							weightp3_reg10(t3)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3));
    							t3=t3+1;
    						else
    							n_reg1_4=find(neighborss(n_reg1_3(i3),:));
    							n_reg1_4=n_reg1_4(n_reg1_4~=n_reg1_2(i2) & n_reg1_4~=n_reg1_3(i3) & n_reg1_4~=loc_of_reg1);
    							sz4=size(n_reg1_4);
    							for i4=1:sz4(2)
    								if (n_reg1_4(i4) == loc_of_sec_reg)
    								display('A path to reg 10 is 4');
    								weightp4_reg10(t4)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4));
    								t4=t4+1;
        							else
    									n_reg1_5=find(neighborss(n_reg1_4(i4),:));
    									n_reg1_5=n_reg1_5(n_reg1_5~=n_reg1_2(i2) & n_reg1_5~=n_reg1_3(i3) & n_reg1_5~=n_reg1_4(i4) & n_reg1_5~=loc_of_reg1);
    									sz5=size(n_reg1_5);
    									for i5=1:sz5(2)
    										if (n_reg1_5(i5) == loc_of_sec_reg)
    										display('A path to reg 10 is 5');
    										weightp5_reg10(t5)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5));
    										t5=t5+1;
    										else
    											n_reg1_6=find(neighborss(n_reg1_5(i5),:));
    											n_reg1_6=n_reg1_6(n_reg1_6~=loc_of_reg1 & n_reg1_6~=n_reg1_2(i2) & n_reg1_6~=n_reg1_3(i3) & n_reg1_6~=n_reg1_4(i4) & n_reg1_6~=n_reg1_5(i5));
    											sz6=size(n_reg1_6);
    											for i6=1:sz6(2)
    												if (n_reg1_6(i6) == loc_of_sec_reg)
    													display('A path to reg 10 is 6');
    													weightp6_reg10(t6)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6));
    													t6=t6+1;
    													else
															n_reg1_7=find(neighborss(n_reg1_6(i6),:));
															n_reg1_7=n_reg1_7(n_reg1_7~=loc_of_reg1 & n_reg1_7~=n_reg1_2(i2) & n_reg1_7~=n_reg1_3(i3) & n_reg1_7~=n_reg1_4(i4) & n_reg1_7~=n_reg1_5(i5) & n_reg1_7~=n_reg1_6(i6));
															sz7=size(n_reg1_7);
															for i7=1:sz7(2)
																if (n_reg1_7(i7) == loc_of_sec_reg)
																	display('A path to reg 10 is 7');
																	weightp7_reg10(t7)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7));
																	t7=t7+1;
																	else
																		n_reg1_8=find(neighborss(n_reg1_7(i7),:));
																		n_reg1_8=n_reg1_8(n_reg1_8~=loc_of_reg1 & n_reg1_8~=n_reg1_2(i2) & n_reg1_8~=n_reg1_3(i3) & n_reg1_8~=n_reg1_4(i4) & n_reg1_8~=n_reg1_5(i5) & n_reg1_8~=n_reg1_6(i6) & n_reg1_8~=n_reg1_7(i7));
																		sz8=size(n_reg1_8);
																		for i8=1:sz8(2)
																			if (n_reg1_8(i8) == loc_of_sec_reg)
																			display('A path to reg 10 is 8');
																			weightp8_reg10(t8)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8));
																			t8=t8+1;
																			else
																				n_reg1_9=find(neighborss(n_reg1_8(i8),:));
																				n_reg1_9=n_reg1_9(n_reg1_9~=loc_of_reg1 & n_reg1_9~=n_reg1_2(i2) & n_reg1_9~=n_reg1_3(i3) & n_reg1_9~=n_reg1_4(i4) & n_reg1_9~=n_reg1_5(i5) & n_reg1_9~=n_reg1_6(i6) & n_reg1_9~=n_reg1_7(i7) & n_reg1_9~=n_reg1_8(i8));
																				sz9=size(n_reg1_9);
																				for i9=1:sz9(2)
																					if (n_reg1_9(i9) == loc_of_sec_reg)
																						display('A path to reg 10 is 9');
																						weightp9_reg10(t9)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9));
																						t9=t9+1;
																						else
																							n_reg1_10=find(neighborss(n_reg1_9(i9),:));
																							n_reg1_10=n_reg1_10(n_reg1_10~=loc_of_reg1 & n_reg1_10~=n_reg1_2(i2) & n_reg1_10~=n_reg1_3(i3) & n_reg1_10~=n_reg1_4(i4) & n_reg1_10~=n_reg1_5(i5) & n_reg1_10~=n_reg1_6(i6) & n_reg1_10~=n_reg1_7(i7) & n_reg1_10~=n_reg1_8(i8) & n_reg1_10~=n_reg1_9(i9));
																							sz10=size(n_reg1_10);
																							for i10=1:sz10(2)
																								if (n_reg1_10(i10) == loc_of_sec_reg)
																									display('A path to reg 10 is 10');
																									weightp10_reg10(t10)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10));
																									t10=t10+1;
																									else
																										n_reg1_11=find(neighborss(n_reg1_10(i10),:));
																										n_reg1_11=n_reg1_11(n_reg1_11~=loc_of_reg1 & n_reg1_11~=n_reg1_2(i2) & n_reg1_11~=n_reg1_3(i3) & n_reg1_11~=n_reg1_4(i4) & n_reg1_11~=n_reg1_5(i5) & n_reg1_11~=n_reg1_6(i6) & n_reg1_11~=n_reg1_7(i7) & n_reg1_11~=n_reg1_8(i8) & n_reg1_11~=n_reg1_9(i9) & n_reg1_11~=n_reg1_10(i10));
																										sz11=size(n_reg1_11);
																										for i11=1:sz11(2)
																											if (n_reg1_11(i11) == loc_of_sec_reg)
																												display('A path to reg 10 is 11');
																												weightp11_reg10(t11)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11));
																												t11=t11+1;
																												else
																													n_reg1_12=find(neighborss(n_reg1_11(i11),:));
																													n_reg1_12=n_reg1_12(n_reg1_12~=loc_of_reg1 & n_reg1_12~=n_reg1_2(i2) & n_reg1_12~=n_reg1_3(i3) & n_reg1_12~=n_reg1_4(i4) & n_reg1_12~=n_reg1_5(i5) & n_reg1_12~=n_reg1_6(i6) & n_reg1_12~=n_reg1_7(i7) & n_reg1_12~=n_reg1_8(i8) & n_reg1_12~=n_reg1_9(i9) & n_reg1_12~=n_reg1_10(i10) & n_reg1_12~=n_reg1_11(i11));
																													sz12=size(n_reg1_12);
																													for i12=1:sz12(2)
																														if (n_reg1_12(i12) == loc_of_sec_reg)
																															display('A path to reg 10 is 12');
																															weightp12_reg10(t12)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12));
																															t12=t12+1;
																															else
																																n_reg1_13=find(neighborss(n_reg1_12(i12),:));
																																n_reg1_13=n_reg1_13(n_reg1_13~=loc_of_reg1 & n_reg1_13~=n_reg1_2(i2) & n_reg1_13~=n_reg1_3(i3) & n_reg1_13~=n_reg1_4(i4) & n_reg1_13~=n_reg1_5(i5) & n_reg1_13~=n_reg1_6(i6) & n_reg1_13~=n_reg1_7(i7) & n_reg1_13~=n_reg1_8(i8) & n_reg1_13~=n_reg1_9(i9) & n_reg1_13~=n_reg1_10(i10) & n_reg1_13~=n_reg1_11(i11) & n_reg1_13~=n_reg1_12(i12));
																																sz13=size(n_reg1_13);
																																for i13=1:sz13(2)
																																	if (n_reg1_13(i13) == loc_of_sec_reg)
																																		display('A path to reg 10 is 13');
																																		weightp13_reg10(t13)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13));
																																		t13=t13+1;
																																		% else
																																		% 	n_reg1_14=find(neighborss(n_reg1_13(i13),:));
																																		% 	n_reg1_14=n_reg1_14(n_reg1_14~=loc_of_reg1 & n_reg1_14~=n_reg1_2(i2) & n_reg1_14~=n_reg1_3(i3) & n_reg1_14~=n_reg1_4(i4) & n_reg1_14~=n_reg1_5(i5) & n_reg1_14~=n_reg1_6(i6) & n_reg1_14~=n_reg1_7(i7) & n_reg1_14~=n_reg1_8(i8) & n_reg1_14~=n_reg1_9(i9) & n_reg1_14~=n_reg1_10(i10) & n_reg1_14~=n_reg1_11(i11) & n_reg1_14~=n_reg1_12(i12) & n_reg1_14~=n_reg1_13(i13));
																																		% 	sz14=size(n_reg1_14);
																																		% 	for i14=1:sz14(2)
																																		% 		if (n_reg1_14(i14) == loc_of_sec_reg)
																																		% 			display('A path to reg 10 is 14');
																																		% 			weightp14_reg10(t14)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14));
																																		% 			t14=t14+1;
																												 					% 				else
																												 					% 					n_reg1_15=find(neighborss(n_reg1_14(i14),:));
																												 					% 					n_reg1_15=n_reg1_15(n_reg1_15~=loc_of_reg1 & n_reg1_15~=n_reg1_2(i2) & n_reg1_15~=n_reg1_3(i3) & n_reg1_15~=n_reg1_4(i4) & n_reg1_15~=n_reg1_5(i5) & n_reg1_15~=n_reg1_6(i6) & n_reg1_15~=n_reg1_7(i7) & n_reg1_15~=n_reg1_8(i8) & n_reg1_15~=n_reg1_9(i9) & n_reg1_15~=n_reg1_10(i10) & n_reg1_15~=n_reg1_11(i11) & n_reg1_15~=n_reg1_12(i12) & n_reg1_15~=n_reg1_13(i13) & n_reg1_15~=n_reg1_14(i14));
																												 					% 					sz15=size(n_reg1_15);
																												 					% 					for i15=1:sz15(2)
																												 					% 						if (n_reg1_15(i15) == loc_of_sec_reg)
																												 					% 							display('A path to reg 10 is 15');
																												 					% 							weightp15_reg10(t15)=lev(loc_of_reg1,n_reg1(i1))*lev(n_reg1(i1),n_reg1_2(i2))*lev(n_reg1_2(i2),n_reg1_3(i3))*lev(n_reg1_3(i3),n_reg1_4(i4))*lev(n_reg1_4(i4),n_reg1_5(i5))*lev(n_reg1_5(i5),n_reg1_6(i6))*lev(n_reg1_6(i6),n_reg1_7(i7))*lev(n_reg1_7(i7),n_reg1_8(i8))*lev(n_reg1_8(i8),n_reg1_9(i9))*lev(n_reg1_9(i9),n_reg1_10(i10))*lev(n_reg1_10(i10),n_reg1_11(i11))*lev(n_reg1_11(i11),n_reg1_12(i12))*lev(n_reg1_12(i12),n_reg1_13(i13))*lev(n_reg1_13(i13),n_reg1_14(i14))*lev(n_reg1_14(i14),n_reg1_15(i15));
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
    
% weightp1_reg2_m=mean(weightp1_reg2);
% weightp2_reg2_m=mean(weightp2_reg2/1e2);
%weightp3_reg2_m=mean(weightp3_reg2/1e4);
weightp4_reg2_m=mean(weightp4_reg2/1e6);
weightp5_reg2_m=mean(weightp5_reg2/1e8);
weightp6_reg2_m=mean(weightp6_reg2/1e10);
weightp7_reg2_m=mean(weightp7_reg2/1e12);
weightp8_reg2_m=mean(weightp8_reg2/1e14);
weightp9_reg2_m=mean(weightp9_reg2/1e16);
weightp10_reg2_m=mean(weightp10_reg2/1e18);
weightp11_reg2_m=mean(weightp11_reg2/1e20);
weightp12_reg2_m=mean(weightp12_reg2/1e22);
weightp13_reg2_m=mean(weightp13_reg2/1e24);

 %to_reg_2_w=mean([weightp1_reg2_m/1 weightp2_reg2_m/2 weightp3_reg2_m/3 weightp4_reg2_m/4 weightp5_reg2_m/5 weightp6_reg2_m/6 weightp7_reg2_m/7 weightp8_reg2_m/8 weightp9_reg2_m/9 weightp10_reg2_m/10 weightp11_reg2_m/11 weightp12_reg2_m/12 weightp13_reg2_m/13]);

%weightp1_reg3_m=mean(weightp1_reg3);
%weightp2_reg3_m=mean(weightp2_reg3/1e2);
%weightp3_reg3_m=mean(weightp3_reg3/1e4);
% weightp4_reg3_m=mean(weightp4_reg3/1e6);
% weightp5_reg3_m=mean(weightp5_reg3/1e8);
% weightp6_reg3_m=mean(weightp6_reg3/1e10);
%weightp7_reg3_m=mean(weightp7_reg3/1e12);
weightp8_reg3_m=mean(weightp8_reg3/1e14);
weightp9_reg3_m=mean(weightp9_reg3/1e16);
weightp10_reg3_m=mean(weightp10_reg3/1e18);
weightp11_reg3_m=mean(weightp11_reg3/1e20);
weightp12_reg3_m=mean(weightp12_reg3/1e22);
weightp13_reg3_m=mean(weightp13_reg3/1e24);
 %to_reg_3_w=mean([weightp7_reg3_m/7 weightp8_reg3_m/8 weightp9_reg3_m/9 weightp10_reg3_m/10 weightp11_reg3_m/11 weightp12_reg3_m/12 weightp13_reg3_m/13]);

% % weightp1_reg4_m=mean(weightp1_reg4);
% % weightp2_reg4_m=mean(weightp2_reg4/1e2);
% % weightp3_reg4_m=mean(weightp3_reg4/1e4);
% weightp4_reg4_m=mean(weightp4_reg4/1e6);
% weightp5_reg4_m=mean(weightp5_reg4/1e8);
%weightp6_reg4_m=mean(weightp6_reg4/1e10);
%weightp7_reg4_m=mean(weightp7_reg4/1e12);
%weightp8_reg4_m=mean(weightp8_reg4/1e14);
weightp9_reg4_m=mean(weightp9_reg4/1e16);
weightp10_reg4_m=mean(weightp10_reg4/1e18);
weightp11_reg4_m=mean(weightp11_reg4/1e20);
weightp12_reg4_m=mean(weightp12_reg4/1e22);
weightp13_reg4_m=mean(weightp13_reg4/1e24);
% to_reg_4_w=mean([weightp4_reg4_m/4 weightp5_reg4_m/5 weightp6_reg4_m/6 weightp7_reg4_m/7 weightp8_reg4_m/8 weightp9_reg4_m/9 weightp10_reg4_m/10 weightp11_reg4_m/11 weightp12_reg4_m/12 weightp13_reg4_m/13]);


% % weightp1_reg5_m=mean(weightp1_reg5);
% % weightp2_reg5_m=mean(weightp2_reg5/1e2);
% % weightp3_reg5_m=mean(weightp3_reg5/1e4);
% % weightp4_reg5_m=mean(weightp4_reg5/1e6);
% % weightp5_reg5_m=mean(weightp5_reg5/1e8);
% % weightp6_reg5_m=mean(weightp6_reg5/1e10);
%  weightp7_reg5_m=mean(weightp7_reg5/1e12);
% weightp8_reg5_m=mean(weightp8_reg5/1e14);
% weightp9_reg5_m=mean(weightp9_reg5/1e16);
% weightp10_reg5_m=mean(weightp10_reg5/1e18);
%weightp11_reg5_m=mean(weightp11_reg5/1e20);
weightp12_reg5_m=mean(weightp12_reg5/1e22);
weightp13_reg5_m=mean(weightp13_reg5/1e24);
% to_reg_5_w=mean([ weightp8_reg5_m/8 weightp9_reg5_m/9 weightp10_reg5_m/10 weightp11_reg5_m/11 weightp12_reg5_m/12 weightp13_reg5_m/13]);


% % weightp1_reg6_m=mean(weightp1_reg6);
% % weightp2_reg6_m=mean(weightp2_reg6/1e2);
% % weightp3_reg6_m=mean(weightp3_reg6/1e4);
% % weightp4_reg6_m=mean(weightp4_reg6/1e6);
% % weightp5_reg6_m=mean(weightp5_reg6/1e8);
% weightp6_reg6_m=mean(weightp6_reg6/1e10);
%weightp7_reg6_m=mean(weightp7_reg6/1e12);
%weightp8_reg6_m=mean(weightp8_reg6/1e14);
%weightp9_reg6_m=mean(weightp9_reg6/1e16);
%weightp10_reg6_m=mean(weightp10_reg6/1e18);
weightp11_reg6_m=mean(weightp11_reg6/1e20);
weightp12_reg6_m=mean(weightp12_reg6/1e22);
weightp13_reg6_m=mean(weightp13_reg6/1e24);
% to_reg_6_w=mean([weightp11_reg6_m/11 weightp12_reg6_m/12 weightp13_reg6_m/13]);

% % weightp1_reg7_m=mean(weightp1_reg7);
% % weightp2_reg7_m=mean(weightp2_reg7/1e2);
% % weightp3_reg7_m=mean(weightp3_reg7/1e4);
% % weightp4_reg7_m=mean(weightp4_reg7/1e6);
% weightp5_reg7_m=mean(weightp5_reg7/1e8);
% weightp6_reg7_m=mean(weightp6_reg7/1e10);
%weightp7_reg7_m=mean(weightp7_reg7/1e12);
weightp8_reg7_m=mean(weightp8_reg7/1e14);
weightp9_reg7_m=mean(weightp9_reg7/1e16);
weightp10_reg7_m=mean(weightp10_reg7/1e18);
weightp11_reg7_m=mean(weightp11_reg7/1e20);
weightp12_reg7_m=mean(weightp12_reg7/1e22);
weightp13_reg7_m=mean(weightp13_reg7/1e24);
% to_reg_7_w=mean([weightp5_reg7_m/5 weightp6_reg7_m/6 weightp7_reg7_m/7 weightp8_reg7_m/8 weightp9_reg7_m/9 weightp10_reg7_m/10 weightp11_reg7_m/11 weightp12_reg7_m/12 weightp13_reg7_m/13]);


% % weightp1_reg8_m=mean(weightp1_reg8);
% % weightp2_reg8_m=mean(weightp2_reg8/1e2);
% % weightp3_reg8_m=mean(weightp3_reg8/1e4);
% % weightp4_reg8_m=mean(weightp4_reg8/1e6);
% % weightp5_reg8_m=mean(weightp5_reg8/1e8);
% % weightp6_reg8_m=mean(weightp6_reg8/1e10);
% % weightp7_reg8_m=mean(weightp7_reg8/1e12);
% % weightp8_reg8_m=mean(weightp8_reg8/1e14);
% weightp9_reg8_m=mean(weightp9_reg8/1e16);
% weightp10_reg8_m=mean(weightp10_reg8/1e18);
% weightp11_reg8_m=mean(weightp11_reg8/1e20);
%weightp12_reg8_m=mean(weightp12_reg8/1e22);
weightp13_reg8_m=mean(weightp13_reg8/1e24);
 %to_reg_8_w=mean([weightp9_reg8_m/9 weightp10_reg8_m/10 weightp11_reg8_m/11 weightp12_reg8_m/12 weightp13_reg8_m/13]);


% % weightp1_reg9_m=mean(weightp1_reg9);
% % weightp2_reg9_m=mean(weightp2_reg9/1e2);
% weightp3_reg9_m=mean(weightp3_reg9/1e4);
%weightp4_reg9_m=mean(weightp4_reg9/1e6);
%weightp5_reg9_m=mean(weightp5_reg9/1e8);
%weightp6_reg9_m=mean(weightp6_reg9/1e10);
weightp7_reg9_m=mean(weightp7_reg9/1e12);
weightp8_reg9_m=mean(weightp8_reg9/1e14);
weightp9_reg9_m=mean(weightp9_reg9/1e16);
weightp10_reg9_m=mean(weightp10_reg9/1e18);
weightp11_reg9_m=mean(weightp11_reg9/1e20);
weightp12_reg9_m=mean(weightp12_reg9/1e22);
weightp13_reg9_m=mean(weightp13_reg9/1e24);
%to_reg_9_w=mean([weightp4_reg9_m/4 weightp5_reg9_m/5 weightp6_reg9_m/6 weightp7_reg9_m/7 weightp8_reg9_m/8 weightp9_reg9_m/9 weightp10_reg9_m/10 weightp11_reg9_m/11 weightp12_reg9_m/12 weightp13_reg9_m/13]);

  
% % weightp1_reg10_m=mean(weightp1_reg10);
% % weightp2_reg10_m=mean(weightp2_reg10/1e2);
% % weightp3_reg10_m=mean(weightp3_reg10/1e4);
% % weightp4_reg10_m=mean(weightp4_reg10/1e6);
% % weightp5_reg10_m=mean(weightp5_reg10/1e8);
% % weightp6_reg10_m=mean(weightp6_reg10/1e10);
% % weightp7_reg10_m=mean(weightp7_reg10/1e12);
%weightp8_reg10_m=mean(weightp8_reg10/1e14);
weightp9_reg10_m=mean(weightp9_reg10/1e16);
weightp10_reg10_m=mean(weightp10_reg10/1e18);
weightp11_reg10_m=mean(weightp11_reg10/1e20);
weightp12_reg10_m=mean(weightp12_reg10/1e22);
weightp13_reg10_m=mean(weightp13_reg10/1e24);
%to_reg_10_w=mean([weightp8_reg10_m/8 weightp9_reg10_m/9 weightp10_reg10_m/10 weightp11_reg10_m/11 weightp12_reg10_m/12 weightp13_reg10_m/13]);

 %all_reg_w=[to_reg_2_w to_reg_3_w to_reg_4_w to_reg_5_w to_reg_6_w to_reg_7_w to_reg_8_w to_reg_9_w to_reg_10_w];
 all_reg_shortest_path_withouth_sfr=[weightp4_reg2_m weightp8_reg3_m weightp9_reg4_m weightp12_reg5_m weightp11_reg6_m weightp8_reg7_m weightp13_reg8_m weightp7_reg9_m weightp9_reg10_m];
all_reg_shortest_path_withouth_sfr2=[weightp4_reg2_m/4 weightp8_reg3_m/8 weightp9_reg4_m/9 weightp12_reg5_m/12 weightp11_reg6_m/11 weightp8_reg7_m/8 weightp13_reg8_m/13 weightp7_reg9_m/7 weightp9_reg10_m/9]

