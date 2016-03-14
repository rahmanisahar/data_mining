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

% count=2;
% for i=1:17
%     i_name=int2str(i);
%     m1=0;
%     count=count+2;
%     figure(count)
%     for ind=18:33
%         m1=m1+1;
%         subplot(4,4,m1)
%         for h1=n_1:-1:1
%             for h2=1:1:n_2
%                 check_s=CAT_1{h1,h2};
%                 size_ch=size(check_s);
%                 if (size_ch(2) > 0) 
%                     params=CAT_1{h1,h2};
%                     s=scatter(params(ind,:),params(i,:));
%                     s.MarkerFaceColor = [0 h1/20 h2/20];
%                     if (size_ch(2) > 2) 
%                         corr_mat(h1,h2)=corr(params(i,:)',params(ind,:)');
%                         if (abs(corr_mat(h1,h2)) > 0.5)
%                             strmax = ['p=',num2str(corr_mat(h1,h2))];
%                              text(mean(params(ind,:)),mean(params(i,:)),strmax,'Color',[h1*h2/4000 h1*5/2000 h2*4/2000]);
%                         end
%                     end
%                     hold on
%                 end
%             end
%         end
%     end
%     name1 = strcat(dir,i_name,'vs_raws_18_to_33_for_',n1st,'by',n2st,'.jpeg');
%     saveas(figure(count),name1,'jpeg')    
%     m1=0;
%     figure(count+1)
%     for ind=34:39
%         m1=m1+1;
%         subplot(2,3,m1)
%         for h1=n_1:-1:1
%             for h2=1:1:n_2
%                 check_s=CAT_1{h1,h2};
%                 size_ch=size(check_s);
%                 if (size_ch(2) > 0) 
%                     params=CAT_1{h1,h2};
%                     s=scatter(params(ind,:),params(i,:));
%                     s.MarkerFaceColor = [0 h1/20 h2/20];
%                     if (size_ch(2) > 2) 
%                         corr_mat(h1,h2)=corr(params(i,:)',params(ind,:)');
%                         if (abs(corr_mat(h1,h2)) > 0.5)
%                             strmax = ['p=',num2str(corr_mat(h1,h2))];
%                             text(mean(params(ind,:)),mean(params(i,:)),strmax,'Color',[h1*h2/4000 h1*5/2000 h2*4/2000]);
%                         end
%                     end
%                     hold on
%                 end
%             end
%         end
%     end
%     name2 = strcat(dir,i_name,'vs_raws_34_to_43_for_',n1st,'by',n2st,'.jpeg');
%     saveas(figure(count+1),name2,'jpeg')
% end
% 
% for i=1:17
%     for ind=18:40
%     i_name=int2str(i);
%     j_name=int2str(ind);
%     corr_mat=zeros(n_1,n_2)+0.000000006;
%         for h1=n_1:-1:1
%             for h2=1:1:n_2
%                 check_s=CAT_1{h1,h2};
%                 size_ch=size(check_s);
%                 if (size_ch(2) > 2) 
%                     params=CAT_1{h1,h2};
%                     corr_mat(h1,h2)=corr(params(i,:)',params(ind,:)');
%                 end
%             end
%         end
%         II= find(abs(corr_mat) >=0.5) ;
%         if (II > 0)
%           name1 = strcat(dir,i_name,'_',j_name,'_',n1st,'_by',n2st,'.csv');   
%           % table1 = cell2table(corr_mat);
%          csvwrite(name1,corr_mat);
%         end
%     end
% end