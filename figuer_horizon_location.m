%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成坐标%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IRS_location_center = [200 0 25]';
AP_location = [0 40 10]';

K = 5; %%用户坐标数
center = [ 185 45 1];
% r = 2.5;
% for user_index = 1: K
%     IU_location(:,user_index) = [center(1)+r*cos(rand(1)*2*pi); center(2)+r*sin(rand(1)*2*pi); 1]';
% end

r = 30;
theta1 = 0:2*pi/(r-1):2*pi;                                    %%分成N个theta
IU_location = [    
    168  180  205  210  164
     25   41   59   42   53
      1    1    1    1    1];  
center = sum(IU_location,2)/K;
center = [185 45 1];
Q1_r = [center(1)+r*cos(theta1); center(2)+r*sin(theta1)]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%画平面图%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(Q1_r(1,1:30),Q1_r(2,1:30),'r--','LineWidth',1.5,'Color',[0.75 0.75 0.75])
hold on
plot(AP_location(1),AP_location(2),'r.','MarkerSize',20)
text(AP_location(1)+10,AP_location(2),['BS'],'color','black','FontSize',15)
hold on
plot(IRS_location_center(1),IRS_location_center(2),'k*','LineWidth',2,'MarkerSize',10)
text(IRS_location_center(1)-10,IRS_location_center(2)+6,['IRS'],'color','black','FontSize',15)
hold on
for i=1:K
    plot( IU_location(1,i),IU_location(2,i),'b.', 'LineWidth',2,'MarkerSize',15)
    text(IU_location(1,i)-13,IU_location(2,i)-3,['{\itu}_',num2str(i)],'color','black','FontSize',15)
    hold on
end
%plot( Q_initial(1),Q_initial(2),'+', 'LineWidth',8,'MarkerSize',8 ,'Color',[230 75 90]/256)
%text(Q_initial(1)+5,Q_initial(2)+70,['Initial position'],'color','black','FontSize',13)
% plot( Q1_r(1,T),Q1_r(2,T),'+', 'LineWidth',8,'MarkerSize',8 ,'Color',[230 75 90]/256)
% text(Q1_r(1,T)-180,Q1_r(2,T)+20,['Destination'],'color','black','FontSize',13)
axis([0,300,0,90])

%legend('UAV operation time minimization','Power efficient','Algorithm 1 in stage one', 'Location','best','FontSize',13)
grid on
toc;