load('new_IRS_N400_p1.mat');
SNR_list_400_1 = SNR_list(1:6,2);

load('new_IRS_N300_p1.mat');
SNR_list_300_1 = SNR_list(1:7,3);

load('new_IRS_N400_p13.mat');
SNR_list_400_13 = SNR_list(1:6,2);

figure      
plot( 1:6,SNR_list_400_1,'r-d','LineWidth',1.3);
hold on
plot( 1:7,SNR_list_300_1 ,'b-o','LineWidth',1.3);
hold on
% plot( 1:6,SNR_list_400_13,'b->','LineWidth',1.3,'Color',[0 0.4 0]);
% hold on
%axis([1,7,35,70]);
grid on
h=legend('$N=400$','$N=300$','Interpreter','latex','Location','best');
%h=legend('$N=400, P_{max}=1$','$N=300, P_{max}=1$','$N=400, P_{max}=1.3$','Interpreter','latex','Location','best');
set(h,'FontSize',15);
xlabel({'Iteration number'},'FontSize',19);
ylabel('Max-min SNR (dB)','FontSize',19);
set(gca, 'YTick', [35 40 45 50 55 60 65 70 75],'FontSize',13); 









