load('new_IRS_N400_p13.mat');
unit_modules_400 = zeros(1,400);
for i = 1:400
    unit_modules_400(i)=norm(v_all(i,1))^2;
end

load('new_IRS_N300_p1.mat');
unit_modules_300 = zeros(1,300);
for i = 1:300
    unit_modules_300(i)=norm(v_all(i,1))^2;
end

figure      
subplot(2,1,1);
plot( 1:400,unit_modules_400(1:400),'r-','LineWidth',1.2);
axis([1,400,0.9999,1.0001]);
legend('N=400','Location','northeast')
hold on
subplot(2,1,2);
plot( 1:300,unit_modules_300(1:300),'b-','LineWidth',1.2);
axis([1,300,0.9999,1.0001]);
legend('N=300','Location','northeast')
hold on
xlabel({'Reflecting element index'},'FontSize',19);
t1=ylabel('The value of phase shift modulus','FontSize',19);

% axis([1,500,0.95,1.01]);
% grid on

% set(gca, 'XTick', [1 100 200 3300 400 500], 'YTick', [0.95 1 1.01],'FontSize',13); 
%legend('Proposed Algorithm 1', 'Location','northeast','FontSize',13)