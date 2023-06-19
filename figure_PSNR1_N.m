noIRS = 8.0139*ones(1,7);
IRS =     [10.0106		14.8537	19.1104	26.2808	34.7315		47.102	65.3631];
randIRS = [8.2489		8.4035	8.4719	8.6108	8.7867		8.9012	9.1218];
IRSb2 =   [9.6788		13.2704	16.4307	22.3475	27.839		36.767	51.8386];
IRSb3 =   [9.8339		14.4672	18.2693	25.0433	32.841		44.001	61.1921];
N_number = [10			50	100	150	200	300 400];

alpha = 0;
c1 = 0.905;
c2 = 1.34;
B = 200;

theta1 = 13870;
beta1 = 493.2;
theta2 = 2876;
beta2 = 23.6;

Q1_noIRS = -10*log10( theta1./( c1*B*log2(1+noIRS./c2)-beta1 ) -alpha ) + 20*log10(255);
Q1_IRS = -10*log10( theta1./( c1*B*log2(1+IRS./c2)-beta1 ) -alpha ) + 20*log10(255);
Q1_randIRS = -10*log10( theta1./( c1*B*log2(1+randIRS./c2)-beta1 ) -alpha ) + 20*log10(255);
Q1_IRSb2 = -10*log10( theta1./( c1*B*log2(1+IRSb2./c2)-beta1 ) -alpha ) + 20*log10(255);
Q1_IRSb3 = -10*log10( theta1./( c1*B*log2(1+IRSb3./c2)-beta1 ) -alpha ) + 20*log10(255);

Q2_noIRS = -10*log10( theta2./( c1*B*log2(1+noIRS./c2)-beta2 ) -alpha ) + 20*log10(255);
Q2_IRS = -10*log10( theta2./( c1*B*log2(1+IRS./c2)-beta2 ) -alpha ) + 20*log10(255);
Q2_randIRS = -10*log10( theta2./( c1*B*log2(1+randIRS./c2)-beta2 ) -alpha ) + 20*log10(255);
Q2_IRSb2 = -10*log10( theta2./( c1*B*log2(1+IRSb2./c2)-beta2 ) -alpha ) + 20*log10(255);
Q2_IRSb3 = -10*log10( theta2./( c1*B*log2(1+IRSb3./c2)-beta2 ) -alpha ) + 20*log10(255);


%Q1=-10*log10( theta1/( c1*B*log2(1+average_muu/c2)-beta1 ) -alpha ) + 20*log10(255)

figure    
plot( N_number,Q1_IRS,'r-o','LineWidth',1.3);
hold on
plot( N_number,Q1_randIRS,'b->','LineWidth',1.3);
hold on
plot( N_number,Q1_noIRS,'k-d','LineWidth',1.3);
hold on
plot( N_number,Q1_IRSb2,'-*','LineWidth',1.3);
hold on
plot( N_number,Q1_IRSb3,'-s','LineWidth',1.3);
hold on
axis([10,400,18,35]);
grid on
xlabel({'The number of reflecting elements $N$'},'Interpreter','latex','FontSize',19);
ylabel('Max-min average PSNR (dB)','FontSize',19);
set(gca, 'XTick', N_number,'YTick', [18 20 25 30 35],'FontSize',15); 
h=legend('Algorithm 1','Random phase shifts','Without IRS','Discrete phase shifts with $b=2$','Discrete phase shifts with $b=3$','Interpreter','latex', 'Location','east');
set(h,'FontSize',13);


figure    
plot( N_number,Q2_IRS,'r-o','LineWidth',1.3);
hold on
plot( N_number,Q2_randIRS,'b->','LineWidth',1.3);
hold on
plot( N_number,Q2_noIRS,'k-d','LineWidth',1.3);
hold on
plot( N_number,Q2_IRSb2,'-*','LineWidth',1.3);
hold on
plot( N_number,Q2_IRSb3,'-s','LineWidth',1.3);
hold on
axis([10,400,40,44]);
grid on
xlabel({'The number of reflecting elements $N$'},'Interpreter','latex','FontSize',19);
ylabel('Max-min average PSNR (dB)','FontSize',19);
set(gca, 'XTick', N_number,'YTick', [40 41 42 43 44],'FontSize',15); 
h=legend('Algorithm 1','Random phase shifts','Without IRS','Discrete phase shifts with $b=2$','Discrete phase shifts with $b=3$','Interpreter','latex', 'Location','east');
set(h,'FontSize',13);