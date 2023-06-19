noIRS =   [8.0139	8.827771598	9.639041732	10.42975569	11.23157531	12.03643808];
IRS =     [65.3631	71.8014	    78.5913	    84.859	    91.3289	    98.2893];
randIRS = [9.1218	10.03496286	10.95300602	11.83810382	12.83595269	13.6237];
IRSb2 =   [51.8386	56.6697	    62.254	    67.105	    72.6246	    77.5701];
IRSb3 =   [61.1921	67.3501	    73.7453	    79.3452	    85.5629	    92.0691];
POWER = [1	1.1	1.2	1.3	1.4	1.5];


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


figure    
plot( POWER,Q1_IRS,'r-o','LineWidth',1.3);
hold on
plot( POWER,Q1_randIRS,'b->','LineWidth',1.3);
hold on
plot( POWER,Q1_noIRS,'k-d','LineWidth',1.3);
hold on
plot( POWER,Q1_IRSb2,'-*','LineWidth',1.3);
hold on
plot( POWER,Q1_IRSb3,'-s','LineWidth',1.3);
hold on
axis([1,1.5,15,35]);
grid on
xlabel({'The transmit power budget $P_{max}$ (watt)'},'Interpreter','latex','FontSize',19);
ylabel('Max-min average PSNR (dB)','FontSize',19);
set(gca, 'XTick', POWER,'FontSize',15); 
h=legend('Algorithm 1','Random phase shifts','Without IRS','Discrete phase shifts with $b=2$','Discrete phase shifts with $b=3$', 'Location','best','Interpreter','latex');
set(h,'FontSize',15);


figure    
plot( POWER,Q2_IRS,'r-o','LineWidth',1.3);
hold on
plot( POWER,Q2_randIRS,'b->','LineWidth',1.3);
hold on
plot( POWER,Q2_noIRS,'k-d','LineWidth',1.3);
hold on
plot( POWER,Q2_IRSb2,'-*','LineWidth',1.3);
hold on
plot( POWER,Q2_IRSb3,'-s','LineWidth',1.3);
hold on
%axis([1,1.5,15,35]);
grid on
xlabel({'The transmit power budget $P_{max}$ (watt)'},'Interpreter','latex','FontSize',19);
ylabel('Max-min average PSNR (dB)','FontSize',19);
set(gca, 'XTick', POWER,'FontSize',15); 
h=legend('Algorithm 1','Random phase shifts','Without IRS','Discrete phase shifts with $b=2$','Discrete phase shifts with $b=3$','Interpreter','latex', 'Location','best');
set(h,'FontSize',15);
