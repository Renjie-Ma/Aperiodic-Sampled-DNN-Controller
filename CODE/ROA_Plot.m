%%%%%%%%% ROA_Plot
%%%%%%%%% Renjie Ma, Harbin Institute of Technology
%%%%%%%%% Dec 2023

clear
clc

n =2;

folder = 'PlotROA/'; 
load([folder 'P1.csv']) % 0.25 0.97
load([folder 'P2.csv']) % 0.35 0.97
load([folder 'P3.csv']) % 0.45 0.97
load([folder 'P4.csv']) % 0.25 1
load([folder 'P5.csv']) % 0.35 1
load([folder 'P6.csv']) % 0.45 1
load([folder 'P7.csv']) % 0.25 1.03
load([folder 'P8.csv']) % 0.35 1.03
load([folder 'P9.csv']) % 0.45 1.03

Ellipsoid1=ellipsoid(inv(P1(1:n,1:n)),zeros(n,1));
Ellipsoid2=ellipsoid(inv(P2(1:n,1:n)),zeros(n,1));
Ellipsoid3=ellipsoid(inv(P3(1:n,1:n)),zeros(n,1));
Ellipsoid4=ellipsoid(inv(P4(1:n,1:n)),zeros(n,1));
Ellipsoid5=ellipsoid(inv(P5(1:n,1:n)),zeros(n,1));
Ellipsoid6=ellipsoid(inv(P6(1:n,1:n)),zeros(n,1));
Ellipsoid7=ellipsoid(inv(P7(1:n,1:n)),zeros(n,1));
Ellipsoid8=ellipsoid(inv(P8(1:n,1:n)),zeros(n,1));
Ellipsoid9=ellipsoid(inv(P9(1:n,1:n)),zeros(n,1));

bar_varphi = 0.73; 
bbar_varphi = -bar_varphi;

figure (15)
xline(bar_varphi,'r')
hold on 
xline(bbar_varphi,'b')
hold on
plot(Ellipsoid1,[1,2],'-');
hold on 
plot(Ellipsoid2,[1,2],'-.','linewidth',1.0);
hold on 
plot(Ellipsoid3,[1,2],'--');
hold on 
plot(Ellipsoid4,[1,2],'-');
hold on 
plot(Ellipsoid5,[1,2],'-.','linewidth',1.0);
hold on 
plot(Ellipsoid6,[1,2],'--');
hold on 
plot(Ellipsoid7,[1,2],'-');
hold on 
plot(Ellipsoid8,[1,2],'-.','linewidth',1.0);
hold on 
plot(Ellipsoid9,[1,2],'--');
hold on 
plot(0,0,'s');
hold off
hh = legend('$\underline{x}_{1}(\mathrm{m}) $','$\overline{x}_{1}(\mathrm{m})  $',...
    '$\delta_{\rho}=0.25, \delta_{\beta}=0.97$','$\delta_{\rho}=0.35, \delta_{\beta}=0.97$',...
    '$\delta_{\rho}=0.45, \delta_{\beta}=0.97$',...
    '$\delta_{\rho}=0.25, \delta_{\beta}=1.00$',...
    '$\delta_{\rho}=0.35, \delta_{\beta}=1.00$','$\delta_{\rho}=0.45, \delta_{\beta}=1.00$',...
    '$\delta_{\rho}=0.25, \delta_{\beta}=1.03$','$\delta_{\rho}=0.35, \delta_{\beta}=1.03$',...
    '$\delta_{\rho}=0.45, \delta_{\beta}=1.03$', 'Equilibrium','NumColumns', 2);
set(hh,'Interpreter','latex');
title('\fontsize{11} The ellipsoidal inner approximations of robust RoA (STC)') 
xl=xlabel('Angular position $x_{1}(\mathrm{m})$');
set(xl,'Interpreter','latex');
yl=ylabel('Angular velocity $x_{2}(\mathrm{m/s})$');
set(yl,'Interpreter','latex');
ylim([-5,8]);

%% 

folder = 'PlotROA/'; 
load([folder 'PP1.csv']) % vartheta = 1;
load([folder 'PP2.csv']) % vartheta = 2
load([folder 'PP3.csv']) % vartheta = 3;
load([folder 'PP4.csv']) % vartheta = 4;
load([folder 'PP5.csv']) % vartheta = 5;


Ellipsoid10=ellipsoid(inv(PP1(1:n,1:n)),zeros(n,1));
Ellipsoid11=ellipsoid(inv(PP2(1:n,1:n)),zeros(n,1));
Ellipsoid12=ellipsoid(inv(PP3(1:n,1:n)),zeros(n,1));
Ellipsoid13=ellipsoid(inv(PP4(1:n,1:n)),zeros(n,1));
Ellipsoid14=ellipsoid(inv(PP5(1:n,1:n)),zeros(n,1));
%% 

figure (16)
xline(bar_varphi,'r')
hold on 
xline(bbar_varphi,'b')
hold on
plot(Ellipsoid14,[1,2],'-','Color', [0.4660 0.6740 0.1880]);
hold on 
plot(Ellipsoid13,[1,2],'-.','Color', [0 0.4470 0.7410]);
hold on 
plot(Ellipsoid12,[1,2],'-','Color', [0.8500 0.3250 0.0980]);
hold on 
plot(Ellipsoid10,[1,2],'-.','Color', [1 0 1]);
hold on 
plot(Ellipsoid11,[1,2],'-','Color', [0.4940 0.1840 0.5560]);
hold on 
plot(0,0,'s');
hold off
hh = legend('$\underline{x}_{1}(\mathrm{rad}) $','$\overline{x}_{1}(\mathrm{rad}) $',... 
    '$\vartheta = 5$','$\vartheta =4$',...
    '$\vartheta = 3$',...
    '$\vartheta = 2$',...
    '$\vartheta = 1$', 'Equilibrium','NumColumns', 3);
set(hh,'Interpreter','latex');
title('\fontsize{11} The ellipsoidal inner approximations of robust RoA (ETC)') 
xl=xlabel('Angular position $x_{1}(\mathrm{rad})$');
set(xl,'Interpreter','latex');
yl=ylabel('Angular velocity $x_{2}(\mathrm{rad/s})$');
set(yl,'Interpreter','latex');
ylim([-4 5])


%% 

figure (166)
xline(bar_varphi,'r')
hold on 
xline(bbar_varphi,'b')
hold on
plot(Ellipsoid14,[1,2],'-','Color', [0.4660 0.6740 0.1880]);
hold on 
plot(Ellipsoid13,[1,2],'-.','Color', [0 0.4470 0.7410]);
hold on 
plot(Ellipsoid12,[1,2],'-','Color', [0.8500 0.3250 0.0980]);
hold on 
plot(Ellipsoid10,[1,2],'-.','Color', [1 0 1]);
hold on 
plot(Ellipsoid11,[1,2],'-','Color', [0.4940 0.1840 0.5560]);
hold on 
plot(0,0,'s');
hold off
hh = legend('$\underline{x}_{1} $','$\overline{x}_{1} $',... 
    '$\vartheta = 5$','$\vartheta =4$',...
   '$\vartheta = 3$',...
    '$\vartheta = 2$',...
    '$\vartheta = 1$', 'Equilibrium');
set(hh,'Interpreter','latex');
title('\fontsize{11} The ellipsoidal inner approximations of robust RoA (ETC)') 
xl=xlabel('Angular position $x_{1}$');
set(xl,'Interpreter','latex');
yl=ylabel('Angular velocity $x_{2}$');
set(yl,'Interpreter','latex');
axis([-0.1, 0.1, 3.32, 3.38])








