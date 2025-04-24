%%%%%%%%% ROA_Plot
%%%%%%%%% Renjie Ma, Harbin Institute of Technology
%%%%%%%%% Dec 2023

clear
clc

n =2;

folder = 'PlotROA/'; 
load([folder 'PS1.csv']) % 0.25 0.97
load([folder 'PS2.csv']) % 0.35 0.97
load([folder 'PS3.csv']) % 0.45 0.97
load([folder 'PS4.csv']) % 0.25 1
load([folder 'PS5.csv']) % 0.35 1
load([folder 'PS6.csv']) % 0.45 1
load([folder 'PS7.csv']) % 0.25 1.03
load([folder 'PS8.csv']) % 0.35 1.03
load([folder 'PS9.csv']) % 0.45 1.03

Ellipsoid1=ellipsoid(inv(PS1(1:n,1:n)),zeros(n,1));
Ellipsoid2=ellipsoid(inv(PS2(1:n,1:n)),zeros(n,1));
Ellipsoid3=ellipsoid(inv(PS3(1:n,1:n)),zeros(n,1));
Ellipsoid4=ellipsoid(inv(PS4(1:n,1:n)),zeros(n,1));
Ellipsoid5=ellipsoid(inv(PS5(1:n,1:n)),zeros(n,1));
Ellipsoid6=ellipsoid(inv(PS6(1:n,1:n)),zeros(n,1));
Ellipsoid7=ellipsoid(inv(PS7(1:n,1:n)),zeros(n,1));
Ellipsoid8=ellipsoid(inv(PS8(1:n,1:n)),zeros(n,1));
Ellipsoid9=ellipsoid(inv(PS9(1:n,1:n)),zeros(n,1));

bar_varphi = 0.73; 
bbar_varphi = -bar_varphi;

figure (25)
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
hh = legend('$\underline{x}_{1}(\mathrm{m}) $','$\overline{x}_{1}(\mathrm{m}) $',...
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
ylim([-5 8])

figure (26)
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
hh = legend('$\underline{x}_{1} $','$\overline{x}_{1} $',...
    '$\delta_{\rho}=0.25, \delta_{\beta}=0.97$','$\delta_{\rho}=0.35, \delta_{\beta}=0.97$',...
    '$\delta_{\rho}=0.45, \delta_{\beta}=0.97$',...
    '$\delta_{\rho}=0.25, \delta_{\beta}=1.00$',...
    '$\delta_{\rho}=0.35, \delta_{\beta}=1.00$','$\delta_{\rho}=0.45, \delta_{\beta}=1.00$',...
    '$\delta_{\rho}=0.25, \delta_{\beta}=1.03$','$\delta_{\rho}=0.35, \delta_{\beta}=1.03$',...
    '$\delta_{\rho}=0.45, \delta_{\beta}=1.03$', 'Equilibrium');
set(hh,'Interpreter','latex');
title('\fontsize{11} The ellipsoidal inner approximations of robost RoA (STC)') 
axis([-0.1, 0.1, 1, 5])

xl=xlabel('Angular position $x_{1}$');
set(xl,'Interpreter','latex');
yl=ylabel('Angular velocity $x_{2}$');
set(yl,'Interpreter','latex');




