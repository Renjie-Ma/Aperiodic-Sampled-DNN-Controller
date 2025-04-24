clear
clc

% Initial Parameters
Nstep = 800;                % Number of time steps
Ts = 0.01;                  % Sampling time
for k = 1: Nstep % 时域上界为Nstep
    time(k) = Ts*(k-1); 
end
up_bound = 0.73*ones(1,Nstep);
x_0c1 = [0.19, 3.5];          % Initial values case 1
folder = sprintf('%.2f_%.2f/', x_0c1(1), x_0c1(2)); 
folder = strrep(folder, '-', 'm');  % Replace '-' with 'minus'

% Load data from CSV files
etaa11 = load([folder 'staa1.csv']);
etaa12 = load([folder 'staa2.csv']);
etaa13 = load([folder 'staa3.csv']);
inti1E = load([folder 'STintiE.csv']) ;
seqi1 = load([folder 'STseqi.csv']) ;

exp1 = load([folder 'expert_state.csv']);

x_0c2 = [0.43, 3.0]; % Initial values case 2

folder = sprintf('%.2f_%.2f/', x_0c2(1), x_0c2(2)); 
folder = strrep(folder, '-', 'm');  % Replace '-' with 'minus

% Load data from CSV files
etaa21 = load([folder 'staa1.csv']);
etaa22 = load([folder 'staa2.csv']);
etaa23 = load([folder 'staa3.csv']);
inti2E = load([folder 'STintiE.csv']) ;
seqi2 = load([folder 'STseqi.csv']) ;
exp2 = load([folder 'expert_state.csv']);

x_0c3 = [-0.33, -3.3]; % Initial values case 3

folder = sprintf('%.2f_%.2f/', x_0c3(1), x_0c3(2)); 
folder = strrep(folder, '-', 'm');  % Replace '-' with 'minus

% Load data from CSV files
etaa31 = load([folder 'staa1.csv']);
etaa32 = load([folder 'staa2.csv']);
etaa33 = load([folder 'staa3.csv']);
inti3E = load([folder 'STintiE.csv']) ;
seqi3 = load([folder 'STseqi.csv']) ;
exp3 = load([folder 'expert_state.csv']);
%% 

figure(11)
subplot(3,1,1)


plot(time,etaa11,'-','linewidth',1);
hold on 
plot(time,etaa21,'-o','MarkerIndices',1:50:length(etaa21),'linewidth',1);
hold on
plot(time,etaa31,'-.','linewidth',1);
hold on
plot (time,up_bound,'--', 'linewidth', 1);
hold off
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex');
current_ticks = yticks;

% 添加 0.73 到默认刻度中
if ~ismember(0.73, current_ticks) % 检查 0.73 是否已存在
    new_ticks = sort([current_ticks, 0.73]); % 排序
    yticks(new_ticks); % 更新刻度
    yticklabels(arrayfun(@num2str, new_ticks, 'UniformOutput', false)); % 更新标签
end
hh = legend('$x_{1}$ (Case 1)','$x_{1}$ (Case 2)','$x_{1}$ (Case 3)','FontSize',10,'Orientation','horizon');

if false
    plot(time, exp1(1, 1:Nstep), '-', 'linewidth', 1.5)
    hold on
    plot(time, exp2(1, 1:Nstep), '-', 'linewidth', 1.5)
    hold on
    plot(time, exp3(1, 1:Nstep), '-', 'linewidth', 1.5)
    hold off
    hh = legend('$x_{1}Expert$ (Case 1)','$x_{1}Expert$ (Case 2)','$x_{1}Expert$ (Case 3)','$x_{1}$ (Case 1)','$x_{1}$ (Case 2)','$x_{1}$ (Case 3)','FontSize',10,'Orientation','horizon');
end

set(hh, 'Interpreter', 'latex');
t='The Self-Triggered State Response of Angular Position $x_{1}(\mathrm{rad})$';
%title('\fontsize{11}(a) ');
title(t, 'Interpreter', 'latex');
%title('\fontsize{11}(a) The self-triggered state response of angular position') 

subplot(3,1,2)
plot(time,etaa12,'-','linewidth',1);
hold on 
plot(time,etaa22,'-o','MarkerIndices',1:50:length(etaa22),'linewidth',1);
hold on
plot(time,etaa32,'-.','linewidth',1);
hold off
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex');
hh = legend('$x_{2}$ (Case 1)','$x_{2}$ (Case 2)','$x_{2}$ (Case 3)','FontSize',10,'Orientation','horizon');
if false
    plot(time, exp1(2, 1:Nstep), '-', 'linewidth', 1.5)
    hold on
    plot(time, exp2(2, 1:Nstep), '-', 'linewidth', 1.5)
    hold on
    plot(time, exp3(2, 1:Nstep), '-', 'linewidth', 1.5)
    hold off
    hh = legend('$x_{1}Expert$ (Case 1)','$x_{1}Expert$ (Case 2)','$x_{1}Expert$ (Case 3)','$x_{1}$ (Case 1)','$x_{1}$ (Case 2)','$x_{1}$ (Case 3)','FontSize',10,'Orientation','horizon');
end
set(hh, 'Interpreter', 'latex');
t='The Self-Triggered State Response of Angular Velocity $x_{2}(\mathrm{rad/s})$';
%title('\fontsize{11}(a) ');
title(t, 'Interpreter', 'latex');
%title('\fontsize{11}(b) The self-triggered state response of angular velocity')  

subplot(3,1,3)
plot(time,etaa13,'-','linewidth',1);
hold on 
plot(time,etaa23,'-o','MarkerIndices',1:50:length(etaa23),'linewidth',1);
hold on
plot(time,etaa33,'-.','linewidth',1);
hold off
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex');
hh = legend('$\xi$ (Case 1)','$\xi$ (Case 2)',' $\xi$ (Case 3)','FontSize',10,'Orientation','horizon');
set(hh, 'Interpreter', 'latex');
t='The Self-Triggered State Response of IQC Virtual-Filter State $\xi(\mathrm{rad})$';
%title('\fontsize{11}(a) ');
title(t, 'Interpreter', 'latex');
%title('\fontsize{11}(c) The self-triggered state response of IQC virtual filter ') 
%% 


figure(12)
subplot(3,1,1)
stem(seqi1,inti1E,'color','(0.0,0.45,0.74)','linewidth',1.0);
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex'); 
% yl=ylabel('Inter-event intervals');
% yl=ylabel('$d_{j}$');
% set(yl,'Interpreter','latex','fontsize',13); 

title('(a) The self-triggered instants and durations (Case 1)', 'Interpreter', 'latex')
axis([0,200,1,12])
current_xticks = get(gca, 'XTick');
new_xticks = current_xticks * 0.02;
set(gca, 'XTickLabel', new_xticks);
current_yticks = get(gca, 'YTick');
new_yticks = current_yticks * 0.02;
set(gca, 'YTickLabel', new_yticks);  
subplot(3,1,2)


stem(seqi2,inti2E,'color','(0.0,0.45,0.74)','linewidth',1.0);
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex'); 
% yl=ylabel('Inter-event intervals');
% yl=ylabel('$d_{j}$');
% set(yl,'Interpreter','latex','fontsize',13); 
title('(b) The self-triggered instants and durations (Case 2)', 'Interpreter', 'latex') 
axis([0,200,1,12])
current_xticks = get(gca, 'XTick');
new_xticks = current_xticks * 0.02;
set(gca, 'XTickLabel', new_xticks);
current_yticks = get(gca, 'YTick');
new_yticks = current_yticks * 0.02;
set(gca, 'YTickLabel', new_yticks);  


subplot(3,1,3)
stem(seqi3,inti3E,'color','(0.0,0.45,0.74)','linewidth',1.0);
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex'); 
% yl=ylabel('Inter-event intervals');
% yl=ylabel('$d_{j}$');
% set(yl,'Interpreter','latex','fontsize',13); 
title('(c) The self-triggered instants and durations (Case 3)', 'Interpreter', 'latex')
axis([0,200,1,12])
current_xticks = get(gca, 'XTick');
new_xticks = current_xticks * 0.02;
set(gca, 'XTickLabel', new_xticks);
current_yticks = get(gca, 'YTick');
new_yticks = current_yticks * 0.02;
set(gca, 'YTickLabel', new_yticks);  