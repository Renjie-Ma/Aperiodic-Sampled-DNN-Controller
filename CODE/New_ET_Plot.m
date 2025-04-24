%%%%%%%%% ET_Plot
%%%%%%%%% Renjie Ma, Harbin Institute of Technology
%%%%%%%%% Dec 2023

clear
clc

% Initial Parameters
Nstep = 800;                % Number of time steps
Ts = 0.01;                  % Sampling time
x_0 = [0.19, 3.5];          % Initial values
folder = sprintf('%.2f_%.2f/', x_0(1), x_0(2)); 
folder = strrep(folder, '-', 'minus');  % Replace '-' with 'minus'

% Time vector
time = Ts * (0:Nstep-1); 
up_bound = 0.73*ones(1,Nstep);
% Load data from CSV files
load([folder 'etaa11.csv']) 
load([folder 'etaa12.csv']) 
load([folder 'etaa13.csv']) 
load([folder 'intiE1.csv']) 
load([folder 'seqi1.csv']) 

load([folder 'etaa21.csv']) 
load([folder 'etaa22.csv']) 
load([folder 'etaa23.csv']) 
load([folder 'intiE2.csv']) 
load([folder 'seqi2.csv']) 

load([folder 'etaa31.csv']) 
load([folder 'etaa32.csv']) 
load([folder 'etaa33.csv']) 
load([folder 'intiE3.csv']) 
load([folder 'seqi3.csv']) 

load([folder 'etaa41.csv']) 
load([folder 'etaa42.csv']) 
load([folder 'etaa43.csv']) 
load([folder 'intiE4.csv']) 
load([folder 'seqi4.csv']) 

load([folder 'expert_state.csv']) 
load([folder 'expert_xi.csv']) 
%% 

% Plotting Figure 4 (The event-triggered state responses)
figure(9)
set(gcf, 'Position', [100, 100, 800, 600])  % Resize figure

% Subplot 1: Angular Position x1
subplot(3, 1, 1)

plot(time, expert_state(1, 1:Nstep), '-', 'linewidth', 1)
hold on
plot(time, etaa11, '-', 'linewidth', 1);
hold on
plot(time, etaa21, '-o','MarkerIndices',1:50:length(etaa21), 'linewidth', 1);
hold on
plot(time, etaa31, '-', 'linewidth', 1);
hold on
plot(time, etaa41, '-.', 'linewidth', 1);
hold on
plot (time,up_bound,'--', 'linewidth', 1);
hold off
xlabel('Time (Sec.)', 'Interpreter', 'latex')
current_ticks = yticks;
ylim([-0.3,0.8])
% 添加 0.73 到默认刻度中
if ~ismember(0.73, current_ticks) % 检查 0.73 是否已存在
    new_ticks = sort([current_ticks, 0.73]); % 排序
    yticks(new_ticks); % 更新刻度
    yticklabels(arrayfun(@num2str, new_ticks, 'UniformOutput', false)); % 更新标签
end
legend('$x_{1}$ (MPC)', '$x_{1}$ (Case 1)', '$x_{1}$ (Case 2)', '$x_{1}$ (Case 3)', '$x_{1}$ (Case 4)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
title('The Event-Triggered and MPC State Response of Angular Position $x_{1}(\mathrm{rad})$', 'Interpreter', 'latex')

% Subplot 2: Angular Velocity x2
subplot(3, 1, 2)

plot(time, expert_state(2, 1:Nstep), '-', 'linewidth', 1)
hold on
plot(time, etaa12, '-', 'linewidth', 1);
hold on
plot(time, etaa22, '-o','MarkerIndices',1:50:length(etaa21), 'linewidth', 1);
hold on
plot(time, etaa32, '-', 'linewidth', 1);
hold on
plot(time, etaa42, '-.', 'linewidth', 1);
hold off
xlabel('Time (Sec.)', 'Interpreter', 'latex')
ylim([-6,6])
legend('$x_{2}$ (MPC)', '$x_{2}$ (Case 1)', '$x_{2}$ (Case 2)', '$x_{2}$ (Case 3)', '$x_{2}$ (Case 4)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
title('The Event-Triggered and MPC State Response of Angular Velocity $x_{2}(\mathrm{rad/s})$', 'Interpreter', 'latex')

% Subplot 3: IQC Virtual-Filter State xi
subplot(3, 1, 3)


plot(time, expert_xi(1:Nstep),'-', 'linewidth', 1);
hold on
plot(time, etaa13, '-', 'linewidth', 1);
hold on
plot(time, etaa23, '-o', 'MarkerIndices',1:50:length(etaa23),'linewidth', 1);
hold on
plot(time, etaa33, '-', 'linewidth', 1);
hold on
plot(time, etaa43, '-.', 'linewidth', 1);
hold off
xlabel('Time (Sec.)', 'Interpreter', 'latex')
legend('$\xi$ (MPC)','$\xi$ (Case 1)', '$\xi$ (Case 2)', '$\xi$ (Case 3)', '$\xi$ (Case 4)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
title('The Event-Triggered and MPC State Response of IQC Virtual-Filter State $\xi(\mathrm{rad})$', 'Interpreter', 'latex')
%% 

if true
    % Plotting Figure 5 (Inter-event intervals and timings)
    figure(5)
    set(gcf, 'Position', [100, 100, 1000, 800])  % Resize figure
    
    % Subplot 1: Inter-event intervals for Case 1
    subplot(2, 2, 1)
    stem(seqi1, intiE1, 'color', [0.0, 0.45, 0.74], 'linewidth', 1)
    xlabel('Time(Sec.)', 'Interpreter', 'latex')
    legend('$\vartheta = 4$ (Case 1)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
    axis([0, 150, 0, 50])
    title('Inter-event Instants (Case 1)', 'Interpreter', 'latex')
    current_xticks = get(gca, 'XTick');
    new_xticks = current_xticks * 0.02;
    set(gca, 'XTickLabel', new_xticks);
    current_yticks = get(gca, 'YTick');
    new_yticks = current_yticks * 0.02;
    set(gca, 'YTickLabel', new_yticks);
    % Subplot 2: Inter-event intervals for Case 2
    subplot(2, 2, 2)
    stem(seqi2, intiE2, 'color', [0.0, 0.45, 0.74], 'linewidth', 1)
    xlabel('Time(Sec.)', 'Interpreter', 'latex')
    legend('$\vartheta = 3$ (Case 2)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
    axis([0, 150, 0, 40])
    title('Inter-event Instants (Case 2)', 'Interpreter', 'latex')
    current_xticks = get(gca, 'XTick');
    new_xticks = current_xticks * 0.02;
    set(gca, 'XTickLabel', new_xticks);
    current_yticks = get(gca, 'YTick');
    new_yticks = current_yticks * 0.02;
    set(gca, 'YTickLabel', new_yticks);    
    % Subplot 3: Inter-event intervals for Case 3
    subplot(2, 2, 3)
    stem(seqi3, intiE3, 'color', [0.0, 0.45, 0.74], 'linewidth', 1)
    xlabel('Time(Sec.)', 'Interpreter', 'latex')
    legend('$\vartheta = 2$ (Case 3)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
    axis([0, 150, 0, 30])
    title('Inter-event Instants (Case 3)', 'Interpreter', 'latex')
    current_xticks = get(gca, 'XTick');
    new_xticks = current_xticks * 0.02;
    set(gca, 'XTickLabel', new_xticks);
    current_yticks = get(gca, 'YTick');
    new_yticks = current_yticks * 0.02;
    set(gca, 'YTickLabel', new_yticks);    
    % Subplot 4: Inter-event intervals for Case 4
    subplot(2, 2, 4)
    stem(seqi4, intiE4, 'color', [0.0, 0.45, 0.74], 'linewidth', 1)
    xlabel('Time(Sec.)', 'Interpreter', 'latex')
    legend('$\vartheta = 1$ (Case 4)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
    axis([0, 150, 0, 20])
    title('Inter-event Instants (Case 4)', 'Interpreter', 'latex')
    current_xticks = get(gca, 'XTick');
    new_xticks = current_xticks * 0.02;
    set(gca, 'XTickLabel', new_xticks);
    current_yticks = get(gca, 'YTick');
    new_yticks = current_yticks * 0.02;
    set(gca, 'YTickLabel', new_yticks);
end

if false

    figure(51)
    set(gcf, 'Position', [100, 100, 1000, 800])  % Resize figure
    
    % Subplot 1: Inter-event intervals for Case 1
    subplot(4, 1, 1)
    stem(seqi1, intiE1, 'color', [0.0, 0.45, 0.74], 'linewidth', 1)
    xlabel('Sampling Time', 'Interpreter', 'latex')
    legend('$\vartheta = 4$ (Case 1)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
    axis([0, 400, 0, 50])
    title('Inter-event Instants (Case 1)', 'Interpreter', 'latex')
    current_xticks = get(gca, 'XTick');
    new_xticks = current_xticks * 0.02;
    set(gca, 'XTickLabel', new_xticks);
    current_yticks = get(gca, 'YTick');
    new_yticks = current_yticks * 0.02;
    set(gca, 'YTickLabel', new_yticks);
    % Subplot 2: Inter-event intervals for Case 2
    subplot(4, 1, 2)
    stem(seqi2, intiE2, 'color', [0.0, 0.45, 0.74], 'linewidth', 1.5)
    xlabel('Sampling Time', 'Interpreter', 'latex')
    legend('$\vartheta = 3$ (Case 2)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
    axis([0, 400, 0, 50])
    title('Inter-event Instants (Case 2)', 'Interpreter', 'latex')
    
    % Subplot 3: Inter-event intervals for Case 3
    subplot(4, 1, 3)
    stem(seqi3, intiE3, 'color', [0.0, 0.45, 0.74], 'linewidth', 1.5)
    xlabel('Sampling Time', 'Interpreter', 'latex')
    legend('$\vartheta = 2$ (Case 3)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
    axis([0, 400, 0, 50])
    title('Inter-event Instants (Case 3)', 'Interpreter', 'latex')
    
    % Subplot 4: Inter-event intervals for Case 4
    subplot(4, 1, 4)
    stem(seqi4, intiE4, 'color', [0.0, 0.45, 0.74], 'linewidth', 1.5)
    xlabel('Sampling Time', 'Interpreter', 'latex')
    legend('$\vartheta = 1$ (Case 4)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
    axis([0, 400, 0, 50])
    title('Inter-event Instants (Case 4)', 'Interpreter', 'latex')

end


if true

    figure(52)
    set(gcf, 'Position', [100, 100, 1000, 800])  % Resize figure
    
    % Subplot 1: Inter-event intervals for Case 1
    subplot(4, 1, 1)
    stem(seqi1, intiE1, 'color', [0.0, 0.45, 0.74], 'linewidth', 1.5)
    xlabel('Time (Sec.)', 'Interpreter', 'latex')
    legend('$\vartheta = 4$ (Case 1)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
    axis([0, 400, 0, 50])
   % title('Inter-event Instants (Case 1)', 'Interpreter', 'latex')
    current_xticks = get(gca, 'XTick');
    new_xticks = current_xticks * 0.02;
    set(gca, 'XTickLabel', new_xticks);


    % Subplot 2: Inter-event intervals for Case 2
    subplot(4, 1, 2)
    stem(seqi2, intiE2, 'color', [0.0, 0.45, 0.74], 'linewidth', 1.5)
    xlabel('Time (Sec.)', 'Interpreter', 'latex')
    legend('$\vartheta = 3$ (Case 2)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
    axis([0, 400, 0, 40])
   % title('Inter-event Instants (Case 2)', 'Interpreter', 'latex')
    current_xticks = get(gca, 'XTick');
    new_xticks = current_xticks * 0.02;
    set(gca, 'XTickLabel', new_xticks);

    % Subplot 3: Inter-event intervals for Case 3
    subplot(4, 1, 3)
    stem(seqi3, intiE3, 'color', [0.0, 0.45, 0.74], 'linewidth', 1.5)
    xlabel('Time(Sec.)', 'Interpreter', 'latex')
    legend('$\vartheta = 2$ (Case 3)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
    axis([0, 400, 0, 40])
   % title('Inter-event Instants (Case 3)', 'Interpreter', 'latex')
    current_xticks = get(gca, 'XTick');
    new_xticks = current_xticks * 0.02;
    set(gca, 'XTickLabel', new_xticks);
    % Subplot 4: Inter-event intervals for Case 4
    subplot(4, 1, 4)
    stem(seqi4, intiE4, 'color', [0.0, 0.45, 0.74], 'linewidth', 1.5)
    xlabel('Time (Sec.)', 'Interpreter', 'latex')
    legend('$\vartheta = 1$ (Case 4)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
    axis([0, 400, 0, 20])
    title('Inter-event Instants (Case 4)', 'Interpreter', 'latex')
    current_xticks = get(gca, 'XTick');
    new_xticks = current_xticks * 0.02;
    set(gca, 'XTickLabel', new_xticks);
end