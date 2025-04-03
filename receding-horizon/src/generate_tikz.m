clc; 

% To generate tikz files. 
addpath('matlab2tikz/src')
addpath('matlab2tikz/src/dev/')
addpath('matlab2tikz/src/private/')

netload_MC = csvread('simulation_data/MC_netload.csv');
% e_MC = csvread('simulation_data/MC_e_1.csv');
% e_ref_MC = csvread('simulation_data/MC_e_ref_1.csv');
% zeta_MC = csvread('simulation_data/MC_q_4.csv');
lower_envelope = zeros(50, 1);
upper_envelope = zeros(50, 1);
avg_envelope = zeros(50, 1);

for i = 1:50 
    lower_envelope(i) = min(netload_MC(i,:));
    upper_envelope(i) = max(netload_MC(i, :));
    % avg_envelope(i) = sum(zeta_MC(i, :)) / MONTE_CARLO_SIMULATIONS; 
end

x_ticks = [0 10 20 30 40 50]
% y_ticks = [0 10 20 30 40 50]

figure(1)
yline(185, 'LineStyle', '-.', 'DisplayName', 'Max', 'color', 'black', 'LineWidth', 1.1);
text(30, 190, 'L_{max} = 185 kWh');
% yline(-7.3182, 'LineStyle', '-.', 'DisplayName', 'Min', 'color', 'black', 'LineWidth', 1.1);
% yline(26, 'LineStyle', '-.', 'DisplayName', 'Min', 'color', 'red', 'LineWidth', 1.1);
hold on; 
plot(0:49, lower_envelope, 'DisplayName', 'LowerBound', 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.1);
plot(0:49, upper_envelope, 'DisplayName', 'UpperBound', 'LineStyle', '--', 'Color', 'g', 'LineWidth', 1.1);
% legend; 
ylim([50 200]);
xlim([0 50]);
NUM_FIGS = 20;
colors = lines(NUM_FIGS);
ppp = 0; 
for i = 1:ceil(200/NUM_FIGS):200 
    ppp = ppp + 1;
    plot(0:49, netload_MC(:,i), 'Color', colors(ppp, :), 'DisplayName', 'Actual', 'LineWidth', 1.1);
end
% title('Agent 4');
xlabel('Time (in hr)')
ylabel('Net Load (in kWh)');
xticks(x_ticks);
% yticks(y_ticks);
hold off; 

% matlab2tikz('figures/p1_zeta.tex', 'figurehandle', figure(1), 'height', '0.5\textwidth');