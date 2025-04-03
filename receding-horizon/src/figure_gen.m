PLOT = true; 
LAST_FIG = 1; 
if PLOT 
    for i = 1:N 
      figure(LAST_FIG)
      plot(0:T, zeta_over_time{i}, 'LineWidth', 5, 'color', colors(progress, :), 'DisplayName', 'State Trajectory');
      hold on; 
      yline(players{i}.ss_state(1), 'LineStyle', '--', 'LineWidth', 5, 'color', 'r', 'DisplayName', 'Steady State');

      xxx = players{i}.Z.extreme([1;0]);
      yline(xxx.x(1), 'LineWidth', 5, 'color', 'r', 'DisplayName', 'Max Disturbance');
      
      xxx = players{i}.Z.extreme([-1;0]);
      yline(xxx.x(1), 'LineStyle', '-.', 'LineWidth', 5, 'color', 'r', 'DisplayName', 'Min Disturbance');
      
      xlabel('Time (in hr)', 'FontWeight', 'bold');
      ylabel('$\zeta$ (in kWh)', 'FontWeight', 'bold');

      ylim([-20 20]);
      % hold off; 
      set(gca, "FontSize", 65, 'FontWeight', 'bold');
      set(groot, 'defaultTextFontWeight', 'bold');
      set(groot, 'defaultAxesFontName', 'Times New Roman'); % or 'Times'
      set(groot, 'defaultTextFontName', 'Times New Roman'); % or 'Times'
      set(groot, 'defaultAxesFontSize', 65, 'defaultAxesFontWeight', 'bold'); % Match LaTeX font size (10-12pt)
      set(groot, 'defaultTextInterpreter', 'latex'); % For LaTeX symbols like $\zeta$
      % legend("FontWeight", "normal", 'FontSize', 40); 
      LAST_FIG = LAST_FIG + 1;
    end
    for i = 1:N 
      figure(LAST_FIG) 
      plot(0:T, q_over_time{i}, 'LineWidth', 6, 'DisplayName', 'State Trajectory', 'color', colors(progress, :));
      yline(players{i}.ss_state(2) + q_ref{i}, 'LineStyle', '--', 'color', 'r', 'LineWidth', 6, 'DisplayName', 'Steady State');
      xlabel('Time (in hr)');
      ylabel('$q$ (in kWh)');
      title(sprintf('Agent %d',i));
      ylim([0 35]);
      set(gca, "FontSize", 65, 'FontWeight', 'bold');
      set(groot, 'defaultTextFontWeight', 'bold');
      set(groot, 'defaultAxesFontName', 'Times New Roman'); % or 'Times'
      set(groot, 'defaultTextFontName', 'Times New Roman'); % or 'Times'
      set(groot, 'defaultAxesFontSize', 65, 'defaultAxesFontWeight', 'bold'); % Match LaTeX font size (10-12pt)
      set(groot, 'defaultTextInterpreter', 'latex'); % For LaTeX symbols like $\zeta$
      % legend("FontWeight", "normal"); 
      LAST_FIG = LAST_FIG + 1;
      
      hold on; 
    end
%     figure(LAST_FIG)
%     l_over_time = e_over_time{1} + s_over_time{1};
%     for j = 2:N 
%        l_over_time = l_over_time + e_over_time{j} + s_over_time{j};
%     end
%     bar(0:T-1, l_over_time, 'LineWidth', 6, 'DisplayName', 'NetLoad');
%     hold on; 
%     yline(L_max, 'color', 'r', 'LineStyle', '--', 'LineWidth', 6, 'DisplayName', 'Maximum Capacity');
%     text(T-9, L_max-10.5, '$L_{max}$ = 190 kWh', 'HorizontalAlignment', 'center', 'FontSize', 40, 'FontWeight', 'bold');
%     % yline(L_min, 'color', 'r', 'LineStyle', '--', 'LineWidth', 5);
%     % text(T-10, L_min-1.5, 'Minimum Capacity', 'HorizontalAlignment', 'center');
%     xlabel('Time (in hr)');
%     ylabel('Net Load (kWh)');
%     % title('Net Load over time');
%     % hold off;
%     set(gca, 'Padding', 'loose', 'FontSize', 45);
%     set(groot, 'defaultTextFontWeight', 'bold');
%     set(groot, 'defaultAxesFontName', 'Times New Roman'); % or 'Times'
%     set(groot, 'defaultTextFontName', 'Times New Roman'); % or 'Times'
%     set(groot, 'defaultAxesFontSize', 45, 'defaultAxesFontWeight', 'bold'); % Match LaTeX font size (10-12pt)
%     set(groot, 'defaultTextInterpreter', 'latex'); % For LaTeX symbols like $\zeta$
%     % legend("FontWeight", "normal"); 
%     LAST_FIG = LAST_FIG + 1;
 end

if PLOT 
    for i = 1:N 
      figure(LAST_FIG) 
      plot(1:T, e_over_time{i}, 'LineWidth', 6, 'color', colors(progress, :), 'DisplayName', 'Actual');
      hold on; 
      plot(1:T, e_ref_over_time{i}, 'LineWidth', 6, 'color', 'r', 'LineStyle', '--', 'DisplayName', 'Desired');
      ylim([0 47]);
      % hold off; 
      title(sprintf('Agent %d', i));
      % set(gca, "FontSize", 45, 'FontWeight', 'bold');
      set(groot, 'defaultTextFontWeight', 'bold');
      set(groot, 'defaultAxesFontName', 'Times New Roman'); % or 'Times'
      set(groot, 'defaultTextFontName', 'Times New Roman'); % or 'Times'
      set(groot, 'defaultAxesFontSize', 65, 'defaultAxesFontWeight', 'bold'); % Match LaTeX font size (10-12pt)
      set(groot, 'defaultTextInterpreter', 'latex'); % For LaTeX symbols like $\zeta$
      % legend("FontWeight", "normal", 'FontSize', 40);
      xlabel('Time (in hr)');
      ylabel('Consumption (in kWh)', 'FontSize', 55);
      LAST_FIG = LAST_FIG + 1;
    end
    for i = 1:N 
      figure(LAST_FIG)
      plot(1:T, s_over_time{i}, 'LineWidth', 6, 'color', colors(progress, :), 'DisplayName', 'Input Trajectory');
      % yline(players{i}.ss_state(4), 'LineStyle', '--', 'LineWidth', 5, 'color', 'r', 'DisplayName', 'SteadyState');

      KZ = players{i}.K*players{i}.Z;
      xxx = KZ.extreme([0;1]); 
      yline(xxx.x(2), 'LineWidth', 6, 'color', 'r',  'LineStyle', '--', 'DisplayName', 'Steady State');
      ylim([-2 7.5]);
      % xxx = KZ.extreme([0;-1]); 
      % yline(xxx.x(2), 'LineStyle', '-.', 'LineWidth', 2, 'color', 'r');
      %   set(gca, "FontSize", 45, 'FontWeight', 'bold');
      set(groot, 'defaultTextFontWeight', 'bold');
      set(groot, 'defaultAxesFontName', 'Times New Roman'); % or 'Times'
      set(groot, 'defaultTextFontName', 'Times New Roman'); % or 'Times'
      set(groot, 'defaultAxesFontSize', 65, 'defaultAxesFontWeight', 'bold'); % Match LaTeX font size (10-12pt)
      set(groot, 'defaultTextInterpreter', 'latex'); % For LaTeX symbols like $\zeta$
      % legend("FontWeight", "normal");
      xlabel('Time (in hr)');
      ylabel('$s$ (kWh/hr)');
      title(sprintf('Agent %d',i));
      LAST_FIG = LAST_FIG + 1;
      hold on; 
    end
end