function plot_chemical_states(time, x_hist, plotTitle, timeType)
    % Indices: H1(1), H2(5), H3(9), XB3(11)
    indices = [1, 5, 9, 11]; 
    titles = {'H1 (Reactor 1 Level) Deviation', 'H2 (Reactor 2 Level) Deviation', ...
              'H3 (Separator Level) Deviation', 'XB3 (Product Concentration) Deviation'};
    y_labels = {'Level Deviation [m]', 'Level Deviation [m]', 'Level Deviation [m]', 'Conc. Deviation [wt.%]'};
    
    % Steady-State (Nominal) Values provided by the user
    h_ss = [29.8, 30.0, 3.27]; 

    figure('Name', plotTitle, 'Color', 'w', 'Position', [150 150 1100 750]);
    
    for k = 1:length(indices)
        subplot(2, 2, k);
        signal = x_hist(indices(k), :);
        final_val = signal(end); % Deviation final value (usually 0)
        
        % --- PERFORMANCE CALCULATIONS ---
        % 1. Settling Time (2% of peak deviation from deviation signal)
        peak_dev_total = max(abs(signal - final_val));
        if peak_dev_total < 1e-6
            ts_val = 0;
        else
            idx_settle = find(abs(signal - final_val) > 0.02 * peak_dev_total, 1, 'last');
            if isempty(idx_settle)
                ts_val = 0;
            else
                ts_val = time(idx_settle);
            end
        end
        
        % 2. Overshoot (OS%) Calculation relative to Steady-State
        if k <= 3 
            max_v = max(signal);
            min_v = min(signal);
            
            % Find the maximum deviation from the equilibrium (0)
            if signal(1) > final_val
                % If starting positive, we look for the lowest dip below final_val
                actual_peak_deviation = abs(min_v - final_val);
            else
                % If starting negative, we look for the highest peak above final_val
                actual_peak_deviation = abs(max_v - final_val);
            end
            
            % Calculate OS% relative to PHYSICAL STEADY STATE
            % Formula: (Deviation / Nominal Value) * 100
            overshoot = (actual_peak_deviation / h_ss(k)) * 100;
            
            % Format string with 3 decimal places as OS will be very small
            header_str = sprintf('%s\nTs = %.2fs | OS_{ss} = %.3f%%', titles{k}, ts_val, overshoot);
        else
            % Title with ONLY Ts for XB3
            header_str = sprintf('%s\nTs = %.2fs', titles{k}, ts_val);
        end
        
        % --- PLOTTING ---
        if strcmpi(timeType, 'Discrete')
            stairs(time, signal, 'LineWidth', 1.5, 'Color', [0.8 0.1 0.1]);
        else
            plot(time, signal, 'LineWidth', 1.5, 'Color', [0.1 0.2 0.8]);
        end
        
        title(header_str);
        grid on;
        xlabel('Time (s)');
        ylabel(y_labels{k});
        
        % Reference line for equilibrium (0 in deviation variables)
        hold on;
        line([time(1) time(end)], [final_val final_val], 'Color', 'k', 'LineStyle', ':');
    end
    sgtitle([plotTitle ' - ' timeType ' Simulation']);
end