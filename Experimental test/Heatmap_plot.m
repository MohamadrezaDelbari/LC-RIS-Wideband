%--------------------------------------------------------------------------
% FINAL MATLAB SCRIPT
% - Performs background subtraction and time-gating.
% - Exports data as a 3-column vector: [Power(dB), Angle, Frequency].
% - Plots Power vs. Frequency for a specific angle.
%--------------------------------------------------------------------------
clear variables
clc
close all

% --- Parameters ---
% Time-gating window (same as Python script)
t_start_gate = 4.5e-9;  % 4.5 nanoseconds
t_end_gate   = 7.5e-9;  % 7.5 nanoseconds

% Define file angles
angles = (-40:70);

% Define paths to data folders
main_data_path = 'Data/Heatmap_folder/steeringangle_46/'; % Example for a specific steering angle
ref1_data_path = 'Data/Heatmap_folder/steeringangle_0/';   % Reference measurement 1
ref2_data_path = 'Data/Heatmap_folder/steeringangle_180/'; % Reference measurement 2

% --- Initialize Data Storage ---
freqs = [];
gated_abs_values_db = []; % This will store the dB values directly

% --- Main Loop to Process Each Angle ---
for angle = angles
    % Construct filenames for main and reference data
    main_filename = sprintf('%s/az=%d_el=0.s2p', main_data_path, angle);
    ref1_filename = sprintf('%s/az=%d_el=0.s2p', ref1_data_path, angle);
    ref2_filename = sprintf('%s/az=%d_el=0.s2p', ref2_data_path, angle);
    
    try
        % Read the main measurement and the two references
        data_main = readmatrix(main_filename, 'FileType', 'text', 'NumHeaderLines', 3, 'CommentStyle', '#');
        data_ref1 = readmatrix(ref1_filename, 'FileType', 'text', 'NumHeaderLines', 3, 'CommentStyle', '#');
        data_ref2 = readmatrix(ref2_filename, 'FileType', 'text', 'NumHeaderLines', 3, 'CommentStyle', '#');
    catch ME
        warning('Cannot read one or more files for angle %d. Skipping. Error: %s', angle, ME.message);
        continue;
    end
    
    % Store frequency vector (only once)
    if isempty(freqs)
        freqs = data_main(:, 1); % Frequency in GHz
    end
    
    % Create complex S21 data for all three measurements
    S21_main_complex = data_main(:, 4) + 1i * data_main(:, 5);
    S21_ref1_complex = data_ref1(:, 4) + 1i * data_ref1(:, 5);
    S21_ref2_complex = data_ref2(:, 4) + 1i * data_ref2(:, 5);

    % Perform background subtraction
    S21_calibrated_complex = S21_main_complex - 0.5 * (S21_ref1_complex + S21_ref2_complex);
    
    % IFFT (Frequency -> Time Domain)
    N = length(S21_calibrated_complex);
    bw = (freqs(end) - freqs(1)) * 1e9;
    time_axis = linspace(0, (N-1)/bw, N);
    S21_time = ifft(S21_calibrated_complex);
    
    % Apply the Time Gate
    gating_window = zeros(size(S21_time));
    [~, idx_start] = min(abs(time_axis - t_start_gate));
    [~, idx_end] = min(abs(time_axis - t_end_gate));
    gating_window(idx_start:idx_end) = 1;
    S21_time_gated = S21_time .* gating_window;
    
    % FFT (Time -> Cleaned Frequency Domain)
    S21_freq_gated = fft(S21_time_gated);
    
    % Convert to dB and store the column of power values
    power_db_column = 20*log10(abs(S21_freq_gated));
    gated_abs_values_db = [gated_abs_values_db, power_db_column];
end

% --- Step 1: Create the 3-column output vector ---
%%% ADDED: Reshape the data into the desired format %%%
num_freqs = length(freqs);
num_angles = length(angles);

% Create meshgrid for angles and frequencies
[angle_grid, freq_grid] = meshgrid(angles, freqs);

% Reshape the matrices into column vectors
power_vector = reshape(gated_abs_values_db, [], 1);
angle_vector = reshape(angle_grid, [], 1);
freq_vector  = reshape(freq_grid, [], 1);

% Combine into the final 3-column matrix
final_data_vector = [power_vector, angle_vector, freq_vector];

% Display the size and first few rows of the final data vector
fprintf('Successfully created final_data_vector with size %d x %d\n', size(final_data_vector, 1), size(final_data_vector, 2));
disp('First 5 rows of final_data_vector: [Power(dB), Angle(deg), Frequency(GHz)]');
disp(final_data_vector(1:5, :));


% --- Step 2: Plot Power vs. Frequency for a specific angle ---
%%% ADDED: New plot as requested %%%
target_angle = 40;

% Find the column index corresponding to the target angle
angle_index = find(angles == target_angle);

% Power + noise power for normalization
gated_abs_values_db=gated_abs_values_db+65;

if isempty(angle_index)
    warning('Target angle %d not found in the processed angles.', target_angle);
else
    % Extract the power data for the specific angle
    power_at_target_angle = gated_abs_values_db(:, angle_index);
    
    % Create the plot
    figure('Name', sprintf('Power vs. Frequency at %d degrees', target_angle));
    plot(freqs, power_at_target_angle, 'LineWidth', 1.5);
    grid on;
    xlabel('Frequency (GHz)');
    ylabel('Power (dB)');
    title(sprintf('Calibrated and Gated Power vs. Frequency at Angle = %d°', target_angle));
    legend(sprintf('Angle %d°', target_angle));
end

% --- Optional: Plot the original heatmap for visual confirmation ---
figure('Name', 'Final Heatmap of Gated |S21|');
pcolor(angles, freqs, gated_abs_values_db);
shading interp
view(0,90)
colormap(hot)
colorbar
colorbar('FontSize',15)
xlabel('Azimuth angle (degree)','FontSize',22,'FontWeight','bold');
ylabel('Frequency (GHz)','FontSize',22,'FontWeight','bold');
title('SNR (dB)','FontSize',30,'FontWeight','bold');
%ylim([53 67])
set(gca, 'YDir', 'normal');
clim([12 25]); % Adjust these limits as needed
ax.Fontsize = 30;
exportgraphics(gca,'Heatmap_0.pdf')
