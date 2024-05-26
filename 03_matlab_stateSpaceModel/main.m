%clear cache and variables
clc;clear

% Open a file dialog window for selecting the CSV file
[filename, filepath] = uigetfile('.csv', 'Select CSV file', '../logs/');

% Check if user canceled file selection
if isequal(filename,0) || isequal(filepath,0)
    disp('File selection canceled.');
    return;
end

% Load the CSV file
fullfile_path = fullfile(filepath, filename);
fileID = fopen(fullfile_path);
headerLine = fgetl(fileID); % Read the first line (header line) of the file
fclose(fileID);

% Parse the header line to get signal names and find the TIME column
header_cells = strsplit(headerLine, ',');
time_column_index = find(strcmp(header_cells, 'TIME'));
valve_column_index = find(strcmp(header_cells, 'VALVE'));
signal_column_indexes = setdiff(1:1:length(header_cells),[time_column_index valve_column_index]);

% Load the CSV file excluding the first row (header line)
data = readmatrix(fullfile_path);

delta_t = 0.01; % Average time interval

% Extract TIME and Signal columns
original_time = data(:, time_column_index);
% time = (min(original_time)+.5:delta_t:(max(original_time))-4)';
valve = data(:, valve_column_index);
signals = data(:, signal_column_indexes);

nan_rows = isnan(valve);
valve_time = original_time(~nan_rows);
time = (min(valve_time)+.5:delta_t:(max(valve_time))-2 )';

valve = interp1(original_time(~nan_rows), valve(~nan_rows),time, "spline");



% Calculate sampling frequency based on the entire time series
fs = 1 / delta_t; % Sampling frequency in Hz

% Design a band-pass filter (Butterworth filter)
cutoff_frequency_high = 10; % Cutoff frequency in Hz
cutoff_frequency_low = 1; % Cutoff frequency in Hz

order = 2; % Filter order
[b, a] = butter(order, [(cutoff_frequency_low / (fs/2)) (cutoff_frequency_high / (fs/2))]);

% Apply the filter to each signal
filtered_signals = zeros(size(time, 1), size(signals, 2));
for i = 1:size(signals, 2)
    nan_rows = isnan(signals(:, i));
    signal_time = original_time(~nan_rows);
    signal_time = (signal_time - min(signal_time) + min(valve_time)) * (max(valve_time) - min(valve_time)) / (max(signal_time) - min(signal_time));
    filterable = interp1(signal_time, signals(~nan_rows, i),time, "spline");
    
    filtered_signals(:, i) = filtfilt(b, a, filterable).*9.81;
end

% Shift time so it starts from 0
time = time - min(time);

acc_in = filtered_signals(:,4)-mean(filtered_signals(:,4), 1);
acc_out = filtered_signals(:,1)-mean(filtered_signals(:,1), 1);

% Velocity from accelaration
in  = cumtrapz(time, acc_in);
out = cumtrapz(time, acc_out);

% Avoiding drift
% apply a high-pass filter to the input signal and the output signal
% Design a high-pass filter (Butterworth filter)
cutoff_frequency = 0.5; % Cutoff frequency in Hz
order = 2; % Filter order
[b, a] = butter(order, cutoff_frequency / (fs/2), 'high');

% % Apply the filter to each signal
in = filtfilt(b, a, in);
in = in-mean(in);
out = filtfilt(b, a, out);
out = out-mean(out);

figure(2);
clf;
hold on;
plot(time, out, 'LineWidth', 1.5);
plot(time, in,  'LineWidth', 1.5);
plot(time, valve);
% plot(time, acc_in);
% plot(time, acc_out);
xlabel('Time');
ylabel('Velocity:m/s');
title('Filtered Signal');
legend('Output', 'Input', 'Valve', 'ACC in', 'ACC out');
grid on;
hold off;
drawnow;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state space model estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = [valve cumtrapz(in) in]; % Input signal (valve, pos, vel)
var = 0.1; % Variance for mutation
G = 50; % Number of generations (increase for better results)
G_2 = 10; % Number of generations inside one learning step
P = 5000; % Population size  
p = 2500; % Number of top mutations to keep


% m1, m2, b1, b2, b3, k1, k2
initial_guess = [0.1 0.03604634 -0.1 -0.1 7.725485  0.1 831.7944]; % Initial guess for the parameters
% initial_guess = [0 3 0 0  8 0 300]; % Initial guess for the parameters
params_b = repmat(initial_guess, P, 1); % Initialize the parameters
error_b = zeros(1, P) + inf; % Initialize the error
best_mutations = 0;
best_outs = 0;
bests_error = 0;

for g = 1:G
    tic
    [best_mutations, best_outs, bests_error] = learn(time, u, params_b, out, var, G_2, P, p, error_b, delta_t);
    toc
    
    params_b = repmat(best_mutations, P/p, 1);
    error_b = repmat(bests_error, 1, P/p);
    disp(['Generation ' num2str(g) ' completed. Best error: ' num2str(min(bests_error)) ' Worst error: ' num2str(max(bests_error))]);
    disp(['Best mutations: ' num2str(best_mutations(1, :))]);
    
    % clearvars -except time u out var G G_2 P p big_variance best_mutations bests_out bests_error params_b error_b
    
    
    figure(3);
    clf;
    hold on;
    % Plot the outputs
    plot(time, out, 'LineWidth', 1.5);
    plot(time, in, 'LineWidth', 1.5);
    for i = 1:5
        plot(time, best_outs(:, i));
    end
    title('Best Mutations');
    xlabel('Time');
    ylabel('ACC:m/s^2');
    legend('output', 'input', 'Mutation 1', 'Mutation 2', 'Mutation 3', 'Mutation 4', 'Mutation 5');
    grid on;
    drawnow;
    hold off;
end