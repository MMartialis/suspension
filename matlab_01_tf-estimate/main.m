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


idata = iddata(acc_out, acc_in, delta_t);
tf_model = tfest(idata, 2, 1);


%% Plot the results

figure(1);
title('Acceleration Input vs. Output');
plot(time, acc_in, 'b', 'LineWidth', 1.5);
hold on;
plot(time, acc_out, 'r', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
legend('Input', 'Output');
grid on;
drawnow;
hold off;