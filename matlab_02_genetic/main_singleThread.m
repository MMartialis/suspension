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


%%%%%%%%%%%%%%%%%%
% stupid algorithm
%%%%%%%%%%%%%%%%%%

theta_b = [12 19 500]; %m D k
offset_b = -10;
Error_b = inf;
var_1 = 0.1;

G = 10000;

for g = 1:G
    theta_g = normrnd(theta_b, var_1*theta_b*(1-(g/G)));
    offset_g = normrnd(offset_b, var_1*abs(offset_b)*(1-(g/G)));

    sys1 = tf([theta_g(2) theta_g(3)],[theta_g(1) theta_g(2) theta_g(3)]);

    A_out = lsim(sys1,in,time)+offset_g;

    Error_g = sum((A_out-out).^2);
    if Error_g < Error_b
        theta_b = theta_g;
        offset_b = offset_g;
        Error_b = Error_g;

        figure(3)
        plot(time,out,time,A_out)
        pause(0.1)

    end

end


%%
%%%%%%%%%%%%%%%%%%
% clever algorithm
%%%%%%%%%%%%%%%%%%

% theta_bs = zeros(100, 4); % Initialize theta_bs matrix
% for i = 1:100
%     theta_bs(i, :) = [randi([1, 20]), -randi([15, 25]), randi([400, 800]), randi([1,25])]; % Generate random values for m, D, k, offset and error
% end

% error_bs = zeros(100, 1)+inf; % Initialize error_bs matrix

% G_1 = 100; % Number of generations
% var_1 = 0.1; % Variance
% G_2 = 100; % Number of generations
% var_2 = 0.2; % Variance

% for g_1 = 1:G_1
%     parfor i = 1:size(theta_bs, 1)
%         for g_2 = 1:G_2
%             theta_g = normrnd(theta_bs(i, :), var_1*abs(theta_bs(i, :))*(1-(g_2/G_2)));

%             sys1 = tf([theta_g(2) theta_g(3)],[theta_g(1) theta_g(2) theta_g(3)]);

%             A_out = lsim(sys1,in,time)+theta_g(4);

%             error_g = sum((A_out-out).^2);
%             if error_g < error_bs(i)
%                 theta_bs(i, :) = theta_g;
%                 error_bs(i) = error_g;
%             end
%         end
%     end
%     % Sort the error_bs matrix in ascending orders
%     [sorted_errors, sorted_indices] = sort(error_bs);

%     % Get the five best mutations
%     best_mutations = theta_bs(sorted_indices(1:5), :);

%     % Plot the best mutations
%     for i = 1:5
%         figure(66);
%         clf;
%         hold on;
%         theta_g = best_mutations(i, :);
%         sys1 = tf([theta_g(2) theta_g(3)],[theta_g(1) theta_g(2) theta_g(3)]);
%         A_out = lsim(sys1,in,time)+theta_g(4);
%         plot(time, A_out);
%     end
%     title('Best Mutations');00
%     xlabel('Time');
%     ylabel('ACC:m/s^2');
%     legend('Mutation 1', 'Mutation 2', 'Mutation 3', 'Mutation 4', 'Mutation 5');
%     grid on;

%     % Update the theta_bs matrix with the best mutations
%     theta_bs = repmat(best_mutations, 20, 1) + var_2*randn(100, 4);
% end
