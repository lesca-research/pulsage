function manualTagQRS(ecgSignal, fs)
    % ecgSignal: The ECG signal array
    % fs: Sampling frequency of the ECG signal

    figure; % Create a new figure window
    plot((1:length(ecgSignal))/fs, ecgSignal);
    title('Click on the QRS peaks, then press Enter');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    hold on;

    % User manually selects points on the plot (QRS peaks)
    [x, ~] = ginput; % Returns the x-coordinates (time) of the points clicked. Ignore y-coordinates.

    % Convert time back to sample indices since ginput returns time
    qrsPeaksIndices = round(x * fs);

    % Plot the selected points on the graph
    plot(qrsPeaksIndices/fs, ecgSignal(qrsPeaksIndices), 'ro', 'MarkerSize', 8, 'LineWidth', 2);

    % Optionally, save the selected QRS peaks indices to a file or workspace
    % save('QRS_Peaks.mat', 'qrsPeaksIndices');

    hold off;
end
