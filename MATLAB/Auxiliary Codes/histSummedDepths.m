function max_value = histSummedDepths(spTimes, spDepths, Depth1, Depth2, n_bins, recParams, idx_chosen, sec_before_trigger, type_seq, color)
%
% Function that represents a histogram with all the spikes occured among
% two depths in a sequence, and the 8 different triggers of every trial.
%
% Inputs
%     spTimes: cell, with the vectors of the # trials containing the times 
%     of the spikes per trial.
%     spDephts: cell, with the vectors of the # trials containing the
%     depths of the clusters to which the spikes belong.
%     Depth1: int, from which depth you want to integrate the spikes.
%     Depth2: int, to which depth you want to integrate the spikes.
%     n_bins: num, number of bins for the histogram.
%     recParams: struct, extracted from function recDataForSubject.
%     idx_chosen: num, index of the sequence we want (1st, 2nd ... time     
%     that the sequence is performed in the recording).
%     sec_before_trigger: num, how many seconds before the triggers we want
%     to see to perform the t-test.
%     type_seq: num, ID of which type of sequence you using.
%     color: vector, RGB values. 
%
% Outputs
%     max_value: double, number of maximum histogram count, for plotting
%     purposes.

% Get the time between triggers
ISI = recParams.ISITrial(idx_chosen);
triggers = 0: ISI: (ISI*7+ISI/2); % we get 8 Triggers this way

spikes = [];

% See how many spikes between those Depths per Trial.
for iTrial = 1:length(spTimes)
    
    trial_sp = spTimes{iTrial}(spDepths{iTrial}<Depth2 & spDepths{iTrial}>Depth1)';
    spikes = [spikes, trial_sp]; %#ok<AGROW>
end

% Performing the histogram and getting the maximum count for plotting.
hist_values = histogram(spikes, n_bins, 'FaceColor', color);
max_value = max(hist_values.Values);

% Plot the trigger times in red vertical lines.
for i=1:length(triggers)
    xline(triggers(i),'blue')
end

% Labeling of the axes
xlabel('Time [s]')
ylabel('# spikes ');

xlim([-sec_before_trigger*(ISI*10),ISI*8])

% Create an appropiate title.
title("Spikes depth " + num2str(Depth1) + " to "+ num2str(Depth2)+ ". "+ string(idx_chosen) + "° " + typeOfSequence(type_seq) + " " + ISI);
