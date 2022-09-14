function max_value = histClusterSpaced(spTimes, spSites, cl_idx, n_bins, recParams, idx_chosen, sec_before_trigger, type_seq, color)
%
% Function that represents a histogram with all the spikes accuring in a 
% in a specified block of a sequence, along with the triggers as vertical
% lines.
%
% Inputs
%     spTimes: cell, with the vectors of the # trials containing the times 
%     of the spikes per trial.
%     spSites: cell, with the vectors of the # trials containing the
%     IDs of the clusters to which the spikes belong.
%     cl_idx: num, cluster ID.
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
    
    trial_sp = spTimes{iTrial}(spSites{iTrial}==cl_idx)';
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
title("Spike cluster: " + num2str(cl_idx) + ". "+ string(idx_chosen) + "° " + typeOfSequence(type_seq) + " " + ISI);
