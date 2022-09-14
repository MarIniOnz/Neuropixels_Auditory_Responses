function rasterplot(spTimes, spSites, recParams, idx_chosen, cluster_chosen, type_seq)
%
% Function that represents a raster plot of the trials done in a sequence
% along with the 8 different trigger times.
%
% Inputs
%     spTimes: cell, with the vectors of the # trials containing the times 
%     of the spikes per trial. Should come from soi_trials.
%     spSites: cell, with the vectors of the # trials containing the ID's 
%     of the clusters to which the spikes belong. Should come from 
%     soi_trials.
%     recParams: struct, extracted from function recDataForSubject.
%     idx_chosen: num, which # of sequence in the block was chosen.
%     cluster_chosen: num, ID of the chosen cluster.
%     type_seq: num, ID of the sequence. 
    
hold on

% Get the time between triggers
ISI = recParams.ISITrial(idx_chosen);
triggers = 0: ISI: (ISI*7+ISI/2); % we get 8 Triggers this way

% For all trials...
for iTrial = 1:length(spTimes)
    
    % Check which spikes pertain to that cluster
    trial_sp = spTimes{iTrial}(spSites{iTrial}==cluster_chosen)';
    spks            = trial_sp;         % Get all spikes of respective trial    
    xspikes         = repmat(spks,3,1);         % Replicate array
    yspikes      	= nan(size(xspikes));       % NaN array
    
    if ~isempty(yspikes)
        yspikes(1,:) = iTrial-1;                % Y-offset for raster plot
        yspikes(2,:) = iTrial;
    end
    
    plot(xspikes, yspikes, 'Color', 'k')
end

% Plot the trigger times in red vertical lines.
for i=1:length(triggers)
    xline(triggers(i),'red')
end

% Setting the limits to the axes.
xlim([-ISI/3 ISI*8]);
ylim([0 length(spTimes)]);

% Proper labeling of the axes.
xlabel('Time [s]')
ylabel('Trials');

% Create an appropiate title.
title("Cluster chosen = " + string(cluster_chosen) + ". " + string(idx_chosen) + "° " + typeOfSequence(type_seq));
