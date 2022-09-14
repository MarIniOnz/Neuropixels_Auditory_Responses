function hist4ProtocolsCluster(cl_idx, n_bins, sp, blockEvents, recParams, sec_before_trigger, type_sequence)
%
% Function that plots 4 histograms from the first 4 blocks of a type of
% sequence (for example, first 4 spaced recordings), summing up the number
% of spikes in all trials of that cluster.
%
% Inputs
%     cl_idx: num, ID of the cluster(s).
%     n_bins: num, number of bins for the histogram.
%     sp: struct, contains information about the spikes.
%     blockEvents: cell, contains all the triggers divided in their
%     different sequences.
%     recParams: struct, extracted from function recDataForSubject.
%     sec_before_trigger: num, how many seconds before the triggers we 
%     want to see to perform the t-test.
%     type_sequence: num, ID of which type of sequence you using.

% Limit of the y axes.
ylimit = 1;
ylimit2 = 1;

% Looking at how many protocols are from each type.
ISIs = recParams.ISITrial(1:4);
ISI_33 = ISIs > 0.3 & ISIs < 0.35;
ISI_10 = ISIs > 0.09 & ISIs < 0.11;

% Counts of how many protocols we have from that ISI.
num_33 = 0;
num_10 = 0;

% Initializing the axes.
ax = gobjects(sum(ISI_33),1);
ax2 = gobjects(sum(ISI_10),1);

for i=1:4
    
    % Choosig the index of the protocol and getting the data from that SOI.
    idx_chosen = i;
    [~, soi_trials, ~] = getDataFromSequence(sp, blockEvents, recParams, idx_chosen, sec_before_trigger);
    
    % Checking whether they belong to the ISI_33 or ISI_10 category.
    if ISI_33(i)
        num_33 = num_33+1;
        color='b'; % Color for ISI_33
        
        % Plotting and getting the maximum value of the histogram, for
        % plotting purposes.
        ax(num_33)= subplot(2,2,i); 
        max_value = histClusterSpaced(soi_trials.spTimes, soi_trials.spSites, cl_idx ,n_bins, recParams, idx_chosen, sec_before_trigger, type_sequence, color);
        ylimit = max(ylimit, max_value);
    else
        num_10 = num_10+1;
        color='r'; % Color for ISI_10
        
        % Plotting and getting the maximum value of the histogram, for
        % plotting purposes.
        ax2(num_10)= subplot(2,2,i);
        max_value = histClusterSpaced(soi_trials.spTimes, soi_trials.spSites, cl_idx, n_bins, recParams, idx_chosen, sec_before_trigger, type_sequence, color);
        ylimit2 = max(ylimit2, max_value);
    end
    
    % If we are at the last plot, set the limit of the axes of each
    % category to the same y limits.
    if num_33 == sum(ISI_33) && num_33 >0
        linkaxes(ax, 'y')
        ylimit = ylimit+ylimit/6; 
        ylim([0, ylimit])
    end
    
    if num_10 == sum(ISI_10) && num_10>0
        linkaxes(ax2, 'y')
        ylimit2 = ylimit2+ylimit2/6; 
        ylim([0, ylimit2])
    end 
end

