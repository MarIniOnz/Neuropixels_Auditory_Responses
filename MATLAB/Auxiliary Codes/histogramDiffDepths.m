function histogramDiffDepths(Depth1, Depth2, n_bins, sp, blockEvents, recParams, sec_before_trigger, type_sequence, num_sequences)
%
% Function that plots x histograms from the first x blocks of a type of
% sequence (for example, first x spaced recordings) for each one of the k
% specified depths. For instance, k = 9 --> x*9 plots. They will show
% the sum of the number of spikes in all trials between specified depths.
%
% Inputs
%     Depth1: int, from which depth you want to integrate the spikes.
%     Depth2: int, to which depth you want to integrate the spikes.
%     n_bins: num, number of bins for the histogram.
%     sp: struct, contains information about the spikes.
%     blockEvents: cell, contains all the triggers divided in their
%     different sequences.
%     recParams: struct, extracted from function recDataForSubject.
%     sec_before_trigger: num, how many seconds before the triggers we want
%     to see to perform the t-test.
%     type_sequence: num, ID of which type of sequence you using.
%     num_sequences: num, number of sequences of that type you want to show.

% Calculating the different depths in which we divide the plots.
k = 9; % # of depths in which we divided the limiting depths.
depths = round(Depth1:(Depth2-Depth1)/k:Depth2);

% Looking at how many protocols are from each type.
ISIs = recParams.ISITrial(1:num_sequences);
ISI_33 = ISIs > 0.3 & ISIs < 0.35;
ISI_10 = ISIs > 0.09 & ISIs < 0.11;

% Initializing the axes.
ax = gobjects(sum(ISI_33),k);
ax2 = gobjects(sum(ISI_10),k);

% Color coding for the histograms
color_hex_33=["d7e1ee", "cbd6e4", "bfcbdb", "b3bfd1", "a4a2a8", "df8879", "c86558", "b04238", "991f17"];
color_hex_10=["0000b3", "0010d9", "0020ff", "0040ff", "0060ff", "0080ff", "009fff", "00bfff", "00ffff"];

for i=1:k % k Times
    
    % Counts of how many protocols we have from that ISI.
    num_33 = 0;
    num_10 = 0;
    
    % Limit of the y axes.
    ylimit = 1;
    ylimit2 = 1;
    
    for j=1:num_sequences % x Spaced sequences
        
        % Choosing the index of the protocol and getting the data from that SOI.
    	idx_chosen = j;
        [~, soi_trials, ~] = getDataFromSequence(sp, blockEvents, recParams, idx_chosen, sec_before_trigger);

        % Checking whether they belong to the ISI_33 or ISI_10 category.
        if ISI_33(idx_chosen)
            num_33 = num_33+1;
            % Selecting the color for each depth.
            color = sscanf(color_hex_33(i),'%2x%2x%2x',[1 3])/255;
            
            % Plotting and getting the maximum value of the histogram, for
            % plotting purposes.
            ax(num_33,k)= subplot(k,num_sequences,j+num_sequences*(i-1));
            max_value = histSummedDepths(soi_trials.spTimes, soi_trials.spDepths, depths(i), depths(i+1), n_bins, recParams, idx_chosen, sec_before_trigger, type_sequence, color);
            ylimit = max(ylimit, max_value);
        else
            num_10 = num_10+1;
            % Selecting the color for each depth.
            color = sscanf(color_hex_10(i),'%2x%2x%2x',[1 3])/255;

            % Plotting and getting the maximum value of the histogram, for
            % plotting purposes.
            ax2(num_10,k)= subplot(k,num_sequences,j+num_sequences*(i-1));
            max_value = histSummedDepths(soi_trials.spTimes, soi_trials.spDepths, depths(i), depths(i+1), n_bins, recParams, idx_chosen, sec_before_trigger, type_sequence, color);
            ylimit2 = max(ylimit2, max_value);
        end

        % If we are at the last plot, set the limit of the axes of each
        % category to the same y limits.
        if num_33 == sum(ISI_33) && num_33>0 && ISI_33(idx_chosen)
            linkaxes(ax(:,k), 'y')
            ylimit = ylimit+ylimit/6; 
            ylim([0, ylimit])
        end

        if num_10 == sum(ISI_10) && num_10>0 && ISI_10(idx_chosen)
            linkaxes(ax2(:,k), 'y')
            ylimit2 = ylimit2+ylimit2/6; 
            ylim([0, ylimit2])
        end 
    end
end

