function plotSpikesFrequenciesSummed(subject, subject_rec, imec, sp, blockEvents, sec_before_trigger, cl_id, idx_chosen)
% 
% Function that plots the mean number of spikes after a trigger with a 
% specified frequency in the sequences Random, and Random Fifty. 
% 
% Performs a mean over all clusters specified and puts them all in the same
% figure.
%
% Inputs:
%     subject: num, subject ID.
%     subjec_rec: num, subject ID adding a 1 or 2 for Lateral or Dorsal
%     recordings.
%     imec: num, 0 or 1: imec meaning the number of probe from which we are
%     looking at the data.
%     sp: struct, contains information about the spikes.
%     blockEvents: cell, where triggers are stored.
%     sec_before_trig: num, seconds to look and integrate before and after
%     trigger.
%     cl_id: num, cluster ID(s).
%     medians: double, vector containing the median response frequency of
%     each cluster, it can be an empty vector of equal length to cl_id.
%     idx_chosen: num, index of the sequence we want (1st, 2nd ... time     
%     that the sequence is performed in the recording).

num_cl = length(cl_id);

% Plotting looping over the clusters.
for m=1:2    
    
    hold on
    
    % Looping over all clusters.
    for i = 1: num_cl 
        
        % 2 = Random, 3 = Random fifty
        type_sequence = m+1;
        
        % Getting the recording parameters and the data from the sequence.
        [direc, ~ , ~ ] = filesForSubject(subject_rec, imec);
        recParams = recDataForSubject(subject, subject_rec, type_sequence);
        [soi, ~, triggersInterest] = getDataFromSequence(sp, blockEvents, recParams, idx_chosen, sec_before_trigger);

        % Loading the data of that sequence (which frequencies were played
        % in each trigger time).
        Sequence(m).data = load("\\charite.de\\centren\\Fakultaet\\MFZ\\NWFZ\\AG-deHoz-Scratch\\Neuropixels\\" + direc.date +"\\" +recParams.folderMat{idx_chosen} + "\\data.mat"); %#ok<*SAGROW>
        
        % Loading the data of that sequence (which frequencies were played
        % in each trigger time).
        Sequence(m).frequencies = [];
        Sequence(m).num_frequencies = [];
        total_frequencies = length(unique(Sequence(m).data.Trseq));
        post_trig=zeros(total_frequencies,num_cl);    
        
        % Spikes pertainig in the SOI to that cluster.
        cl_spikes = soi.spTimes(soi.spSites==cl_id(i))';
        
        for k = 1:length(triggersInterest)
            % Seeing which frequency was played in that trigger.
            freq = Sequence(m).data.Trseq(k);

            % Checking whether that frequency has already been played.
            if ~ismember(freq, Sequence(m).frequencies)
                Sequence(m).frequencies = [Sequence(m).frequencies, freq];
                idx = length(Sequence(m).frequencies);
                Sequence(m).num_frequencies = [Sequence(m).num_frequencies, 1];
            % If it has been played, look for the index in the sequence of
            % frequencies of this current frequency and add a unit to the
            % number of times that frequency was played.
            else 
                idx = (Sequence(m).frequencies == freq);
                Sequence(m).num_frequencies(idx) = Sequence(m).num_frequencies(idx)+1;
            end
            % Summing all the spikes that were played after that frequency
            % and sum it to the previous spikes if that frequency has
            % already been played.
            post_trig(idx) = post_trig(idx)+sum(cl_spikes>triggersInterest(k)& cl_spikes<(triggersInterest(k)+sec_before_trigger));        
        end
            
        % Normalize the number of spikes by the number of times that the
        % frequency was played.
        for n = 1:length(Sequence(m).num_frequencies)
            post_trig(n,i) = post_trig(n,i)/Sequence(m).num_frequencies(n);
        end
    end
    
    % Mean over all clusters.
    plot_spikes = sum(post_trig,2)/num_cl;
    
    % Plot them in different colors.
    if m==1
            scatter(Sequence(m).frequencies, plot_spikes, 'ro', 'DisplayName', 'Random')
        else
            scatter(Sequence(m).frequencies, plot_spikes, 'bo', 'DisplayName', 'Random fifty')
    end
end

% Set a logarithmic scale.
set(gca, 'Xscale',  'log')

% Limits and labelling of the axes.
xlabel('Frequencies');
ylabel('# of spikes (Normalized)')
xlim([min(Sequence(m).frequencies)*0.8 max(Sequence(m).frequencies)*1.2])

% Y limits.
if sum(sum(post_trig>0))>0
     ylim([0, max(max(post_trig)*1.3)])
else
    ylim([0, 5]);
end   

legend()
title(subject_rec +", " +imec +", Random vs Random fifty response.")
