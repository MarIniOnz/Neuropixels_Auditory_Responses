function plotSpikesFrequenciesAllTrials(subject, subject_rec, direc, myKsDir, imec, sp, blockEvents, sec_before_trigger, cl_id, medians, idx_chosen)
% 
% Function that plots the number of spikes after a trigger with a 
% specified frequency in the sequences Random, and Random Fifty. 
% 
% Done in the indicated clusters separately and plots ALL TRIALS,
% highlighting the frequencies that are repeated more than num_big times.
%
% Inputs:
%     subject: num, subject ID.
%     subjec_rec: num, subject ID adding a 1 or 2 for Lateral or Dorsal
%     recordings.
%     direc: directory in which the files of the first recording (a) are 
%     located for that subject.
%     myKsDir: directory where the whole recording in contained for
%     Kylosort.
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

% Figure out the disposition of the graph in rows and cols.
num_cl = length(cl_id);
rows = floor(sqrt(num_cl));
cols = floor(num_cl/rows) + sum(mod(num_cl,rows)~=0);

% Plotting looping over the clusters.
for i = 1: length(cl_id)

    subplot(rows,cols,i)

    for m=1:2   
         
        hold on
        % 2 = Random, 3 = Random fifty
        type_sequence = 4-m;
        
        % Getting the recording parameters and the data from the sequence.
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
        post_trig=cell(total_frequencies,1); 
        
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
            % and storing it in different vectors depending on the trial.
            post_trig{idx}(end+1) = sum(cl_spikes>triggersInterest(k)& cl_spikes<(triggersInterest(k)+sec_before_trigger));        
        end
        
        % Size of the markers.
        small_size=5;
        big_size=20;
        
        % How many times a frequency has to be played to be considered a
        % "common one".
        num_big = 15;
        
        % Plotting all the different trials for each kind of sequence
        % differently.
        if m==1
            for q=1:length(Sequence(m).num_frequencies)
                hold on
                
                % Adding noise to the frequencies so we can see the dots in
                % the plot that have the same frequency.
                frequency = Sequence(m).frequencies(q)*ones(length(post_trig{q}),1);
                noise = 10^(floor(log10(Sequence(m).frequencies(q)))-1) * rand(length(post_trig{q}),1)*2;
                frequency = frequency + noise;
                
                % The frequencies that were "common", plot them in a bigger
                % size.
                if Sequence(m).num_frequencies(q) > num_big
                    s(q) = scatter(frequency, post_trig{q}, big_size, 'bo', 'DisplayName', 'Random fifty');
                else
                    s(q) = scatter(frequency, post_trig{q}, small_size, 'bo', 'DisplayName', '');
                end
            end
        else
            for q=1:length(Sequence(m).num_frequencies)
                hold on
                
                frequency = Sequence(m).frequencies(q)*ones(length(post_trig{q}),1);
                noise = 10^(floor(log10(Sequence(m).frequencies(q)))-1) * rand(length(post_trig{q}),1)*2; % So we can see them
                frequency = frequency + noise;
                if Sequence(m).num_frequencies(q) > num_big
                    s(q) = scatter(frequency, post_trig{q}, big_size, 'ro', 'DisplayName', 'Random');
                else
                    s(q) = scatter(frequency, post_trig{q}, small_size, 'ro', 'DisplayName', '');
                end
            end
        end
    end
    
    % Set the scale, title, and legend.
    [amp, good] = findAmpAndSorting(myKsDir, cl_id(i));
    title(cl_id(i) + " ," + amp +" , " + good)
    set(gca, 'Xscale',  'log'), hold off
    
    % Plot the median
    xline(medians(i),'--g', 'DisplayName', 'Median');
    
    % Limists and labelling of the axes.
    xlabel('Frequencies');
    ylabel('# of spikes (Normalized)')
    % xlim([2400 40000])
    xlim([min(Sequence(m).frequencies)*0.7 max(Sequence(m).frequencies)*1.3])
    
    % Y limits
    if sum(sum(post_trig>0))>0
         ylim([0, max(max(post_trig)*1.5)])
    else
        ylim([0, 5]);
    end   
    
   % Legend only on the last cluster.
   if i == num_cl
     legend('Location','northwest')
   end
   
end

sgtitle(subject_rec +", " +imec +", Random vs Random fifty response.")