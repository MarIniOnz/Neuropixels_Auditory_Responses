function Sequence = plotSpikesFrequenciesDepthsFit(subject, subject_rec, direc, imec, sp, blockEvents, sec_before_trigger, Depth1, Depth2, idx_chosen)
%
% Function that plots the mean number of spikes after a trigger with a 
% specified frequency in the sequences Random, and Random Fifty among 10
% different depths (9 subdivision between Depth1 and Depth2).
%
% Done in the indicated depths separately and plotted the mean over all
% trials, highlighting the frequencies that are repeated more than num_big
% times.
%
% Inputs:
%     subject: num, subject ID.
%     subjec_rec: num, subject ID adding a 1 or 2 for Lateral or Dorsal
%     recordings.
%     direc: directory in which the files of the first recording (a) are 
%     located for that subject.
%     imec: num, 0 or 1: imec meaning the number of probe from which we are
%     looking at the data.
%     sp: struct, contains information about the spikes.
%     blockEvents: cell, where triggers are stored.
%     sec_before_trig: num, seconds to look and integrate before and after
%     trigger.
%     Depth1: int, from which depth you want to integrate the spikes.
%     Depth2: int, to which depth you want to integrate the spikes.
%     idx_chosen: num, index of the sequence we want (1st, 2nd ... time     
%     that the sequence is performed in the recording).
%
% Outputs:
%     Sequence: struct, contains information about the frequencies and the
%     spikes post trigger.
%     Fields:
%         data: struct, data from the acquisition.
%         frequencies: double, frequencies played, sorted from lower to
%         higher.
%         num_frequencies: double, how many times was each frequency
%         actually played.
%         post_trig: double, normalized number of spikes found after the trigger.
%         (Frequency, depth, type of sequence (1= Random fifty, 2=Random)).
%         data: double, frequencies presented.
%         goodness: struct, different metrics of how good is the Gaussian
%         fit to the data (excluding the most repeated frequencies).
%         goodness_with_high_freqs: struct, different metrics of how good
%         is the Gaussian fit to the data (including the most repeated
%         frequencies).

% Calculating the different depths in which we divide the plots.
num_depths = 5; % # of depths in which we divided the limiting depths.
depths = round(Depth1:(Depth2-Depth1)/num_depths:Depth2);

% Plotting looping over the clusters.
for i = 1: num_depths

    a(i) = subplot(num_depths,1,i);

    for m=1:2   
         
        hold on
        % 2 = Random, 3 = Random fifty
        type_sequence = 1+m;
        
        % Getting the recording parameters and the data from the sequence.
        recParams = recDataForSubject(subject, subject_rec, type_sequence);
        [soi, ~, triggersInterest] = getDataFromSequence(sp, blockEvents, recParams, idx_chosen, sec_before_trigger);

        % Loading the data of that sequence (which frequencies were played
        % in each trigger time).
        Sequence(m,i).data = load("\\charite.de\\centren\\Fakultaet\\MFZ\\NWFZ\\AG-deHoz-Scratch\\Neuropixels\\" + direc.date +"\\" +recParams.folderMat{idx_chosen} + "\\data.mat"); %#ok<*SAGROW>

        % Loading the data of that sequence (which frequencies were played
        % in each trigger time).        
        Sequence(m,i).frequencies = [];
        Sequence(m,i).num_frequencies = [];
        total_frequencies = length(unique(Sequence(m,i).data.Trseq));  
        Sequence(m,i).post_trig = zeros(total_frequencies, 1);

        % Spikes pertainig in the SOI to that cluster.
        cl_spikes = soi.spTimes(soi.spDepths<depths(i+1) & soi.spDepths>depths(i))';
        
        for k = 1:length(triggersInterest)
            % Seeing which frequency was played in that trigger.
            freq = Sequence(m,i).data.Trseq(k);
            
            % Checking whether that frequency has already been played.
            if ~ismember(freq, Sequence(m,i).frequencies)
                Sequence(m,i).frequencies = [Sequence(m,i).frequencies, freq];
                idx = length(Sequence(m,i).frequencies);
                Sequence(m,i).num_frequencies = [Sequence(m,i).num_frequencies, 1];
            % If it has been played, look for the index in the sequence of
            % frequencies of this current frequency and add a unit to the
            % number of times that frequency was played.
            else 
                idx = (Sequence(m,i).frequencies == freq);
                Sequence(m,i).num_frequencies(idx) = Sequence(m,i).num_frequencies(idx)+1;
            end
            % Summing all the spikes that were played after that frequency
            % and sum it to the previous spikes if that frequency has
            % already been played.
            Sequence(m,i).post_trig(idx) = Sequence(m,i).post_trig(idx) + sum(cl_spikes>triggersInterest(k)& cl_spikes<(triggersInterest(k)+sec_before_trigger));        
        end

        % Normalize the number of spikes by the number of times that the
        % frequency was played.
        for n = 1:length(Sequence(m,i).num_frequencies)
            Sequence(m,i).post_trig(n) = Sequence(m,i).post_trig(n)/Sequence(m,i).num_frequencies(n);
        end
        
        [Sequence(m,i).frequencies, order]=sort(Sequence(m,i).frequencies);
        Sequence(m,i).num_frequencies = Sequence(m,i).num_frequencies(order);
        Sequence(m,i).post_trig = Sequence(m,i).post_trig(order)';
        
        % Size of the markers.
        small_size=10;
        big_size=70;
        
        % How many times a frequency has to be played to be considered a
        % "common one".
        num_big = 20;
        
        if m==1
            % Checking which frequencies were "common".
            big_freq_1 = Sequence(m,i).num_frequencies > num_big;
            small_freq_1 = Sequence(m,i).num_frequencies <= num_big;
            % The frequencies that were "common", plot them in a bigger
            % size.
            s(1) = scatter(log10(Sequence(m,i).frequencies(small_freq_1)), Sequence(m,i).post_trig(small_freq_1), small_size, 'ro'); hold on
            try
                [f1, Sequence(m,i).goodness] = fit(log10(Sequence(m,i).frequencies(small_freq_1))',Sequence(m,i).post_trig(small_freq_1)','gauss1'); % Gaussian fitting
                [~, Sequence(m,i).goodness_with_high_freq]= fit((Sequence(m,i).frequencies)',Sequence(m,i).post_trig','gauss1'); % Gaussian fitting with high frequencies included
                plot(f1,'Red'); hold on
            catch
            end
          else
            % Checking which frequencies were "common".
            big_freq_2 = Sequence(m,i).num_frequencies > num_big;
            small_freq_2 = Sequence(m,i).num_frequencies <= num_big;
            % The frequencies that were "common", plot them in a bigger
            s(3) = scatter(log10(Sequence(m,i).frequencies(big_freq_2)), Sequence(m,i).post_trig(big_freq_2), big_size, 'bo'); hold on
            try
                [f2, Sequence(m,i).goodness] = fit(log10(Sequence(m,i).frequencies(small_freq_2))',Sequence(m,i).post_trig(small_freq_2)','gauss1'); % Gaussian fitting
                [~, Sequence(m,i).goodness_with_high_freq]= fit((Sequence(m,i).frequencies)',Sequence(m,i).post_trig','gauss1'); % Gaussian fitting with high frequencies included
                plot(f2,'Blue'); hold on
            catch
            end
            s(2) = scatter(log10(Sequence(1,i).frequencies(big_freq_1)), Sequence(1,i).post_trig(big_freq_1), big_size, 'ro', 'DisplayName', ''); hold on
            s(4) = scatter(log10(Sequence(m,i).frequencies(small_freq_2)), Sequence(m,i).post_trig(small_freq_2), small_size, 'bo', 'DisplayName', ''); hold on
            
        end

    end
    
    % Annotation
    try
        text1="R2 Random: "+round(Sequence(1,i).goodness.rsquare,2)+newline+"R2 Random with high freqs: "+round(Sequence(1,i).goodness_with_high_freq.rsquare,2)+newline+"R2 Random Fifty: "+round(Sequence(2,i).goodness.rsquare,2)+newline+"R2 Random Fifty with high freqs: "+round(Sequence(2,i).goodness_with_high_freq.rsquare,2);
        b = annotation('textbox','String',text1,'Position',a(i).Position,'Vert','top','FitBoxToText','on');
    catch
    end
    
    % Set the scale, title, and legend.
    title("Spikes depth " + num2str(depths(i)) + " to "+ num2str(depths(i+1))+ ". "+ string(idx_chosen) + "° time presented." );
%     set(gca,'Xscale',  'log'), hold off
    
    % Limists and labelling of the axes.
    xlabel('Frequencies');
    ylabel('# of spikes (Normalized)')
    % xlim([2400 40000])
%     xlim([min(Sequence(m,i).frequencies)*0.7 max(Sequence(m,i).frequencies)*1.3])
    
    % Y limits.
    if sum(sum(Sequence(m,i).post_trig>0))>0
         ylim([0, max(max(Sequence(m,i).post_trig)*1.5)])
    else
        ylim([0, 5]);
    end 
    
   % Legend only on the last depth.
    if i == num_depths
     legend('Random', 'Fitted Random', 'Random fifty','Fitted Random Fifty','Location','northeast')
    else
        legend('off')
    end
   
end

sgtitle(subject_rec +", " +imec +", Random vs Random fifty response.")