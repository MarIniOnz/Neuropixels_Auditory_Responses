% Very outdated code, can be used to extract small bits of code at the end
% that are commented and might be useful.

%% Restart

clc
clear

%% Add the repositories to your path

addpath(genpath('\\charite.de\centren\Fakultaet\MFZ\NWFZ\AG-deHoz-Scratch\Kilosort3')) % path to kilosort folder
addpath('\\charite.de\centren\Fakultaet\MFZ\NWFZ\AG-deHoz-Scratch\Kilosort3\npy-matlab-master\npy-matlab') % for converting to Phy
addpath('\\charite.de\centren\Fakultaet\MFZ\NWFZ\AG-deHoz-Scratch\Kilosort3\spikes-master\visualization')
addpath(genpath('\\charite.de\centren\Fakultaet\MFZ\NWFZ\AGdeHoz\Martin'))


%% Subjects defined and so on

subjects = [607, 608, 609, 611, 6091, 6092, 6111, 602, 603, 6021, 6031];

Depths_imec0 = [1400, 1900; 1500, 1900; 2000, 3500; 2100, 2500; 2000, 3000; 2000, 2800; 2100, 2400; 2300, 2800; 1300, 2500; 800, 1500; 1250, 2300];
Depths_imec1 = [2500, 3200; 2200, 3200; 2000, 3000; 2100, 2800; 2000, 2800; 1500, 2000; 2400, 2570; 0, 0; 0, 0; 0, 0; 0, 0];

%% FOR loop
% for k = 3
for k = 1:length(subjects)
    %% Get the directories and the recording parameters

    subject_rec = subjects(k);
    
    for m = 1:2
        imec = m-1;
        
        if imec==0
            Depth1 = Depths_imec0(k,1);
            Depth2 = Depths_imec0(k,2);
        else
            Depth1 = Depths_imec1(k,1);
            if Depth1 == 0
                continue
            end
            Depth2 = Depths_imec1(k,2);
        end

        [direc, myKsDir, num_sequences] = filesForSubject(subject_rec, imec);

        % In case we use b, c sequences.
        if subject_rec>1000
            subject=floor(subjects(k)/10);
        else 
            subject = subject_rec;
        end

        %% Getting the spike times and the clusters they belong to

        [sp.spikeTimes, sp.spikeAmps, sp.spikeDepths, sp.spikeSites] = ksDriftmap(myKsDir,1);
        clusters_with_spikes = unique(sp.spikeSites);

        %% Extracting the trigger events of the sequence of our block of interest

        blockEvents = separateTriggersIntoBlocks(direc.folder, num_sequences);
        sec_before_trigger = 0.03; % Default: 0.03 s.

        idx_chosen = 2; % Which time of the block we want the sequence from (1st time we have a spaced, 2nd time we present it...)

        %% Selecting the type of sequence we want and the recording parameters for that type of sequence
        % Select which kind of sequence you want to look at: Legend
            % 1. FRA. 2. Fixed random. 3. Random fifty. 4. Random fifty silence.
            % 5. Fixed fifty silence. 6. Fixed spaced. 7. Fixed compressed.
            % 8. Completely random. 9. Completely fixed.    

        type_sequence = 6; % Do not do FRA, analysis will not work. Specific function for it.
        recParams= recDataForSubject(subject, subject_rec, type_sequence);

        % Limits of the drift map
        triggersInterest = blockEvents{recParams.idxSequence(idx_chosen)};
        t_start = triggersInterest(1);
        t_end = triggersInterest(end);
        
        %% Get FRA info.
        
        min_amp = 0; % mV
        num_min_spikes = 50;

        % We find the clusters in FRA that have significant differences before and
        % after triggers with 70 dB.
        FRA = clustersInFRA(myKsDir, direc, sp, blockEvents, subject, subject_rec, sec_before_trigger);

        %% Finding clusters with Omission Response

        [soi, soi_trials, triggersInterest] = getDataFromSequence(sp, blockEvents, recParams, idx_chosen, sec_before_trigger);
        [cl_omission, pre_trig_omi, post_trig_omi, p_omi] = findSignificantClustersOmission(soi.spTimes, soi.spSites, clusters_with_spikes, triggersInterest, sec_before_trigger);

        if ~isempty(cl_omission)
            for z = 1:length(cl_omission)
                [amp, good] = findAmpAndSorting(myKsDir, cl_omission(z));
                fprintf('Cluster ID: %d, Amp: %.1f, Good: %d\n', cl_omission(z), amp, good);
            end
        end

        %% Select the cluster you want to see.
        
        n_bins = 50;
        
        if ~isempty(cl_omission)
            
            continue_yes = input('Do you want to visualize any cluster? Y/N: ', 's');

            cl_idx = -1;

            while strcmpi(continue_yes, 'Y')
                

                while cl_idx == -1
                    cl_idx = -1; %#ok<NASGU>
                    try
                        cl_idx = input('Which cluster is of your interest? ');
                    catch
                        error('Introduce an index from 0');
                    end
                end
                
                figure_width = 43.5; %cm
                figure_hight = 31.5; %cm
                figure('NumberTitle','off','name', 'figure_size_control', 'units', 'centimeters', 'color','w', 'position', [0, 0, figure_width, figure_hight], 'PaperSize', [figure_width, figure_hight]); % this is the trick!

                subplot(2,2,1)
                plotFRA_tuning(FRA, myKsDir, cl_idx, n_bins)
                
                for p=1:3
                    
                    % Plot the raster plot of 3 spaced sequences.
                    subplot(2,2,p+1)

                    idx_chosen = p;

                    [~, soi_trials, ~] = getDataFromSequence(sp, blockEvents, recParams, idx_chosen, sec_before_trigger);
                    rasterplot(soi_trials.spTimes, soi_trials.spSites, recParams, idx_chosen, cl_idx, type_sequence);
                    
                end
                sgtitle(subject_rec + " imec" +imec)
                
                figure_width = 43.5; %cm
                figure_hight = 31.5; %cm
                figure('NumberTitle','off','name', 'figure_size_control', 'units', 'centimeters', 'color','w', 'position', [0, 0, figure_width, figure_hight], 'PaperSize', [figure_width, figure_hight]); % this is the trick!
                
                hist4ProtocolsCluster(cl_idx, n_bins, sp, blockEvents, recParams, sec_before_trigger, type_sequence)
                sgtitle(subject_rec + " imec" +imec)
                
                cl_idx = -1;
                continue_yes = input('Do you want to visualize any other cluster? Y/N: ', 's');
            end
        end
    end
end

            %% Plot histogram if there is omission responses in that cluster
                
%                 figure
%                 
%                 cluster = cl_idx;
%                 ISI = recParams.ISITrial(idx_chosen);
%                 triggers = 0: ISI: (ISI*7+ISI/2); % we get 8 Triggers this way
% 
%                 spikes = [];
% 
%                 % See how many spikes between those Depths per Trial.
%                 for iTrial = 1:length(soi_trials.spTimes)
% 
%                     trial_sp = soi_trials.spTimes{iTrial}(soi_trials.spSites{iTrial}==cluster)';
%                     spikes = [spikes, trial_sp]; %#ok<AGROW>
%                 end
% 
%                 % Performing the histogram
%                 histogram(spikes, n_bins)
% 
%                 % Plot the trigger times in red vertical lines.
%                 for i=1:length(triggers)
%                     xline(triggers(i),'red')
%                 end
            
        %% Fixed spaced 4 of them raster plot with k different depths

%         % # bins in the histogram
%         n_bins = 60;
%         
%         figure_width = 43.5; %cm
%         figure_hight = 31.5; %cm
%         figure('NumberTitle','off','name', 'figure_size_control', 'units', 'centimeters', 'color','w', 'position', [0, 0, figure_width, figure_hight], 'PaperSize', [figure_width, figure_hight]); % this is the trick!
% 
%         histogramDiffDepths(Depth1, Depth2, n_bins, sp, blockEvents, recParams, sec_before_trigger, type_sequence)
%         sgtitle("M" + subject_rec + " imec" +imec)
% 
%         saveas(gcf,"S:\Fakultaet\MFZ\NWFZ\AGdeHoz\Martin\Figures\Histograms_divided_Depths\Hist_Depths_M"+subject_rec+"_imec"+imec+".fig")

        %% Fixed spaced: histogram of the first 4 between specified depths
% 
%         % % # bins in the histogram
%         % n_bins = 100;
%         
%         figure_width = 43.5; %cm
%         figure_hight = 31.5; %cm
%         figure('NumberTitle','off','name', 'figure_size_control', 'units', 'centimeters', 'color','w', 'position', [0, 0, figure_width, figure_hight], 'PaperSize', [figure_width, figure_hight]); % this is the trick!
% 
%         hist4Protocols(Depth1, Depth2, n_bins, sp, blockEvents, recParams, sec_before_trigger, type_sequence)
%         sgtitle("M" + subject_rec + " imec" +imec)
% 
%         saveas(gcf,"S:\Fakultaet\MFZ\NWFZ\AGdeHoz\Martin\Figures\Histograms_full_Depths\Hist_M"+subject_rec+"_imec"+imec+".fig")
%     

%% Do Histograms of trials
% 
% [soi, soi_trials, triggersInterest, trialsTriggers] = getDataFromSequence(sp, blockEvents, recParams, idx_chosen, sec_before_trigger);
% 
% figure
% 
% histSummedDepths(soi_trials.spTimes, soi_trials.spDepths, Depth1, Depth2, n_bins, recParams, idx_chosen, type_sequence, color)

%% Find best clusters from FRA, with high amplitude and # of spikes and get the tuning curve

% min_amp = 0; % mV
% num_min_spikes = 50;
% 
% % We find the clusters in FRA that have significant differences before and
% % after triggers with 70 dB.
% FRA = clustersInFRA(myKsDir, direc, sp, blockEvents, subject, sec_before_trigger);
% 
% % We select also from those only the ones with a minimum amplitude and
% % number of spikes.
% [cl_FRA_amplitude, FRA.cl_good, FRA.cl_amp] = selectGoodClusters(myKsDir, FRA.spSites, min_amp, num_min_spikes);
% 
% % We take the clusters that belong to both groups and sort them by amp.
% [~, pos1, ~] = intersect(cl_FRA_amplitude, FRA.clusters);
% pos1 = sort(pos1);
% 
% FRA.cl_final = cl_FRA_amplitude(pos1);
% FRA.cl_good_final = FRA.cl_good(pos1);
% FRA.cl_amp_final = FRA.cl_amp(pos1);

%% Plot the tuning curves of the clusters we want

% Number of bins in which we divide the frequencies into.
% n_bins = 50;
% 
% figure
% plotFRA_tuning(FRA, FRA.cl_final, n_bins)

%% Finding good clusters to read from (based on amplitude and # of spikes)
% 
% [cl_chosen, kiloGood] = selectGoodClusters(myKsDir, soi.spSites, min_amp, num_min_spikes);
% [cl_final, pre_trig, post_trig, h, p] = findSignificantCluster(soi.spTimes, soi.spSites, cl_chosen, triggersInterest, sec_before_trigger);

%% Plotting the spikes of that cluster

% Concatenating for representation
% chosen_spikes = [];
% trials_pertenence = [];
% amount_trials = 0;
% for i=1:length(soi_trials.spSites)
%     trial_sp = soi_trials.spTimes{i}(soi_trials.spSites{i}==cl_chosen)';
%     if ~isempty(trial_sp)
%         chosen_spikes = [chosen_spikes, trial_sp];
%         trials_pertenence = [trials_pertenence, i*ones(length(trial_sp),1)'];
%     end
% end

%% Loading data from kilosort/phy easily and loading the trigger events

% There are sometimes a directory name with an extra underline.
% try 
%     sp = loadKSdir(myKsDir{1});
%     myKsDir = myKsDir{1};
% catch
%     sp = loadKSdir(myKsDir{2});
%     myKsDir = myKsDir{2};
% end
    
% sp.st are spike times in seconds
% % sp.clu are cluster identities
% spikes from clusters labeled "noise" have already been omitted

%% Choosing cluster depending on the number of spikes in the SOI

% num_min_spikes = 50;

% % Take the IDs and whether a cluster is good or not
% cl_id = C{1,1}; 
% cl_good = contains(C2{1,2},'good');
% 
% % Open the file where th information of how good a cluster is is contained
% fid2 = fopen(myKsDir+'\cluster_KSLabel.tsv');
% C2 = textscan(fid, '%d %s', 'HeaderLines', 1);
% fclose(fid2);
% 
% % Histogram to see how many spikes a cluster has in the SOI.
% spikes_in_cl = hist(double(spikesSites), double(cl_id));
% [num_spikes, max_counts] = sort(spikes_in_cl,'descend');
% 
% num_good_cl = sum(num_spikes>num_min_spikes);
% cl_chosen = [];
% 
% % Taking only the number of cluster with minimum X spikes that are
% % considered by Kilosort to be good.
% for j = 1:num_good_cl
%     if cl_good(max_counts(j))== true  
%         cl_chosen = [cl_chosen, max_counts(j)-1]; %#ok<AGROW>
%     end
% end