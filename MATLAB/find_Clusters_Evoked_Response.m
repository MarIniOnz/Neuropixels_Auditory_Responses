% Only works directly with GAMZE's data. You will need to slightly change 
% some functions to adapt it to your data (you can check Poppy's adaptation
% for it).

% 1) Get the spike times, clusters and triggers of a sequence, subject, and
% imec.

% 2) Finding the best clusters in the FRA and plotting their frequency
% tunings.

% 3) Choosing a cluster and plotting the frequency tuning from the FRA and
% the raster plot of the first 3 recordings of the specified type of
% sequence (like spaced). It also plots the first 4 recordings of the
% specified type of sequence in histogram shape.

%% Restart

clc
clear

%% Add the repositories to your path

addpath(genpath('\\charite.de\centren\Fakultaet\MFZ\NWFZ\AG-deHoz-Scratch\Kilosort3')) % path to kilosort folder
addpath('\\charite.de\centren\Fakultaet\MFZ\NWFZ\AG-deHoz-Scratch\Kilosort3\npy-matlab-master\npy-matlab') % for converting to Phy
addpath('\\charite.de\centren\Fakultaet\MFZ\NWFZ\AG-deHoz-Scratch\Kilosort3\spikes-master\visualization')
addpath(genpath('\\charite.de\centren\Fakultaet\MFZ\NWFZ\AGdeHoz\Martin'))

%% Subjects defined and so on

% Being 60900, 60930: 609 High and Wide respectively. 61340, 61440, second block of recordings.
subjects = [607, 608, 609, 611, 6091, 6092, 6111, 602, 603, 6021, 6031];

Depths_imec0 = [1400, 1900; 1500, 1900; 2000, 3500; 2100, 2500; 2000, 3000; 2000, 2800; 2100, 2400; 2300, 2800; 1300, 2500; 800, 1500; 1250, 2300];
Depths_imec1 = [2500, 3200; 2200, 3200; 2000, 3000; 2100, 2800; 2000, 2800; 1500, 2000; 2400, 2570; 0, 0; 0, 0; 0, 0; 0, 0];

%% FOR loop

for k = 8:9
    %% Get the directories and the recording parameters

    subject_rec = subjects(k);
    
    for m = 2
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
        recParams = recDataForSubject(subject, subject_rec, type_sequence);

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
            
        % We select also from those only the ones with a minimum amplitude and
        % number of spikes.
        [cl_FRA_amplitude, FRA.cl_good, FRA.cl_amp] = selectGoodClusters(myKsDir, FRA.spSites, min_amp, num_min_spikes);

        % We take the clusters that belong to both groups and sort them by amp.
        [~, pos1, ~] = intersect(cl_FRA_amplitude, FRA.clusters);
        pos1_sort = sort(pos1);

        FRA.cl_final = cl_FRA_amplitude(pos1_sort);
        FRA.cl_good_final = FRA.cl_good(pos1_sort);
        FRA.cl_amp_final = FRA.cl_amp(pos1_sort);

        cl_chosen = FRA.cl_final;
        
        %% Plot the tuning curves of the clusters we want.

        % Number of bins in which we divide the frequencies into.
        n_bins = 40;

        figure
        medians = plotFRA_tuning(FRA, myKsDir, FRA.cl_final, n_bins);

        %% Open phy if you want to look more in detail at the clusters to decide which want to look at.

        system("start cmd.exe.")
        fprintf('pushd \\' + myKsDir + '\n')
        fprintf("phy template-gui params.py" + "\n")

        %% Select the cluster you want to see.
        
        n_bins = 50;
        
        continue_yes = input('Do you want to visualize any cluster? Y/N: ', 's');
        cl_idx = -1;
        
        % While loop.
        while strcmpi(continue_yes, 'Y')

            while cl_idx == -1
                cl_idx = -1; %#ok<NASGU>
                try
                    cl_idx = input('Which cluster is of your interest? ');
                catch
                    error('Introduce an index from 0');
                end
            end
            
            % Set a big enough figure.
            figure_width = 43.5; %cm
            figure_hight = 31.5; %cm
            figure('NumberTitle','off','name', 'figure_size_control', 'units', 'centimeters', 'color','w', 'position', [0, 0, figure_width, figure_hight], 'PaperSize', [figure_width, figure_hight]); % this is the trick!
            
            % Plot FRA tuning.
            subplot(2,2,1)
            medians = plotFRA_tuning(FRA, myKsDir, cl_idx, n_bins);

            for p=1:3

                % Plot the raster plot of 3 spaced sequences.
                subplot(2,2,p+1)

                idx_chosen = p;
                
                % Plot raster plot of 3 first repetitions of the type of
                % frequency.
                [~, soi_trials, ~] = getDataFromSequence(sp, blockEvents, recParams, idx_chosen, sec_before_trigger);
                rasterplot(soi_trials.spTimes, soi_trials.spSites, recParams, idx_chosen, cl_idx, type_sequence);

            end
            sgtitle(subject_rec + " imec" +imec)
            
            % Set a big enough figure.
            figure_width = 43.5; %cm
            figure_hight = 31.5; %cm
            figure('NumberTitle','off','name', 'figure_size_control', 'units', 'centimeters', 'color','w', 'position', [0, 0, figure_width, figure_hight], 'PaperSize', [figure_width, figure_hight]); % this is the trick!
            
            % Plots 4 histograms from the first 4 blocks of a type of
            % sequence (for example, first 4 spaced recordings), summing up 
            % the number of spikes in all trials of that cluster.
            hist4ProtocolsCluster(cl_idx, n_bins, sp, blockEvents, recParams, sec_before_trigger, type_sequence) % Might be an error because there are no 4 sequences of that type. Change code to reduce the number of sequneces.
            
            sgtitle(subject_rec + " imec" +imec)

            cl_idx = -1;
            continue_yes = input('Do you want to visualize any other cluster? Y/N: ', 's');
        end
    end
end