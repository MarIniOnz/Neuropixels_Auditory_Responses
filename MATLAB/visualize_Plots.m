% Only works directly with GAMZE's data. You will need to slightly change 
% some functions to adapt it to your data (you can check Poppy's adaptation
% for it).

% 1) Get the spike times, clusters and triggers of a sequence, subject, and
% imec.

% 2) Finding the best clusters in the FRA and plotting their frequency
% tunings.

% 3) Choosing a cluster and plot the raster plots of the spikes per trial
% in the 1st, 2nd, and 3rd time we get a spaced sequence, along with the
% frequency tuning extracted from the FRA.

%% Restart

clc
clear

%% add the repositories to your path

addpath(genpath('\\charite.de\centren\Fakultaet\MFZ\NWFZ\AG-deHoz-Scratch\Kilosort3')) % path to kilosort folder
addpath('\\charite.de\centren\Fakultaet\MFZ\NWFZ\AG-deHoz-Scratch\Kilosort3\npy-matlab-master\npy-matlab') % for converting to Phy
addpath('\\charite.de\centren\Fakultaet\MFZ\NWFZ\AG-deHoz-Scratch\Kilosort3\spikes-master\visualization')
addpath(genpath('\\charite.de\centren\Fakultaet\MFZ\NWFZ\AGdeHoz\Martin'))

%% Get the directories and the recording parameters

subject_rec = 609;
imec = 1;

[direc, myKsDir, num_sequences] = filesForSubject(subject_rec, imec);

% In case we use b, c sequences.
if subject_rec>1000 
    subject=floor(subject/10);
else
    subject= subject_rec;
end

%% Getting the spike times and the clusters they belong to

[sp.spikeTimes, sp.spikeAmps, sp.spikeDepths, sp.spikeSites] = ksDriftmap(myKsDir,1);
clusters_with_spikes = unique(sp.spikeSites);

%% Extracting the trigger events of the sequence of our block of interest

blockEvents = separateTriggersIntoBlocks(direc.folder, num_sequences);
sec_before_trigger = 0.03;

%% Find best clusters from FRA, with high amplitude and # of spikes and get the tuning curve

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

%% Plot the tuning curves of the clusters we want.

% Number of bins in which we divide the frequencies into.
n_bins = 40;

figure
medians = plotFRA_tuning(FRA, myKsDir, FRA.cl_final, n_bins);

%% Selecting the type of sequence we want and the recording parameters for that type of sequence

% Select which kind of sequence you want to look at: Legend
    % 1. FRA. 2. Fixed random. 3. Random fifty. 4. Random fifty silence.
    % 5. Fixed fifty silence. 6. Fixed spaced. 7. Fixed compressed.
    % 8. Completely random. 9. Completely fixed.     

type_sequence = 6; % Do not do FRA, analysis will not work. Specific function for it.
recParams= recDataForSubject(subject, subject_rec, type_sequence);
idx_chosen=1;

[soi, ~, triggersInterest] = getDataFromSequence(sp, blockEvents, recParams, idx_chosen, sec_before_trigger);

% Find clusters significantly different before and after the trigger.
[clusters_significant, pre_trig, post_trig, p_values] = findSignificantCluster(soi.spTimes, soi.spSites, clusters_with_spikes, triggersInterest, sec_before_trigger);

% Find clusters with minimum number of spikes and amplitude.
[clusters_good, cl_good, cl_amp] = selectGoodClusters(myKsDir, soi.spSites, min_amp, num_min_spikes);

% We take the clusters that belong to both groups and sort them by
% AMPLITUDE.
[~, pos1, ~] = intersect(clusters_good, clusters_significant);
pos1_sort = sort(pos1);

cl_final = clusters_good(pos1_sort);
cl_good_final = clusters_good(pos1_sort);
cl_amp_final = clusters_good(pos1_sort);

%% Open phy

system("start cmd.exe.")
fprintf('pushd \\' + myKsDir + '\n')
fprintf("phy template-gui params.py" + "\n")

%% Select the cluster you want to see.

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
    
    figure
    
    subplot(2,2,1)
    median_point = plotFRA_tuning(FRA, myKsDir, cl_idx, n_bins);
    
    subplot(2,2,2)
    
%     type_sequence = 6; 
%     recParams= recDataForSubject(subject, type_sequence);
    idx_chosen = 1;
    
    [~, soi_trials, ~] = getDataFromSequence(sp, blockEvents, recParams, idx_chosen, sec_before_trigger);
    rasterplot(soi_trials.spTimes, soi_trials.spSites, recParams, idx_chosen, cl_idx, type_sequence);
    
    subplot(2,2,3)
    
%     type_sequence = 6; 
%     recParams= recDataForSubject(subject, type_sequence);
    idx_chosen = 2;
    
    [~, soi_trials, ~] = getDataFromSequence(sp, blockEvents, recParams, idx_chosen, sec_before_trigger);
    rasterplot(soi_trials.spTimes, soi_trials.spSites, recParams, idx_chosen, cl_idx, type_sequence);
    
    subplot(2,2,4)
    
%     type_sequence = 6; 
%     recParams= recDataForSubject(subject, type_sequence);
    idx_chosen = 3;
    
    [~, soi_trials, ~] = getDataFromSequence(sp, blockEvents, recParams, idx_chosen, sec_before_trigger);
    rasterplot(soi_trials.spTimes, soi_trials.spSites, recParams, idx_chosen, cl_idx, type_sequence);
    
    cl_idx = -1;
    continue_yes = input('Do you want to visualize any other cluster? Y/N: ', 's');
    
end    