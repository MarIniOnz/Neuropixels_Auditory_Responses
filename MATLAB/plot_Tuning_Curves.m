% Only works directly with GAMZE's data. You will need to slightly change 
% some functions to adapt it to your data (you can check Poppy's adaptation
% for it).

% 1) Get the spike times, clusters and triggers of a sequence, subject, and
% imec.

% 2) Finding the best clusters in the FRA and plotting their frequency
% tunings.

% 3) Choosing a cluster and plotting the frequency responses to the Random
% and Random fifty blocks of that cluster.

%% Restart

clc
clear

%% Get the directories and the recording parameters

subject_rec = 602;
imec = 0;

[direc, myKsDir, num_sequences] = filesForSubject(subject_rec, imec);

% In case we use b, c sequences.
if subject_rec>1000 && subject_rec<10000
    subject=floor(subjects(k)/10);
elseif subject_rec > 10000
    subject=floor(subjects(k)/100);
else 
    subject = subject_rec;
end


%% Getting the spike times and the clusters they belong to

[sp.spikeTimes, sp.spikeAmps, sp.spikeDepths, sp.spikeSites] = ksDriftmap(myKsDir,1);
clusters_with_spikes = unique(sp.spikeSites);

%% Extracting the trigger events of the sequence of our block of interest

blockEvents = separateTriggersIntoBlocks(direc.folder, num_sequences);
sec_before_trigger = 0.03;

%% Find best clusters from FRA, with high amplitude and # of spikes and get the tuning curve

min_amp = 100; % mV
num_min_spikes = 100;

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

%% Plot number of spikes versus the frequencies (normalized)

cl_id = FRA.cl_final(12);
median_cl = medians(12);

plotSpikesFrequencies(subject, subject_rec, imec, sp, blockEvents, direc, sec_before_trigger, FRA.cl_final, medians)
