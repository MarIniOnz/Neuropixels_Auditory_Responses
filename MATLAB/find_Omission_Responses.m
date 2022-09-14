% Only works directly with GAMZE's data. You will need to slightly change 
% some functions to adapt it to your data (you can check Poppy's adaptation
% for it).
% 
% 1) Get the spike times, clusters and triggers of a sequence, subject, and
% imec.
% 
% 2) Getting only the data from our Sequence of Interest (SOI), introducing
% see.
% the type of sequence, and which repetition of that sequence we want to
% 
% 3) Checking for any clusters with omission responses.
%
% 4) Does the histogram of spikes of that cluster of all Trials.

%% Restart

clc
clear

%% Get the directories and the recording parameters

subject_rec = 609;
imec = 0;

[direc, myKsDir, num_sequences] = filesForSubject(subject_rec, imec);

% In case we use b, c sequences.
if subject_rec>1000 && subject_rec<10000
    subject=floor(subject_rec/10);
elseif subject_rec > 10000
    subject=floor(subject_rec/100);
else 
    subject = subject_rec;
end

%% Getting the spike times and the clusters they belong to

[sp.spikeTimes, sp.spikeAmps, sp.spikeDepths, sp.spikeSites] = ksDriftmap(myKsDir,1);
clusters_with_spikes = unique(sp.spikeSites);

%% Extracting the trigger events of the sequence of our block of interest

blockEvents = separateTriggersIntoBlocks(direc.folder, num_sequences);
sec_before_trigger = 0.03;

if subject_rec == 61300 || subject_rec == 61400
    blockEvents(1:15)=[];
end

idx_chosen = 1; % Which time of the block we want the sequence from (1st time we have a spaced, 2nd time we present it...)

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

[soi, soi_trials, triggersInterest] = getDataFromSequence(sp, blockEvents, recParams, idx_chosen, sec_before_trigger);

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

%% Finding clusters with Omission Response

[cl_omission, pre_trig_omi, post_trig_omi, p_omi] = findSignificantClustersOmission(soi.spTimes, soi.spSites, clusters_with_spikes, triggersInterest, sec_before_trigger);

%% Plot histogram if there is omission responses in that cluster

n_bins = 60;
cluster = cl_omission(7);
ISI = recParams.ISITrial(idx_chosen);
triggers = 0: ISI: (ISI*7+ISI/2); % we get 8 Triggers this way

spikes = [];

figure

% See how many spikes between those Depths per Trial.
for iTrial = 1:length(soi_trials.spTimes)
    
    trial_sp = soi_trials.spTimes{iTrial}(soi_trials.spSites{iTrial}==cluster)';
    spikes = [spikes, trial_sp]; %#ok<AGROW>
end

% Performing the histogram
histogram(spikes, n_bins)

% Plot the trigger times in red vertical lines.
for i=1:length(triggers)
    xline(triggers(i),'red')
end
