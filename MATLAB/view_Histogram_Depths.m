% Only works directly with GAMZE's data. You will need to slightly change 
% some functions to adapt it to your data (you can check Poppy's adaptation
% for it).

% Script to 

% 1) Get the spike times, clusters and triggers of a sequence, subject, and
% imec.

% 2) Plotting a drift map of the desired sequence (if you uncoment it).

% 3) Perform histograms with the number of spikes that are contained within
% the speficied depths.

% 3) Choosing a type of sequence (like Spaced) and plotting the raster
% plots of the spikes per trial in the best, most significant clusters.

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
subjects = [607, 608, 609, 611, 6091, 6092, 6111, 602, 603, 6021, 6031, 60900, 60930, 612, 613, 61340, 614, 61440, 612]; % Took out 608, problem with recParams (it does not read the idx_chosen well)

Depths_imec0 = [1400, 1900; 1500, 1900; 2500, 3500; 2100, 2500; 2000, 3000; 2000, 2800; 2100, 2400; 2300, 2700; 1800, 2500; 800, 1500; 1250, 2300; 0,0; 0,0; 0,0; 0,0; 0,0; 0, 0; 0, 0;0,0];
Depths_imec1 = [2500, 3200; 2200, 3200; 2000, 3000; 2100, 2800; 2000, 2800; 1500, 2000; 2400, 2570; 0, 0; 0, 0; 0, 0; 0, 0; 1, 1; 2,2; 1, 1; 1900, 3000; 1800, 2800; 2850, 3250; 2850, 3250; 1400,2800];

%% Get the directories and the recording parameters

k = 3; % Select the subject. 
subject_rec = subjects(k);
imec = 1;

if imec==0
    Depth1 = Depths_imec0(k,1);
    Depth2 = Depths_imec0(k,2);
else
    Depth1 = Depths_imec1(k,1);
    Depth2 = Depths_imec1(k,2);
end

if subject_rec==60900
    direc.folder = "\\charite.de\\centren\\Fakultaet\\MFZ\\NWFZ\\AG-deHoz-Scratch\\Neuropixels\\15_10_2020\\GG_M609__g0_1_high_processed\\";
    direc.date = "15_10_2020";
    num_sequences = 18;
    myKsDir = "\\charite.de\\centren\\Fakultaet\\MFZ\\NWFZ\\AG-deHoz-Scratch\\Neuropixels\\15_10_2020\\GG_M609__g0_1_high_processed\\GG_M609__g0_imec" +imec+"\\imec"+imec+"_ks2_orig";      
elseif subject_rec==60930
    direc.folder = "\\charite.de\\centren\\Fakultaet\\MFZ\\NWFZ\\AG-deHoz-Scratch\\Neuropixels\\15_10_2020\\GG_M609__g0_1_wide_processed\\";
    direc.date = "15_10_2020";
    num_sequences = 18;
    myKsDir = "\\charite.de\\centren\\Fakultaet\\MFZ\\NWFZ\\AG-deHoz-Scratch\\Neuropixels\\15_10_2020\\GG_M609__g0_1_wide_processed\\GG_M609__g0_imec" +imec+"\\imec"+imec+"_ks2_orig";
elseif subject_rec==61340
    direc.folder = "\\charite.de\\centren\\Fakultaet\\MFZ\\NWFZ\\AG-deHoz-Scratch\\Neuropixels\\15_12_2020\\GG_M613_g0_t21,49\\catgt_GG_M613_g0\\";
    direc.date = "15_12_2020";
    num_sequences = 21;
    % Sequences start at T24 and end at T49,

    myKsDir = "\\charite.de\\centren\\Fakultaet\\MFZ\\NWFZ\\AG-deHoz-Scratch\\Neuropixels\\15_12_2020\\GG_M613_g0_t21,49\\catgt_GG_M613_g0\\GG_M613_g0_imec" +imec+"\\imec"+imec+"_ks2_orig";
elseif subject_rec==61440
    direc.folder = "\\charite.de\\centren\\Fakultaet\\MFZ\\NWFZ\\AG-deHoz-Scratch\\Neuropixels\\16_12_2020\\GG_M614_g0_t21,49\\catgt_GG_M614_g0\\";
    direc.date = "16_12_2020";
    num_sequences = 21;
    % Sequences start at T24 and end at T49,
    myKsDir = "\\charite.de\\centren\\Fakultaet\\MFZ\\NWFZ\\AG-deHoz-Scratch\\Neuropixels\\16_12_2020\\GG_M614_g0_t21,49\\catgt_GG_M614_g0\\GG_M614_g0_imec" +imec+"\\imec"+imec+"_ks2_orig";

else
    [direc, myKsDir, num_sequences] = filesForSubject(subject_rec, imec);
end

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

if subject_rec == 61340 || subject_rec == 61440
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

%% Ploting a drift map if you want

%         figure; plotDriftmap_livia(spikeTimes, spikeAmps, spikeDepths);
%         y=[min(spikeDepths) max(spikeDepths)];

figure;

plotDriftmap_livia(sp.spikeTimes(sp.spikeTimes>t_start & sp.spikeTimes<t_end), sp.spikeAmps(sp.spikeTimes>t_start & sp.spikeTimes<t_end), sp.spikeDepths(sp.spikeTimes>t_start & sp.spikeTimes<t_end));

%         for i=2:2:length(triggersInterest)
%             xline(triggersInterest(i),'red')
%         end
for i=1:1:length(triggersInterest)
    xline(triggersInterest(i),'red')
end

xlabel('Time (s)')
ylabel('Probe Position (um)')
%         saveas(gcf,'X:\Fakultaet\MFZ\NWFZ\AG-deHoz-Scratch\Martin\Figures_Spaced_2305\M'+subject+'_imec'+imec+'.fig')

%% Fixed spaced 4 of them raster plot with k different depths

% # bins in the histogram
n_bins = 60;
num_sequences = 3;

figure_width = 43.5; %cm
figure_hight = 31.5; %cm
figure('NumberTitle','off','name', 'figure_size_control', 'units', 'centimeters', 'color','w', 'position', [0, 0, figure_width, figure_hight], 'PaperSize', [figure_width, figure_hight]); % this is the trick!

histogramDiffDepths(Depth1, Depth2, n_bins, sp, blockEvents, recParams, sec_before_trigger, type_sequence, num_sequences)
sgtitle("M" + subject_rec + " imec" +imec)

% saveas(gcf,"S:\Fakultaet\MFZ\NWFZ\AGdeHoz\Martin\Figures\Histograms_divided_Depths\Hist_Depths_M"+subject_rec+"_imec"+imec+".fig")

%% Fixed spaced: histogram of the first 4 between specified depths

% Among which lengths we want to look at.
% Depth1 = 1400;
% Depth2 = 2000;

% # bins in the histogram
n_bins = 60;

figure_width = 43.5; %cm
figure_hight = 31.5; %cm
figure('NumberTitle','off','name', 'figure_size_control', 'units', 'centimeters', 'color','w', 'position', [0, 0, figure_width, figure_hight], 'PaperSize', [figure_width, figure_hight]); % this is the trick!

hist4ProtocolsDepths(Depth1, Depth2, n_bins, sp, blockEvents, recParams, sec_before_trigger, type_sequence)
sgtitle("M" + subject_rec + " imec" +imec)

% saveas(gcf,"S:\Fakultaet\MFZ\NWFZ\AGdeHoz\Martin\Figures\Histograms_full_Depths\Hist_M"+subject_rec+"_imec"+imec+".fig")