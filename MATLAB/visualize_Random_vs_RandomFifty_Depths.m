% Script to 

% 1) Get the spike times, clusters and triggers of a sequence, subject, and
% imec.

% 2) Finding the best clusters in the FRA and plotting their frequency
% tunings.

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

%% FOR loop

idx_chosen = 1; % Which time of the block we want the sequence from (1st time we have a spaced, 2nd time we present it...)

% 609 does not work for 2nd idx_chosen (not recorded for random fifty), 608
% either

for k = 3
    %% Get the directories and the recording parameters

    subject_rec = subjects(k);
    
    if subject_rec ==608
        continue
    end 
    
    for m = 1:2
        imec = m-1;
        
        if imec==0
            Depth1 = Depths_imec0(k,1);
            if Depth1 == 0
                continue
            end
            Depth2 = Depths_imec0(k,2);
        else
            Depth1 = Depths_imec1(k,1);
            if Depth1 == 0
                continue
            end
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
        sec_before_trigger = 0.08; % Default: 0.03 s.

        min_amp = 30; % mV
        num_min_spikes = 100;

        if subject_rec == 61340 || subject_rec == 61440
            blockEvents(1:15)=[];
        end

       %% Select the cluster you want to see.
        figure_height = 43.5; %cm
        figure_width = 31.5; %cm
        figure('NumberTitle','off','name', 'figure_size_control', 'units', 'centimeters', 'color','w', 'position', [0, 0, figure_width, figure_height], 'PaperSize', [figure_width, figure_height]); % this is the trick!

%         Sequence = plotSpikesFrequenciesDepths(subject, subject_rec, direc, imec, sp, blockEvents, sec_before_trigger, Depth1, Depth2, idx_chosen);
%         hold on 
%         subplot(5,1,1)
%         f = fit(log10(Sequence(2,2).frequencies)',Sequence(2,2).post_trig','gauss2');
%         plot(f);
        Sequence = plotSpikesFrequenciesDepthsFit(subject, subject_rec, direc, imec, sp, blockEvents, sec_before_trigger, Depth1, Depth2, idx_chosen);
        pause

%         saveas(gcf,"Z:\Fakultaet\MFZ\NWFZ\AGdeHoz\Martin\Figures\Random_vs_randomFifty\Depths_"+subject_rec+"_imec"+imec+".fig")
        
    end
end