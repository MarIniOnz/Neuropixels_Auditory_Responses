function FRA = clustersInFRA(myKsDir, direc, sp, blockEvents, subject, subject_rec, sec_before_trigger)
%
% Function to get the cluster with highest significant difference between
% spikes appearing right before and after the trigger in the FRA sequence.
%
% Inputs:
%     myKsDir: str, where is the corresponding ks folder located (given by
%     getTheFilesForSubject.m).
%     direc: struct, directory where the files are stored, you can obtain 
%     it from filesForSubject function.
%     blockEvents: cell, where triggers are stored.
%     subject: num, subject ID.
%     subjec_rec: num, subject ID adding a 1 or 2 for Lateral or Dorsal
%     recordings.
%     sec_before_trig: num, seconds to look and integrate before and after
%     trigger.
%
% Outputs:
%     FRA: struct, contains the trigger times and spike fields from the FRA
%     recording sequence.
%     Fields:
%         triggers: double, time of the triggers
%         spTimes: double, time of the spikes
%         spSites: uint16, cluster to which the spikes belong
%         clusters: int32, clusters that have significant difference between spikes
%         appearing right before and after the trigger in the FRA sequence
%         cl_good: str, kilosort sorting of those cluster (good, mua, bad(?))
%         p_FRA: double, p values of the difference between spikes before 
%         and after the trigger
%         data: struct, contains data of how the recording was performed 
%         (which triggers were activated, their frequencies and dBs)
%         desiredTrigs: double, indexes of the triggers that are 70 dB
%         freqTrig: double, frequencies of those triggers
%         pre_trig: double, number of spikes just before the trigger
%         post_trig: double, number of spikes just after the trigger

% Open the file where th information of how good a cluster is is contained
fid = fopen(myKsDir+"\cluster_KSLabel.tsv");
C = textscan(fid, '%d %s', 'HeaderLines', 1);
fclose(fid);

% Take the IDs and whether a cluster is good or mua
cl_id = C{1,1};
cl_good = C{1,2};

FRA.triggers = blockEvents{1};

% Different for that sequence.
if subject_rec == 6031 || subject_rec == 61340 || subject_rec == 61440
    FRA.triggers = blockEvents{end};
end

t_start = FRA.triggers(1) - sec_before_trigger; % check the activity prior to the start of the sequence
t_end = FRA.triggers(end) + sec_before_trigger; % check the activity posterior to the end of the sequence

%% Using only the spikes comprehended in the FRA

FRA.spTimes = sp.spikeTimes(sp.spikeTimes>t_start & sp.spikeTimes<t_end);
FRA.spSites = sp.spikeSites(sp.spikeTimes>t_start & sp.spikeTimes<t_end);

% Loading the parameters used to create the triggers.
recDataFRA = recDataForSubject(subject, subject_rec, 1);

% If it does not work automatically, introduce it by hand.
FRA.data = load("\\charite.de\centren\Fakultaet\MFZ\NWFZ\AG-deHoz-Scratch\Neuropixels\" + direc.date +"\" +recDataFRA.folderMat{1} + "\data.mat");

%% Finding which are the triggers that are 70 dB so we use only those.

trigs = FRA.data.Trseq;
FRA.desiredTrigs = [];
FRA.freqTrig = [];

for i= 1:length(trigs)

    x=num2str(trigs(i));
    out=contains((x(1)),'2');
    if out && length(x)>2 % The trigger was 70 dB 
        FRA.desiredTrigs = [FRA.desiredTrigs, i];
        switch length(x)
            case 2
                freq = FRA.data.freqlist(str2double(x(2)));
            case 3
                freq = FRA.data.freqlist(str2double(x(2:3)));
        end
        FRA.freqTrig = [FRA.freqTrig, freq];
    end
end

%% Finding which clusters have significant differences between the number of spikes they have before and after a trigger.

% Creating variables
FRA.pre_trig = zeros(length(FRA.desiredTrigs), length(cl_id));
FRA.post_trig = zeros(length(FRA.desiredTrigs), length(cl_id));
h = zeros(length(cl_id),1);
p = zeros(length(cl_id),1);

% Performing the sums and t-tests
for i = 1: length(cl_id)
    cl_spikes = FRA.spTimes(FRA.spSites==cl_id(i))';
    
    for k = 1:length(FRA.desiredTrigs)
        FRA.pre_trig(k,i) = sum(cl_spikes<FRA.triggers(FRA.desiredTrigs(k)) & cl_spikes>(FRA.triggers(FRA.desiredTrigs(k))-sec_before_trigger));        
        FRA.post_trig(k,i) = sum(cl_spikes>FRA.triggers(FRA.desiredTrigs(k)) & cl_spikes < (FRA.triggers(FRA.desiredTrigs(k))+sec_before_trigger));
    end
    
    [h(i),p(i)] = ttest2(FRA.pre_trig(:,i),FRA.post_trig(:,i));

end

% Storing only the clusters with p values below 0.05 in the t-test.
FRA.clusters = cl_id(h==1);
FRA.p_FRA = p(h==1);
FRA.cl_good = cl_good(cl_id(h==1)+1);