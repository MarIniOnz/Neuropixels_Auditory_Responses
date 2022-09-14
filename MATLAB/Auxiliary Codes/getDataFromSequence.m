function [soi, soi_trials, triggersInterest] = getDataFromSequence(sp, blockEvents, recParams, idx_chosen, sec_before_trigger)
%
% Function that gets the spikes, the blocks of triggers, the recording parameters,
% and the index of the sequence we want. With those variables, it slices the spike times
% and triggers so that we only have our region of interest.
%    
% Inputs:
%     sp: struct, contains all the spikes' information from the whole
%     recording.
%     blockEvents: cell, contains all the triggers divided in their 
%     different sequences.
%     recParamas: struct, recording parameters.
%     idx_chosen: num, index of the sequence we want (1st, 2nd ... time     
%     that the sequence is performed in the recording).
%     sec_before_trigger: num, how many seconds before the triggers we want
%     to see to perform the t-test.
%
% Outputs:        
%     soi: struct, containing spikes of the Sequence of Interest (SOI).
%     Fields:
%         idx: double, index of the Sequence of Interest chosen (from the
%         number of total sequences in blockEvents.
%         spTimes: double, time of the spikes.
%         spSites: uint16, cluster to which the spikes belong.
%         spAmps: double, amplitude of the spikes.
%         spDepths: double, depth of the spikes.
%     soi_trials: struct, containing spikes of the Sequence of Interest 
%     divided into trials.
%     Fields:
%         trialsTriggers: double, trigger times arranged in trials.
%         spTimes: double, time of the spikes.
%         spSites: uint16, cluster to which the spikes belong.
%         spAmps: double, amplitude of the spikes.
%         spDepths: double, depth of the spikes.
%     triggersInterest: triggers appearing in the SOI.

% Defining the start and end of the SOI, and the triggers within it.
soi.idx = recParams.idxSequence(idx_chosen); % Sequence of Interest Chosen
triggersInterest = blockEvents{soi.idx};

t_start = triggersInterest(1)-sec_before_trigger; % check the activity prior to the start of the sequence.
t_end = triggersInterest(end)+ recParams.ISITrial(idx_chosen); % check the activity posterior to the end of the sequence.

%% Using only the spikes comprehended in the SOI

% Storing the spike information only pertaining to that SOI.
soi.spTimes = sp.spikeTimes(sp.spikeTimes>t_start & sp.spikeTimes<t_end);
soi.spAmps = sp.spikeAmps(sp.spikeTimes>t_start & sp.spikeTimes<t_end);
soi.spDepths = sp.spikeDepths(sp.spikeTimes>t_start & sp.spikeTimes<t_end);
soi.spSites = sp.spikeSites(sp.spikeTimes>t_start & sp.spikeTimes<t_end);

%% Separate the SOI in trials

if recParams.ISITrial(idx_chosen)== 0.1
    % Problem with 1st Trial in the case of 0.1 ISI, so we take it out in  
    % this case and correct for it (displacing everything by one, and thus,
    % taking out also the last trial.
    triggersInterest = triggersInterest(2:end-7);
    
    % Restructure the triggers into trials.
    soi_trials.trialsTriggers = reshape(triggersInterest,[8,length(triggersInterest)/8]);
    
    % Taking the start of all trials.
    trialsStart = [soi_trials.trialsTriggers(1,2:end), triggersInterest(end-6)]; % We do not look at the first trial because there is problems with it.
else
    soi_trials.trialsTriggers = reshape(triggersInterest,[8,length(triggersInterest)/8]);
    trialsStart = [soi_trials.trialsTriggers(1,2:end), t_end];
end

% We introduce a shift at the end and at the start so we can see a
% bit before the triggers how the activity is.
for i=1:length(trialsStart)-1
    soi_trials.spTimes{i} = soi.spTimes(soi.spTimes>(trialsStart(i)-sec_before_trigger) & soi.spTimes<(trialsStart(i+1)-sec_before_trigger))-trialsStart(i);
    soi_trials.spAmps{i} = soi.spAmps(soi.spTimes>(trialsStart(i)-sec_before_trigger) & soi.spTimes<(trialsStart(i+1)-sec_before_trigger));
    soi_trials.spDepths{i} = soi.spDepths(soi.spTimes>(trialsStart(i)-sec_before_trigger) & soi.spTimes<(trialsStart(i+1)-sec_before_trigger));
    soi_trials.spSites{i} = soi.spSites(soi.spTimes>(trialsStart(i)-sec_before_trigger) & soi.spTimes<(trialsStart(i+1)-sec_before_trigger));
end
