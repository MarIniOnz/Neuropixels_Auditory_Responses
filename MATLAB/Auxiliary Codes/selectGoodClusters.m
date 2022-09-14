function [cl_chosen, kiloGood, kiloAmp] = selectGoodClusters(myKsDir, spSites, min_amp, num_min_spikes)
%
% Function that find the top N clusters with minimum X spikes in the
% sequence of interest (SOI) and with a minimum amplitude Y (mV) on that
% cluster. They are ordered by amplitude.
% 
% Inputs:
%     myKsDir: str, where is the corresponding ks folder located (given by 
%     getTheFilesForSubject.m).
%     spSites: uint16, clusters to which the spikes pertain.
%     min_amp: num, minimum amplitude of the cluster (mV) to be considered.
%     num_min_spikes: num, minimum number of spikes in SOI to be
%     considered.
%
% Outputs: 
%     cl_chosen: double, number id's of the clusters chosen with minimum
%     "min_amp" ampltiude in the SOI. It can contain ID num = 0 (py 
%     originated)
%     kiloGood: bool, vector containing whether the clusters are classified
%     by kiloSort as good or not.
%     kiloAmp: amplitude of the cluster given by kilosort (mV)

% Open the file where th information of how good a cluster is is contained
fid = fopen(myKsDir+"\cluster_Amplitude.tsv");
C = textscan(fid, '%d %f', 'HeaderLines', 1);
fclose(fid);

% Take the IDs and whether a cluster is good or not
cl_id = C{1,1};
cl_amp = C{1,2};

% Open the file where th information of how good a cluster is is contained
fid2 = fopen(myKsDir+"\cluster_KSLabel.tsv");
C2 = textscan(fid, '%d %s', 'HeaderLines', 1);
fclose(fid2);

% Take the IDs and whether a cluster is good or not
cl_good = contains(C2{1,2},'good');

% Histogram to see how many spikes a cluster has in the SOI.
spikes_in_cl = hist(double(spSites), double(cl_id));
[num_spikes, max_counts] = sort(spikes_in_cl,'descend');

% Sorting which clusters have the biggest amplitude.
[amp, cl_id_amps] = sort(cl_amp, 'descend');
num_good_cl = sum(amp > min_amp);
cl_id_amps = cl_id_amps(1:num_good_cl); % The clusters start at 0.

% From the clusters with biggest amplitude, take the ones that have a
% minimum number of spikes. Store in a vector whether they are "good".
cl_chosen =[];
kiloGood = [];
kiloAmp = [];

% Only taking the clusters with a minimum number of spikes AND a minimum
% amplitude, ordering them by amplitude.
for j = 1:num_good_cl
    idx = max_counts==cl_id_amps(j);
    if num_spikes(idx)>num_min_spikes 
        cl_chosen = [cl_chosen, cl_id_amps(j)-1]; %#ok<AGROW>
        kiloGood = [kiloGood, cl_good(cl_id_amps(j))]; %#ok<AGROW>
        kiloAmp = [kiloAmp, amp(j)]; %#ok<AGROW>
    end
end