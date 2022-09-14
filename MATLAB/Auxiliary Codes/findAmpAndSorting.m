function [amp, sorting]= findAmpAndSorting(myKsDir, cluster)
%
% Function that takes the ID of a cluster and returns its amplitude and
% sorting from KiloSort.
%
% Inputs:
%     myKsDir: str, where the KiloSort files are located.
%     cluster: num, ID of the cluster.
%     
% Outputs:
%     amp: double, amplitude of the cluster.
%     sorting: logical, whether the cluster is classified as good or not.


% Open the file where th information of how good a cluster is is contained
fid = fopen(myKsDir+"\cluster_Amplitude.tsv");
C = textscan(fid, '%d %f', 'HeaderLines', 1);
fclose(fid);

% Take the IDs and whether a cluster is good or not
cl_amp = C{1,2};

% Open the file where th information of how good a cluster is is contained
fid2 = fopen(myKsDir+"\cluster_KSLabel.tsv");
C2 = textscan(fid, '%d %s', 'HeaderLines', 1);
fclose(fid2);

% Take the IDs and whether a cluster is good or not
cl_good = contains(C2{1,2},'good');

% Take the data.
amp = cl_amp(cluster+1);
sorting = cl_good(cluster+1);
