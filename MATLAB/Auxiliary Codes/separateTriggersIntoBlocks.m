function [blockEvents] = separateTriggersIntoBlocks(directory, numSequences)
%
% This function extracts the trigger times from a .adj file and separates
% in columns the blocks of triggers corresponding to different sequences.
%
% Inputs:
%     directory: str, which is the main directory where the files are stored.
%     numSequences: num, total number of sequences in the recording.
%
% Outputs:
%     blockEvents: double, trigger times in each sequence block, ordered
%     chronologically.

% Read the trigger times from that file.
file_loc = dir(directory+"*.adj.txt.");
if isempty(file_loc)
    file_loc = dir(directory+"*.XA_0_0.txt.");
end

% Reading the file that contains the triggers.
filename = file_loc.name;
eventTimes = dlmread(directory+filename);

% Variables initialization.
n = 1000;
timeDiffs = diff(eventTimes);
incrmTime = timeDiffs(1);
diffTime = incrmTime;

% Check how big the difference between triggers in time there is: between
% sequences, triggers were not used, so we can divide the triggers in
% different sequences once we know how many blocks there were
% (numSequences).
while n > numSequences-0.1
    locsChangeSequence= timeDiffs > diffTime;
    n = sum(locsChangeSequence);
    diffTime = diffTime + incrmTime;
end

% Adding a final locations, so we can use a +1 at the last loop.
locs = [0 transpose(find(locsChangeSequence)) length(eventTimes)]; 

% Locate the trigger times in different cells.
blockEvents = cell(n+1,1);
for i=1:n+1
    blockEvents{i,1} = eventTimes(locs(i)+1:locs(i+1));
    if mod(length(blockEvents{i,1}),2)~= 0 && i~=1
        blockEvents{i,1}(end) =[];
    end
end