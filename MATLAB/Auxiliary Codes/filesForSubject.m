function [direc, myKsDir, num_sequences] = filesForSubject(subject_rec, imec)
%
% Function that gets the directories for the recordings of a subject, once
% its data has been introduced in the Excel sheet subjects_characteristics.
%    
% Inputs:
%     subjec_rec: num, subject ID adding a 1 or 2 for Lateral or Dorsal
%     imec: num, probe number.
%     
% Outputs:        
%     direc: directory in which the files are located for that subject.
%     myKsDir: directory where the whole recording in contained for
%     Kylosort.
%     num_sequences: number of sequences in recording.

imec = num2str(imec);

% Read the Excel Sheet
[nums,info] = xlsread('\\charite.de\centren\Fakultaet\MFZ\NWFZ\AGdeHoz\Martin\Documents\subjects_characteristics.xlsx');
[~, colFolder] = find(contains(info,'Folder of interest'));
[~, colDate] = find(contains(info,'date'));
[~, colSequences] = find(contains(info,'Sequences in file'));
[row, ~] = find(nums==subject_rec);

% Extract the information from the Excel Sheet
date0 = info{row+1, colDate};
date = convertCharsToStrings(date0);
folder0 = info{row+1, colFolder};
folder = convertCharsToStrings(folder0);
num_sequences = nums(row, colSequences);

% Store the file path and date in the direc structure
direc.date = date;
direc.folder = '\\charite.de\\centren\\Fakultaet\\MFZ\\NWFZ\\AG-deHoz-Scratch\\Neuropixels\\'+date+'\\'+folder+'\\';

myKsDir = direc.folder + dir(direc.folder+'*imec'+ imec).name;

% Check whether we have a folder within that folder with the files from KS2
if exist(myKsDir + '\\imec' + imec +'_ks2_orig', 'dir')
    myKsDir = myKsDir + '\\imec' + imec +'_ks2_orig';
end