function  recParams = recDataForSubject(subject, subject_rec, type_sequence)
%
% Function that get the recording parameters under which the data for a
% specfic type of sequence was extracted. It has A LOT OF EXCEPTIONS at the
% end, introduced by hand, so you might need to introduce more exceptions
% for your data, adequating the data in the database with the data on your
% written files. 
%
% This only works for GAMZE's data (extracted from that database).
%    
% Inputs:
%     subject: num, subject ID.
%     subjec_rec: num, subject ID adding a 1 or 2 for Lateral or Dorsal
%     recordings.
%     type_sequence: num, which type of sequence you want to analyze (fixed
%     space, compressed, random fifty ...)
%         1. FRA. 2. Fixed random. 3. Random fifty. 4. Random fifty silence.
%         5. Fixed fifty silence. 6. Fixed spaced. 7. Fixed compressed.
%         8. Completely random. 9. Completely fixed.
%     
% Outputs:
%     recData: struct which contains parameters defining how the data was obtained
%         Fields:
%         idxSpaced: indexes of the recordings which were Spaced.
%         ISItrial: Interspace interval between triggers of the sequences.
%         freqType: Type of frequencies used in the recording.
%             1 old freqs, 2 low freqs, 3 high freqs, 4 wide freqs for frequency type. 
%             21 new freqs, low long tail
%             22 new freqs, high long tail
%             11 old, short tail
%             12 mixed 
%         playOrder: In which order they were played in their block
%         folderMat: str, folder in which the data.mat are stored
%         containing the information of that specific recording, like
%         which frequencies were played. (For example: "15_32").

% Load the recording database and get the recording data from it.
subject_str = num2str(subject);
load('\\charite.de\centren\Fakultaet\MFZ\NWFZ\AGdeHoz\Martin\Codes\Matlab\GamzeCodes\database_sept20_gamze.mat') %#ok<LOAD>
idxData = find(strcmp({database.mouse}, subject_str)==1);
recParams = struct;

% Depending on the type of sequence, we get different data from database.
switch type_sequence
        case 1
            seq_idx = database(idxData).FRA;
            recParams.ISITrialAll = database(idxData).FRA_ISI/1000; % in seconds
            recParams.folderMatAll = database(idxData).FRA_mat;
            recParams.playOrderAll = database(idxData).FRA_playingorder;
        case 2
            seq_idx = database(idxData).Seq_fixed_random;
            recParams.ISITrialAll = database(idxData).Seq_fixed_random_ISI/1000;
            recParams.freqTypeAll = database(idxData).Seq_fixed_random_freqtype;
            recParams.playOrderAll = database(idxData).Seq_fixed_random_playingorder;
            recParams.folderMatAll = database(idxData).Seq_fixed_random_mat;
        case 3
            seq_idx = database(idxData).Seq_random_fifty;
            recParams.ISITrialAll = database(idxData).Seq_random_fifty_ISI/1000;
            recParams.freqTypeAll = database(idxData).Seq_random_fifty_freqtype;
            recParams.playOrderAll = database(idxData).Seq_random_fifty_playingorder;
            recParams.folderMatAll = database(idxData).Seq_random_fifty_mat;
        case 4
            seq_idx = database(idxData).Seq_random_fiftysilence;
            recParams.ISITrialAll = database(idxData).Seq_random_fiftysilence_ISI/1000;
            recParams.freqTypeAll = database(idxData).Seq_random_fiftysilence_freqtype;
            recParams.playOrderAll = database(idxData).Seq_random_fiftysilence_playingorder;
            recParams.folderMatAll = database(idxData).Seq_random_fiftysilence_mat;
        case 5
            seq_idx = database(idxData).Seq_fixed_fiftysilence;
            recParams.ISITrialAll = database(idxData).Seq_fixed_fiftysilence_ISI/1000;
            recParams.freqTypeAll = database(idxData).Seq_fixed_fiftysilence_freqtype;
            recParams.playOrderAll = database(idxData).Seq_fixed_fiftysilence_playingorder;
            recParams.folderMatAll = database(idxData).Seq_fixed_fiftysilence_mat;
        case 6
            seq_idx = database(idxData).Seq_fixed_spaced;
            recParams.ISITrialAll = database(idxData).Seq_fixed_spaced_ISI/1000;
            recParams.freqTypeAll = database(idxData).Seq_fixed_spaced_freqtype;
            recParams.playOrderAll = database(idxData).Seq_fixed_spaced_playingorder;
            recParams.folderMatAll = database(idxData).Seq_fixed_spaced_mat;
        case 7
            seq_idx = database(idxData).Seq_fixed_compressed;
            recParams.ISITrialAll = database(idxData).Seq_fixed_compressed_ISI/1000;
            recParams.freqTypeAll = database(idxData).Seq_fixed_compressed_freqtype;
            recParams.playOrderAll = database(idxData).Seq_fixed_compressed_playingorder;
            recParams.folderMatAll = database(idxData).Seq_fixed_compressed_mat;
        case 8
            seq_idx = database(idxData).Seq_completelyrandom;
            recParams.ISITrialAll = database(idxData).Seq_completelyrandom_ISI/1000;
            recParams.freqTypeAll = database(idxData).Seq_completelyrandom_freqtype;
            recParams.playOrderAll = database(idxData).Seq_completelyrandom_playingorder;
            recParams.folderMatAll = database(idxData).Seq_completelyrandom_mat;
        case 9
            seq_idx = database(idxData).Seq_completelyfixed;
            recParams.ISITrialAll = database(idxData).Seq_completelyfixed_ISI/1000;
            recParams.freqTypeAll = database(idxData).Seq_completelyfixed_freqtype;
            recParams.playOrderAll = database(idxData).Seq_completelyfixed_playingorder;
            recParams.folderMatAll = database(idxData).Seq_completelyrandom_mat;
end

% Finding which type of subject it is: for Gamze, if it has 3 numbers, it
% is an a) recording. If it has 4 numbers, it will be the b), c)
% recordings. For example: subject_rec = 609 --> a recording. subject_rec =
% 6091 --> b recording. subject_rec = 6092 --> c recording.
num_final = 0;
total_idx = 0;

if subject_rec ~= subject
    x =num2str(subject_rec);
    num_str = x(end);
    num_final = str2double(num_str);
end

% This logic does not work for the recordings of this two kind of subjects
if subject == 602 || subject == 603
    num_final = 0;
end

% Variable initialization. 
recParams.ISITrial =[];
recParams.playOrder = [];
recParams.folderMat = [];
recParams.freqType = [];

% Taking the first, or b, c sequences.
for i=1:length(seq_idx)
    
    string = cell2mat(seq_idx(i));    
    if num_final == 0
        % Only taking initial recordings.
        if ~contains(string,'b') && ~contains(string,'c')
            id = strfind(string,'t');
            total_idx = total_idx + 1;
            
            % There is a +1 in the indexes because the FRA was a T0 recording
            recParams.idxSequence(total_idx) = str2num(string(id+1:end))+1; %#ok<ST2NM>
            recParams.ISITrial = [recParams.ISITrial, recParams.ISITrialAll(i)];
            recParams.playOrder = [recParams.playOrder, recParams.playOrderAll(i)];
            recParams.folderMat{total_idx} = recParams.folderMatAll{i};
            
            if type_sequence~=1
                recParams.freqType = [recParams.freqType, recParams.freqTypeAll(i)];
            end
        end
    elseif num_final == 1
       % Only taking b recordings.
       if contains(string,'b')
            id = strfind(string,'t');
            total_idx = total_idx + 1;
            
            % There is a +1 in the indexes because the FRA was a T0 recording
            recParams.idxSequence(total_idx) = str2num(string(id+1:end))+1; %#ok<ST2NM>
            recParams.ISITrial = [recParams.ISITrial, recParams.ISITrialAll(i)];
            recParams.playOrder = [recParams.playOrder, recParams.playOrderAll(i)];
            recParams.folderMat{total_idx} = recParams.folderMatAll{i};
            
            if type_sequence~=1
                recParams.freqType = [recParams.freqType, recParams.freqTypeAll(i)];
            end
        end
    elseif num_final == 2
        % Only taking c recordings.
       if contains(string,'c')
            id = strfind(string,'t');
            total_idx = total_idx + 1;
            
            % There is a +1 in the indexes because the FRA was a T0 recording
            recParams.idxSequence(total_idx) = str2num(string(id+1:end))+1; %#ok<ST2NM>
            recParams.ISITrial = [recParams.ISITrial, recParams.ISITrialAll(i)];
            recParams.playOrder = [recParams.playOrder, recParams.playOrderAll(i)];
            recParams.folderMat{total_idx} = recParams.folderMatAll{i};
            
            if type_sequence~=1
                recParams.freqType = [recParams.freqType, recParams.freqTypeAll(i)];
            end
       end
    end
end

% A lot of exceptions, in which the extra recordings (not the initial ones)
% have some recParams that do not match with the initial ones and need to
% be adequated by hand.
%     - 6031: 603, recording b.
%     - 6021, 602, recording b.
%     - 608.
%     - 60900: 609, high sequence.
%     - 60930: 609, wide sequence.
%     - 61300: 613, block VIII of the recording.

if subject_rec == 6031 && type_sequence~=1
    recParams.ISITrial = recParams.ISITrialAll(5:end); %#ok<*NODEF>
    recParams.playOrder = recParams.playOrderAll(5:end);
    recParams.idxSequence = recParams.idxSequence(5:end)-21;
    recParams.freqType = recParams.freqTypeAll(5:end);
    for k = 5:length(recParams.freqTypeAll)
        recParams.folderMat{k-4} = recParams.folderMatAll{k};
    end
end

if subject_rec == 6021 && type_sequence~=1
    recParams.ISITrial = recParams.ISITrialAll(5:end);
    recParams.playOrder = recParams.playOrderAll(5:end);
    recParams.idxSequence = recParams.idxSequence(5:end)-20;
    recParams.freqType = recParams.freqTypeAll(5:end);
    for k = 5:length(recParams.freqTypeAll)
        recParams.folderMat{k-4} = recParams.folderMatAll{k};
    end
end

% Problem with the recording.
if subject_rec == 608 && type_sequence~=1
    recParams.idxSequence (2:end) = recParams.idxSequence(2:end)-1;
end

if subject_rec == 60900 && type_sequence~=6
    recParams.idxSequence = recParams.idxSequence-1;
    recParams.freqType = recParams.freqTypeAll(3:end);
    for k = 3:length(recParams.freqTypeAll)
        recParams.folderMat{k-2} = recParams.folderMatAll{k};
    end
end

if subject_rec == 60930 && type_sequence~=6
    recParams.idxSequence = recParams.idxSequence-1;
    recParams.freqType = recParams.freqTypeAll(5:end);
    for k = 5:length(recParams.freqTypeAll)
        recParams.folderMat{k-4} = recParams.folderMatAll{k};
    end
end

if subject_rec == 61340 || subject_rec == 61440
    recParams.idxSequence = recParams.idxSequence(1)-1;
    recParams.folderMat{1} = recParams.folderMatAll{end};
    if type_sequence~=1
        recParams.freqType = recParams.freqTypeAll(end);
    end
end
    
% In recordings of 607 and 609 there is a mismatch between the planned
% recordings and the ones stored in the computer (we are missing
% 607_c_20,28 and 609_c_37,45)

% Something will have to be done about that when we investigate that data.