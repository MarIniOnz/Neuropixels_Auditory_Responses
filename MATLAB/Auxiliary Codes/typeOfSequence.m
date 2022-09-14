function string = typeOfSequence(num_type)
% 
% Function that gets a number referring to a type of sequence and returns
% the corresponding string.
% 
% Input:
%     num_type: num, index of type of sequence.
% Output: 
%     string: str, name of the type of sequence.

switch num_type
    case 1
        string = 'FRA';
    case 2
        string = 'Fixed Random';
    case 3
        string = 'Random fifty';
    case 4
        string = 'Fixed fifty';
    case 5
        string = 'Fixed fifty silence';
    case 6
        string = 'Fixed spaced';
    case 7
        string = 'Fixed compressed';
    case 8
        string = 'Completely random';
    case 9
        string = 'Completely fixed';
end