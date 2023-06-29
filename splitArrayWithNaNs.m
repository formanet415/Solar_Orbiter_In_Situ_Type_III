function [splitArray, rowLengths] = splitArrayWithNaNs(inputArray)
    % Find indices of non-NaN values
    nonNaNIndices = find(~isnan(inputArray));
    
    % Find indices where row transition occurs
    rowTransitionIndices = find(diff(nonNaNIndices) > 1);
    
    % Determine row lengths
    rowLengths = diff([0, rowTransitionIndices, numel(nonNaNIndices)]);
    
    % Preallocate splitArray based on row lengths
    splitArray = NaN(numel(rowLengths), max(rowLengths));
    
    % Fill splitArray with non-NaN values
    startIndex = 1;
    for i = 1:numel(rowLengths)
        endIndex = startIndex + rowLengths(i) - 1;
        splitArray(i, 1:rowLengths(i)) = inputArray(nonNaNIndices(startIndex:endIndex));
        startIndex = endIndex + 1;
    end
end
