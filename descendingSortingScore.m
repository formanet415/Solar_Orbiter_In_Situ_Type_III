function sortingScore = descendingSortingScore(arr)
    % Calculate the number of inversions (detecting dust in the frequency domain)
    numInversions = 0;
    n = length(arr);
    
    for i = 1:n
        for j = i+1:n
            if arr(i) < arr(j)
                numInversions = numInversions + 1;
            end
        end
    end
    
    % Calculate the sorting score
    sortingScore = 1 - (numInversions / (n * (n - 1) / 2));
end
