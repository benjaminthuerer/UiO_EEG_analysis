% permutation Test:
% permutation test for 1D or 2D data. Please note that in case of 2D data,
% the p-values will be corrected for multiple comparisons using the maximum
% permutated statistics across the 2d dimenson (maximum statistics)
% 
% [p_value,obsDiff] = permutationTest(mData,groupA,groupB,n)
% 
% input:
% mData: is a one-dimensional vector of all data or a matrix. If mData is a
%        matrix, the first dimension has to be the group dimension.
% groupA: are the indices of mData which belong to the first data set
% groupB: are the indices of mData which belong to the second data set
% n: is the number of permutations which should be performed
% correction: do you want to correct (1) or not correct (0) for multiple
%             comparisons? (default = 1)
% 
% output:
% p_value is the p-value of the test
% obsDiff is the signed difference of (groupA-groupB)
% ciThresh is the 95% confidence interval over all permutations (as
% correction for multiple comparisons)
% medianThresh is the median over all permutations
% stdThresh is the standarddeviation over all permutations
%
% written by Benjamin Thürer (benjamin.thuerer@kit.edu)

function [p_value, obsDiff, ciThresh, medianThresh, stdThresh] = permutationTest(mData, groupA, groupB, n, correction, cluster)

if nargin < 4
    disp('number of permutations not provided. Taking default = 5000')
    n = 5000;
    correction = 1;
    cluster = 0;
elseif nargin < 5
    correction = 1;
    cluster = 0;
elseif nargin < 6
    cluster = 0;
end

sizeInfo = size(mData);

if sizeInfo(1) == 1 || sizeInfo(2) == 1
    obsDiff = nanmedian(mData(groupA))-nanmedian(mData(groupB));
    singlePermutation = 1;
else
    obsDiff = nanmedian(mData(groupA,:))-nanmedian(mData(groupB,:));
    singlePermutation = 0;
end

percentages = n/10:n/10:n;
percentages2 = 10:10:100;
i = 1;
permDiff = [];
while i < n+1
    if singlePermutation == 1
        permI = randperm(length(mData));
        permA = permI(1:length(groupA));
        permB = permI(length(groupA)+1:end);  
        permDiff(i) = nanmedian(mData(permA))-nanmedian(mData(permB));
    elseif singlePermutation == 0 && correction == 1
        permI = randperm(size(mData,1));
        permA = permI(1:length(groupA));
        permB = permI(length(groupA)+1:end); 
        permDiff(i) = max(abs(nanmedian(mData(permA,:),1)-nanmedian(mData(permB,:),1)));
    else
        permI = randperm(size(mData,1));
        permA = permI(1:length(groupA));
        permB = permI(length(groupA)+1:end); 
        permDiff(i,:) = nanmedian(mData(permA,:),1)-nanmedian(mData(permB,:),1);
    end
    if ~isempty(find(percentages == i)) && cluster == 0
        disp([num2str(percentages2(percentages == i)) ' % done']);
    end
    i = i+1;
end

if singlePermutation == 1
    p_value = sum(abs(permDiff) > abs(obsDiff)) / n;
elseif singlePermutation == 0 && correction == 1
    p_value = sum(abs(repmat(permDiff',1,size(obsDiff,2))) > abs(repmat(obsDiff,n,1))) / n;
else
    p_value = sum(abs(permDiff) > abs(repmat(obsDiff,n,1))) / n;
end

if singlePermutation == 0
    permDiff = sort(reshape(permDiff,1,numel(permDiff)));
end
ciThresh = [permDiff(round(n/100*2.5)),permDiff(round(n/100*97.5))];
medianThresh = nanmedian(permDiff);
stdThresh = nanstd(permDiff);
obsDiff = squeeze(obsDiff);
end
    
