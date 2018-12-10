function [peptideLengthHeight, peptideLength, confidenceMatHighHeight, confidenceMatHigh, confidenceMatLowHeight, confidenceMatLow] = woodsPlot(data,confidenceInterval,timepoint)
% Passes timepoint series into function and outputs list of peptides and
% exposure for all global list, deprotected significant peptides and exposure,
% and protected significant peptides and exposure

exposureAll = data.data(:,4:end);
exposureSize = size(exposureAll',1);
exposureSum = [exposureAll,sum(exposureAll,2)];
exposureMax = max(max(abs(exposureAll)));
exposure = exposureSum(:,timepoint);

peptideStart = [];
for i=1:length(data.data)
    temp = [data.data(i,1):data.data(i,2)];
    peptideStart{i,1} = temp;
end

exposureMat = cellfun(@(x) (x*0)+1,peptideStart,'un',0);

for i=1:length(exposure)
    n = exposure(i);
    exposureMat(i,1) = cellfun(@(x) x*n,exposureMat(i,1),'un',0);
end

longest = [];
for i=1:length(peptideStart)
    num = length(peptideStart{i,1});
    longest = [longest, num];
end

peptideLength = zeros(length(peptideStart),max(longest));

every_n = [1:2:length(peptideStart)*2];

for i=1:length(peptideStart)
    peptideLength(every_n(i),1:longest(i)') = cell2mat(peptideStart(i,1));
    peptideLength(every_n(i)+1,1:longest(i)') = cell2mat(exposureMat(i,1));
end

peptideLength(peptideLength == 0) = NaN;

peptideLengthHeight = size(peptideLength,1);

confidenceMatHigh = [];
confidenceMatLow = [];

for i=2:2:(peptideLengthHeight)
    if sum(confidenceInterval <= peptideLength(i,:)) > 0
        confidenceMatHigh = [confidenceMatHigh; peptideLength(i-1:i,:)];
    end

    if sum(-confidenceInterval >= peptideLength(i,:)) > 0
        confidenceMatLow = [confidenceMatLow; peptideLength(i-1:i,:)];
    end
end

confidenceMatHighHeight = size(confidenceMatHigh,1);
confidenceMatLowHeight = size(confidenceMatLow,1);