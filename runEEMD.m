function [imfsCell, residsCell, numImfs] = runEEMD(data,fs,numEnsemb,maxNumImfs)

if size(data,1) > size(data,2)
    data = data';
end

numChans = size(data,1);
numSamples = size(data,2);

if nargin < 3
    numEnsemb = 20;
end
if nargin < 4
    maxNumImfs = 20;
end

imfsCell = cell(numChans,1);
residsCell = cell(numChans,1);
numImfs = zeros(numChans,1);
for i = 1:numChans
    % adding noises
    wns = zeros(numEnsemb,numSamples);
    preImfs = zeros(maxNumImfs,numSamples);
    preResids = zeros(1,numSamples);
    for j = 1:numEnsemb
        wns(j,:) = 0.2*std(data(i,:))*randn(numSamples,1);
        data_wNoise = data(i,:) + wns(j,:);
        
        [tempImf,tempResid] = emd(data_wNoise,'MaxNumIMF',maxNumImfs);
        preImfs(1:size(tempImf,2),:) = preImfs(1:size(tempImf,2),:) + tempImf';
        preResids = preResids + tempResid';
    end
    
    imfsCell{i} = preImfs / numEnsemb;
    numImfs(i) = length(find(any(imfsCell{i},2)));
    imfsCell{i} = imfsCell{i}(any(imfsCell{i},2),:);
    residsCell{i} = preResids / numEnsemb;
end