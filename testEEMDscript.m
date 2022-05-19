clear fname path imfs resids data numCh numDataPoints

[fname,path] = uigetfile('*.rhd');
rhdStruct = read_Intan_RHD2000_file_szb([path,fname]);
fs = rhdStruct.fs;
data = rhdStruct.amplifier_data;
chanSubset = [7:14];
data = data(chanSubset,:);
data = periodicNoise(data,rhdStruct.fs);
numCh = size(data,1);
numDataPoints = size(data,2);

numEnsemb = 20;
maxNumImfs = 20;

imfs = cell(numCh,1);
resids = cell(numCh,1);
emdInfo = cell(numCh,1);
for i = 1:numCh
    % adding noises
    wns = zeros(numEnsemb,numDataPoints);
    preImfs = zeros(maxNumImfs,numDataPoints);
    preResids = zeros(1,numDataPoints);
    for j = 1:numEnsemb
        wns(j,:) = 0.2*std(data(i,:))*randn(numDataPoints,1);
        data_wNoise = data(i,:) + wns(j,:);
        
        [tempImf,tempResid] = emd(data_wNoise,'MaxNumIMF',maxNumImfs);
        preImfs(1:size(tempImf,2),:) = preImfs(1:size(tempImf,2),:) + tempImf';
        preResids = preResids + tempResid';
    end
    
    imfs{i} = preImfs / numEnsemb;
    numImfs = length(find(any(imfs{i},2)));
    imfs{i} = imfs{i}(any(imfs{i},2),:);
    resids{i} = preResids / numEnsemb;
    
    numFigRows = ceil((numImfs+2)/2);
    
    figure('Name',['EEMD - Channel #',num2str(i),' IMFs']);
    subplot(numFigRows,2,1)
    plot(data(i,:))
    title('Raw data')
    for j = 1:numImfs
        subplot(numFigRows,2,j+1)
        plot(imfs{i}(j,:))
        title(['IMF#',num2str(j),'/',num2str(numImfs)])
    end
    subplot(numFigRows,2,numImfs+2)
    plot(resids{i})
    title('Residual')
    ax2synch = findobj('Parent',gcf,'Type','axes');
    linkaxes(ax2synch,'x')
    
end

clear ax2synch i j numCh numDataPoints numICs numFigRows