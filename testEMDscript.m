% clear fname path imfs resids data numCh numDataPoints

% [fname,path] = uigetfile('*.rhd');
% rhdStruct = read_Intan_RHD2000_file_szb([path,fname]);
% data = rhdStruct.amplifier_data;
% chanSubset = [7:14];
% data = data(chanSubset,:);
numCh = size(data,1);
numDataPoints = size(data,2);
fs = rhdStruct.fs;

maxNumImfs = 20;

imfs = cell(numCh,1);
resids = cell(numCh,1);
emdInfo = cell(numCh,1);
imfICs = cell(numCh,1);
for i = 1:numCh
    [tempImf,tempResid,emdInfo{i}] = emd(data(i,:),'SiftRelativeTolerance',0.05,'MaxNumIMF',maxNumImfs);
    numImfs = length(emdInfo{i}.NumIMF);
    imfs{i} = tempImf';
    imfs{i} = imfs{i}(any(imfs{i},2),:);
    resids{i} = tempResid';
    clear tempImf tempResid
    
    numFigRows = ceil((numImfs+2)/2);
    
    figure('Name',['EMD - Channel #',num2str(i),' IMFs']);
    subplot(numFigRows,2,1)
    plot(data(i,:))
    title('Raw data')
    for j = 1:numImfs
        subplot(numFigRows,2,j + 1)
        plot(imfs{i}(j,:))
        title(['IMF#',num2str(j),'/',num2str(numImfs)])
    end
    subplot(numFigRows,2,j + 2)
    plot(resids{i})
    title('Residual')
    ax2synch = findobj('Parent',gcf,'Type','axes');
    linkaxes(ax2synch,'x')
    
end

clear ax2synch i j numCh numDataPoints numICs