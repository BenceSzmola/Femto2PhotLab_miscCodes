function output = artSuppMaster(inputStruct)

% two ways to call this function:
% inputStruct = artSuppMaster(1) this will initialize and return the input structure the function expects in
%                                the second way to call it
% output = artSuppMaster(inputStruct) this is the call which will then run the algorithms

%% check which type of call was made
if nargin == 0 || ~isstruct(inputStruct)
    % initialize input struct
    output.data = [];
    output.tAxis = [];
    output.fs = [];
    output.data_perio = [];
    output.doPeriodFilt = [];
    output.decompType = [];
    output.flagType = [];
    output.bssType = [];
    output.slideWinSize = [];
    output.iecThr = [];
    output.autoCorrThr = [];
    output.autoCorrLagNum = [];
    
    return
else
    % extract parameters from input
    data = inputStruct.data;
    tAxis = inputStruct.tAxis;
    fs = inputStruct.fs;
    data_perio = inputStruct.data_perio;
    doPeriodFilt = inputStruct.doPeriodFilt;
    decompType = inputStruct.decompType;
    flagType = inputStruct.flagType;
    bssType = inputStruct.bssType;
    slideWinSize = inputStruct.slideWinSize;
    iecThr = inputStruct.iecThr;
    autoCorrThr = inputStruct.autoCorrThr;
    autoCorrLagNum = inputStruct.autoCorrLagNum;
end

%% handling the input data
if isempty(data)
    rhdStruct = read_Intan_RHD2000_file_szb;
    data = rhdStruct.amplifier_data;
    fs = rhdStruct.fs;
    tAxis = rhdStruct.tAxis;
end

if size(data,1) > size(data,2)
    data = data';
end

%% extract sample and channel number
numSamples = size(data,2);
numChans = size(data,1);

%% run periodic filter if requested
if doPeriodFilt
    data = periodicNoise(data,fs);
elseif ~isempty(data_perio)
    data = data_perio;
end

%% running the requested decomposition
switch decompType
    case 'EEMD'
        %% compute EEMD 
        [compCell, residsCell, numImfs] = runEEMD(data,fs,20,20);
        
    case 'DWT'
        
        %% compute DWT decomposition
        wname = 'db4';
        freqLow = 4;
        decLvl = round(log2(fs/freqLow));
        compCell = cell(numChans,decLvl+1);
        for chan = 1:numChans
            [C,L] = wavedec(data(chan,:),decLvl,wname);
            details = detcoef(C,L,1:decLvl);
            approx = appcoef(C,L,wname);
            compCell(chan,:) = [details, mat2cell(approx,ones(numChans,1))];
        end
        
    case 'SWT'
        %% compute SWT decomposition
        wname = 'db4';
        freqLow = 4;
        decLvl = round(log2(fs/freqLow));
        padLen = numSamples/(2^decLvl);
        padLen = (2^decLvl)*round(padLen) - numSamples;
        dataPad = [data, zeros(numChans,padLen)];
        compCell = cell(numChans,1);
        for chan = 1:numChans
            compCell{chan} = swt(dataPad(chan,:),decLvl,wname);
        end
end
compCellClean = compCell;

%% flag suspicious segments
switch flagType
    case 'IEC'
        %% flagging with inter-electrode-correlation
        %need two branches because of DWT's unequal component lengths
        switch decompType
            case {'EEMD','SWT'}
                if strcmp(decompType,'EEMD')
                    numComps = min(numImfs);
                else
                    numComps = decLvl+1;
                end
                IEC = compWiseIEC(compCell,numComps,slideWinSize);
                flaggedInds = IEC > iecThr;
                
            case 'DWT'
                IEC = compWiseIEC_dwt(compCell,decLvl+1,slideWinSize);                
                flaggedInds = cellfun(@(x) x > iecThr, IEC, 'UniformOutput', false);
        end
        
    case 'autoCorr'
        %% flagging based on autocorrelation        
        %need two branches because of DWT's unequal component lengths
        switch decompType
            case {'EEMD','SWT'}                
                flaggedInds = cell(size(compCell));
                for chan = 1:numChans
                    for compNum = 1:size(compCell{chan},1)
                        temp = slideAutoCorr(compCell{chan}(compNum,:),slideWinSize,autoCorrLagNum);
                        temp = temp < autoCorrThr;
                        flaggedInds{chan}(compNum,:) = temp;
                    end
                end
                
            case 'DWT'
                flaggedInds = cell(size(compCell));
                for chan = 1:numChans
                    for compNum = 1:(decLvl+1)
                        if compNum < (decLvl+1)
                            slideWinSize_mod = slideWinSize / (2^compNum);
                        else
                            slideWinSize_mod = slideWinSize / (2^decLvl);
                        end
                        temp = slideAutoCorr(compCell{chan,compNum},slideWinSize_mod,autoCorrLagNum);
                        temp = temp < autoCorrThr;
                        flaggedInds{chan,compNum} = temp;
                    end
                end
                
        end
        
        
end

%% deal with flagged indices
% either directly on the decomposition components or first use a BSS algorithm and then discard bad sources
% for now dont allow autocorr flagging and BSS
if strcmp(flagType,'autoCorr') && ~strcmp(bssType,'')
    bssType = '';
    warndlg('You cant use autoCorr with a BSS, proceeding without BSS...')
end

switch bssType
    case 'CCA'
        %% use CCA
        switch decompType
            case {'EEMD','SWT'}
                for compNum = 1:numComps % numComps was defined back at flagging stage
                    currComps = cell2mat(cellfun(@(x)x(compNum,:), compCell, 'UniformOutput',false));
                    currCompsCl = currComps;
                    
                    segments = extractLogicSegments(flaggedInds(compNum,:),slideWinSize/2);
                    for segNum = 1:size(segments,1)
                        segInds = segments(segNum,1):segments(segNum,2);
                        currCompSegs = currComps(:,segInds);
                        reconstr = doCCAdiscard(currCompSegs,autoCorrThr,autoCorrLagNum);
                        currCompsCl(:,segInds) = reconstr;
                    end
                    
                    for chan = 1:numChans
                        compCellClean{chan}(compNum,:) = currCompsCl;
                    end
                end
                
            case 'DWT'
                for compNum = 1:decLvl+1
                    currComps = cell2mat(compCell(:,compNum));
                    currCompsCl = currComps;
                    
                    if compNum < decLvl+1
                        slideWinSize_mod = slideWinSize / (2^compNum);
                    else
                        slideWinSize_mod = slideWinSize / (2^decLvl);
                    end
                    segments = extractLogicSegments(flaggedInds{compNum},slideWinSize_mod/2);
                    for segNum = 1:size(segments,1)
                        segInds = segments(segNum,1):segments(segNum,2);
                        currCompSegs = currComps(:,segInds);
                        reconstr = doCCAdiscard(currCompSegs,autoCorrThr,autoCorrLagNum);
                        currCompsCl(:,segInds) = reconstr;
                    end
                    
                    for chan = 1:numChans
                        compCellClean{chan,compNum} = currCompsCl;
                    end
                end
                
        end
    case 'ICA'
        %% use ICA
        
    otherwise
        %% dont use BSS
        for chan = 1:numChans
            switch decompType
                case {'EEMD','SWT'}
                    switch flagType
                        case 'IEC'
                            temp = compCell{chan};
                            temp(flaggedInds) = 0;
                            compCellClean{chan} = temp;
                            
                        case 'autoCorr'
                            temp = compCell{chan};
                            for compNum = 1:size(temp,1)
                                temp(compNum,flaggedInds{chan}(compNum,:)) = 0;
                            end
                            compCellClean{chan} = temp;

                    end

                case 'DWT'
                    switch flagType
                        case 'IEC'
                            temp = compCell(chan,:);
                            for compNum = 1:length(temp)
                                temp{compNum}(flaggedInds{compNum}) = 0;
                            end
                            compCellClean(chan,:) = temp;
                            
                        case 'autoCorr'
                            temp = compCell(chan,:);
                            for compNum = 1:length(temp)
                                temp{compNum}(flaggedInds{chan,compNum}) = 0;
                            end
                            compCellClean(chan,:) = temp;
                    end

            end
        end
end

%% time to reconstruct
dataCleaned = zeros(size(data));
switch decompType
    case 'EEMD'
        dataCleaned = cell2mat(cellfun(@sum, compCellClean, 'UniformOutput', false)) + cell2mat(residsCell);
        
    case 'SWT'
        for chan = 1:numChans
            reconstrPad = iswt(compCellClean{chan},wname);
            dataCleaned(chan,:) = reconstrPad(1:numSamples);
        end
        
    case 'DWT'
        % reconstructing vector c for dwt function
        for chan = 1:numChans
            cClean = zeros(1,sum(L(1:end-1)));
            for i = 1:length(L)-1
                if i == 1
                    cClean(1:L(1)) = compCellClean{chan,decLvl+1};
                else
                    cClean(sum(L(1:i-1))+1:sum(L(chan,1:i))) = compCellClean{chan,length(L)-i};
                end
            end
            dataCleaned(chan,:) = waverec(cClean,L,wname);
        end
        
end

plotBefAft(tAxis,fs,data,dataCleaned,false,decompType,flagType,bssType)

output = dataCleaned;
% end of main function
end

function IEC = compWiseIEC(compCell,numComps,slideWinSize)
    % create IEC placeholder
    numSamples = size(compCell{1},2);
    IEC = zeros(numComps,numSamples);
    
    for compNum = 1:numComps
        currComps = cell2mat(cellfun(@(x) x(compNum,:), compCell,'UniformOutput', false));
        for i = 1:slideWinSize:numSamples
            if (i + slideWinSize - 1) > numSamples
                win = i:numSamples;
            else
                win = i:(i + slideWinSize - 1);
            end
            corrVals = corrcoef(currComps(:,win)');
            corrVals = tril(corrVals,-1);
            corrVals(corrVals == 0) = [];
            IEC(compNum,win) = mean(corrVals);
        end
    end
end

function IEC = compWiseIEC_dwt(compCell,numComps,slideWinSize)
    % create the IECs placeholder, here a cell is needed because of the unequal lengths of components
    IEC = cell(numComps,1);
    
    for compNum = 1:numComps
        currComps = cell2mat(compCell(:,compNum));
        
        currCompsLen = size(currComps,2);
        
        if compNum < numComps
            slideWinSizeComp = round(slideWinSize / (2^compNum));
        else
            slideWinSizeComp = round(slideWinSize / (2^(compNum-1)));
        end
        
        for i = 1:slideWinSizeComp:currCompsLen
            if (i + slideWinSizeComp -1) > currCompsLen
                win = i:currCompsLen;
            else
                win = i:(i + slideWinSizeComp -1);
            end
            
            corrVals = corrcoef(currComps(:,win)');
            corrVals = tril(corrVals,-1);
            corrVals(corrVals == 0) = [];
            IEC{compNum}(win) = mean(corrVals);
        end
    end
end

function compAutoCorr = slideAutoCorr(comp,slideWinSize,lagNum)
    compAutoCorr = zeros(size(comp));
    for i = 1:slideWinSize:length(comp)
        if (i + slideWinSize - 1) > length(comp)
            win = i:length(comp);
        else
            win = i:(i + slideWinSize - 1);
        end
        
        r = xcorr(comp,lagNum,'coeff');
        compAutoCorr(win) = r(end);
    end
end

function segments = extractLogicSegments(inputVec,minSegLen)
    segments = [];
    % find indexes above the threshold for a given amount of time
    aboveThrInds = find(inputVec);

    if isempty(aboveThrInds)
        return
    end
    indDiffs = diff(aboveThrInds);
    disconts = find(indDiffs ~= 1);

    segments = [aboveThrInds([1, disconts+1])', aboveThrInds([disconts, length(aboveThrInds)])'];

    segments(segments(:,2)-segments(:,1) < minSegLen,:) = [];

end

function reconstr = doCCAdiscard(inputMat,discardThr,lagNum)
    % compute CCA, discard sources with the autocorrelations lower than threshold
    xDelay = [zeros(size(inputMat,1),lagNum), inputMat(:,1:end-lagNum)];
    [A,~,r,U,~,~] = canoncorr(inputMat',xDelay');
    sources2discard = r < discardThr;

    Amod = A;
    Amod(:,sources2discard) = 0;
    reconstr = (U*pinv(Amod))';

end

function plotBefAft(tAxis,fs,data,dataCl,spectro,decompType,flagType,bssType)
    for chan = 1:size(data,1)
        if ~spectro
            figure('Name',sprintf('Channel #%d - Before and after cleaning',chan));
            subplot(211)
            plot(tAxis,data(chan,:))
            title(sprintf('Ch#%d - Raw',chan))
            xlabel('Time [s]')
            ylabel('Voltage [\muV]')
            
            subplot(212)
            plot(tAxis,dataCl(chan,:))
            title(sprintf('Ch#%d - %s, %s, %s',chan,decompType,flagType,bssType))
            xlabel('Time [s]')
            ylabel('Voltage [\muV]')
            
            linkaxes(findobj(gcf,'Type','axes'),'xy')
        else
            figure('Name',sprintf('Channel #%d CWT - Before and after cleaning',chan));
            [cfs,f] = cwt(data(chan,:),'amor',fs,'FrequencyLimits',[1,1000]);
            subplot(211)
            imagesc('XData',tAxis,'YData',f,'CData',abs(cfs))
            title(sprintf('Ch#%d CWT - Raw',chan))
            xlabel('Time [s]')
            ylabel('Frequency [Hz]')
            c = colorbar;
            c.Label.String = 'CWT coeff. magnitude';
            clear cfs f
            
            [cfs,f] = cwt(dataCl(chan,:),'amor',fs,'FrequencyLimits',[1,1000]);
            subplot(212)
            imagesc('XData',tAxis,'YData',f,'CData',abs(cfs))
            title(sprintf('Ch#%d CWT - %s, %s, %s',chan,decompType,flagType,bssType))
            xlabel('Time [s]')
            ylabel('Frequency [Hz]')
            c = colorbar;
            c.Label.String = 'CWT coeff. magnitude';
            clear cfs f
            
            linkaxes(findobj(gcf,'Type','axes'),'x')
        end
    end
end