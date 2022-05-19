function output = artSuppMaster(inputStruct)

% two ways to call this function:
% inputStruct = artSuppMaster(1) this will initialize and return the input structure the function expects in
%                                the second way to call it
% output = artSuppMaster(inputStruct) this is the call which will then run the algorithms

%% check which type of call was made
if ~isstruct(inputStruct)
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
            [c,l] = wavedec(data(chan,:),decLvl,wname);
            details = detcoef(c,l,1:decLvl);
            approx = appcoef(c,l,wname);
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
        compCell = cell(numChans,decLvl+1);
        for chan = 1:numChans
             temp = swt(dataPad(chan,:),decLvl,wname);
            compCell(chan,:) = mat2cell(temp,ones(decLvl+1,1))';
        end
end

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
        
        lagNum = 10;
        
        %need two branches because of DWT's unequal component lengths
        switch decompType
            case {'EEMD','SWT'}
                
                flaggedInds = cell(size(compCell));
                for chan = 1:numChans
                    for compNum = 1:size(compCell{chan},1)
                        temp = slideAutoCorr(compCell{chan}(compNum,:),slideWinSize,lagNum);
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
                        temp = slideAutoCorr(compCell{chan,compNum},slideWinSize_mod,lagNum);
                        temp = temp < autoCorrThr;
                        flaggedInds{chan,compNum} = temp;
                    end
                end
                
        end
        
        
end

%% deal with flagged indices
% either directly on the decomposition components or first use a BSS algorithm and then discard bad sources

switch bssType
    case 'CCA'
        %% use CCA
        
    case 'ICA'
        %% use ICA
        
    otherwise
        %% dont use BSS
        switch decompType
            case {'EEMD','SWT'}
                compCellClean = compCell;
                
                switch flagType
                    case 'IEC'
                        
                        
                    case 'autoCorr'
                        
                        
                end
                flaggedIndsUnroll = cell2mat(flaggedInds);
                compCellCleanUnroll = cell2mat(compCellClean);
                
                compCellCleanUnroll(flaggedIndsUnroll) = 0;
                
                compCellClean
                
            case 'DWT'
                
                
        end
end

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