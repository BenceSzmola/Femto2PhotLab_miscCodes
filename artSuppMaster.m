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
    output.doPeriodFilt = false;
    output.decompType = 'SWT';
    output.freqLow = 4;
    output.flagType = 'IEC';
    output.bssType = '';
    output.slideWinSize = 2000;
    output.iecThr = 0.8;
    output.autoCorrThr = 0.99;
    output.autoCorrLagNum = 10;
    output.doRawPlot = false;
    output.doSpectroPlot = false;
    
    return
else
    % extract parameters from input
    data = inputStruct.data;
    tAxis = inputStruct.tAxis;
    fs = inputStruct.fs;
    data_perio = inputStruct.data_perio;
    doPeriodFilt = inputStruct.doPeriodFilt;
    decompType = inputStruct.decompType;
    freqLow = inputStruct.freqLow;
    flagType = inputStruct.flagType;
    bssType = inputStruct.bssType;
    slideWinSize = inputStruct.slideWinSize;
    iecThr = inputStruct.iecThr;
    autoCorrThr = inputStruct.autoCorrThr;
    autoCorrLagNum = inputStruct.autoCorrLagNum;
    doRawPlot = inputStruct.doRawPlot;
    doSpectroPlot = inputStruct.doSpectroPlot;
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
%         freqLow = 4;
        decLvl = round(log2(fs/freqLow));
        compCell = cell(numChans,decLvl+1);
        for chan = 1:numChans
            [C,L] = wavedec(data(chan,:),decLvl,wname);
            details = detcoef(C,L,1:decLvl);
            approx = appcoef(C,L,wname);
            compCell(chan,:) = [details, mat2cell(approx,ones(1,1))];
        end
        
    case 'SWT'
        %% compute SWT decomposition
        wname = 'db4';
%         freqLow = 4;
        decLvl = round(log2(fs/freqLow));
        padLen = numSamples/(2^decLvl);
        padLen = (2^decLvl)*ceil(padLen) - numSamples;
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
    case {'CCA','ICA'}
        %% use CCA or ICA
        
        % create the IC discarding figure if ICA is selected
        if strcmp(bssType,'ICA')
            [fig,spH,dbH] = makeICdiscardFig(numChans);
        end
        
        switch decompType
            case {'EEMD','SWT'}
                for compNum = 1:numComps % numComps was defined back at flagging stage
                    currComps = cell2mat(cellfun(@(x)x(compNum,:), compCell, 'UniformOutput',false));
                    currCompsCl = currComps;
                    
                    segments = extractLogicSegments(flaggedInds(compNum,:),slideWinSize/2);
                    for segNum = 1:size(segments,1)
                        segInds = segments(segNum,1):segments(segNum,2);
                        currCompSegs = currComps(:,segInds);
                        
                        % apply either CCA or ICA
                        if strcmp(bssType,'CCA')
                            reconstr = doCCAdiscard(currCompSegs,autoCorrThr,autoCorrLagNum);
                        elseif strcmp(bssType,'ICA')
                            segtAxis = tAxis(segInds);
                            [reconstr,skip] = doICAdiscard(currCompSegs,decompType,compNum,segNum,size(segments,1),segtAxis,fig,spH,dbH);
                            if skip
                                break
                            end
                        end
                        
                        currCompsCl(:,segInds) = reconstr;
                    end
                    
                    for chan = 1:numChans
                        compCellClean{chan}(compNum,:) = currCompsCl(chan,:);
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
                        
                        if strcmp(bssType,'CCA')
                            reconstr = doCCAdiscard(currCompSegs,autoCorrThr,autoCorrLagNum);
                        elseif strcmp(bssType,'ICA')
                            segtAxis = tAxis(segInds);
                            [reconstr,skip] = doICAdiscard(currCompSegs,decompType,compNum,segNum,size(segments,1),segtAxis,fig,spH,dbH);
                            if skip
                                break
                            end
                        end
                        
                        currCompsCl(:,segInds) = reconstr;
                    end
                    
                    for chan = 1:numChans
                        compCellClean{chan,compNum} = currCompsCl(chan,:);
                    end
                end
                
        end
        
        if strcmp(bssType,'ICA') && ishandle(fig)
            delete(fig)
        end
        
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
            for i = 1:decLvl+1
                if i == 1
                    cClean(1:L(1)) = compCellClean{chan,decLvl+1};
                else
                    cClean(sum(L(1:i-1))+1:sum(L(1:i))) = compCellClean{chan,length(L)-i};
                end
            end
            dataCleaned(chan,:) = waverec(cClean,L,wname);
        end
        
end

if doRawPlot
    plotBefAft(tAxis,fs,data,dataCleaned,false,decompType,flagType,bssType)
end
if doSpectroPlot
    plotBefAft(tAxis,fs,data,dataCleaned,true,decompType,flagType,bssType)
end

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
        IEC{compNum} = zeros(1,currCompsLen);
        
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
    if size(inputMat,2) <= lagNum
        lagNum = 1;
    end
    xDelay = [zeros(size(inputMat,1),lagNum), inputMat(:,1:end-lagNum)];
    [A,~,r,U,~,~] = canoncorr(inputMat',xDelay');
    sources2discard = r < discardThr;

    Amod = A;
    Amod(:,sources2discard) = 0;
    reconstr = (U*pinv(Amod))';

end

function [reconstr,skip] = doICAdiscard(inputMat,decompType,compNum,segNum,numSegs,segtAxis,fig,spH,dbH)
% computing ICA and then asking the user to choose which ICs they want to discard
    [ICs,A,W] = fastica(inputMat);
    
    % calling the function which fills in the plots into the discarding figure
    plot2discardICfig(spH,dbH,inputMat,ICs,decompType,compNum,segNum,numSegs,segtAxis)
    % reveal the discarding figure
    fig.Visible = 'on';
    % halt execution until the user is interacting with the figure, they can resume from the figure
    uiwait
    % extract the selections from the figure, then reset the figure's values
    ICs2discard = fig.UserData.ICs2discard;
    fig.UserData.ICs2discard = false(size(inputMat,1));
    % check whether the user selected the option to skip current component
    skip = fig.UserData.skip;
    fig.UserData.skip = false;
    
    % discard selected ICs, then reconstruct using mixing matrix
    ICs(ICs2discard,:) = 0;
    reconstr = A*ICs;
end

function [fig,spH,dbH] = makeICdiscardFig(numChans)
    UD.ICs2discard = false(numChans,1);
    UD.skip = false;
    fig = figure('Name','Discarding ICs',...
                'NumberTitle','off',...
                'WindowState','maximized',...
                'Visible','off',...
                'UserData',UD);
            
    dbH = gobjects(numChans,1);
    for i = 1:numChans
        subplot(numChans,2,(2*i-1),'Parent',fig)
        sp = subplot(numChans,2,2*i,'Parent',fig);
        pos = sp.Position;
        db = uicontrol(fig,'Style','pushbutton',...
                'Units','normalized',...
                'Position',[pos(1)+pos(3)+0.01, pos(2)+pos(4)/4, 0.05, 0.05],...
                'String','Keep',...
                'BackgroundColor','g',...
                'Callback',@discardCB,...
                'Tag',num2str(i));
        dbH(i) = db;
    end
    
    spH = findobj(fig,'Type','axes');
    spH(1:length(spH)) = spH(length(spH):-1:1);
    
    uicontrol(fig,'Style','pushbutton',...
        'Units','normalized',...
        'Position',[0.4, 0.01, 0.05, 0.05],...
        'String','Continue',...
        'Callback',@carryOn);
    uicontrol(fig,'Style','pushbutton',...
        'Units','normalized',...
        'Position',[0.55, 0.01, 0.05, 0.05],...
        'String','Skip comp.',...
        'Callback',@skipComp);
    
    function discardCB(h,~)
        icNum = str2double(h.Tag);
        if fig.UserData.ICs2discard(icNum)
            fig.UserData.ICs2discard(icNum) = false;
            h.BackgroundColor = 'g';
            h.String = 'Keep';
        else
            fig.UserData.ICs2discard(icNum) = true;
            h.BackgroundColor = 'r';
            h.String = 'Discard';
        end
    end

    function skipComp(~,~)
        fig.UserData.skip = true;
        carryOn
    end

    function carryOn(~,~)
        fig.Visible = 'off';
        [dbH.String] = deal('Keep');
        [dbH.BackgroundColor] = deal('g');
        [dbH.Visible] = deal('on');
        [spH.Visible] = deal('on');
        
        for spNum = 1:2*numChans
            cla(spH(spNum))
        end
        uiresume
    end
end

function plot2discardICfig(spH,dbH,inputMat,ICs,decompType,compNum,segNum,numSegs,segtAxis)
    numICs = size(ICs,1);
    for j = 1:size(inputMat,1)
        plot(spH(2*j-1),segtAxis,inputMat(j,:))
        title(spH(2*j-1),sprintf('Ch #%d - %s comp. #%d, flagged segment #%d/%d',j,decompType,compNum,segNum,numSegs))
        if j <= numICs
            plot(spH(2*j),segtAxis,ICs(j,:))
            title(spH(2*j),sprintf('Extracted IC #%d',j))
        else
            spH(2*j).Visible = 'off';
            dbH(j).Visible = 'off';
        end
    end
end

function plotBefAft(tAxis,fs,data,dataCl,spectro,decompType,flagType,bssType)
    if strcmp(bssType,'')
        bssType = 'No BSS';
    end
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
            imagesc(tAxis,log2(f),abs(cfs))
            axis tight
            ax = gca;
            ax.YDir = 'normal';
            oldYTicks = ax.YTick;
            newYTicks = 2.^oldYTicks;
            newYtickLabs = string(num2str(newYTicks'));
            ax.YTickLabel = newYtickLabs;
            title(sprintf('Ch#%d CWT - Raw',chan))
            xlabel('Time [s]')
            ylabel('Frequency [Hz]')
            c = colorbar;
            c.Label.String = 'CWT coeff. magnitude';
            clear cfs f ax
            
            [cfs,f] = cwt(dataCl(chan,:),'amor',fs,'FrequencyLimits',[1,500]);
            subplot(212)
            imagesc(tAxis,log2(f),abs(cfs))
            axis tight
            ax = gca;
            ax.YDir = 'normal';
            oldYTicks = ax.YTick;
            newYTicks = 2.^oldYTicks;
            newYtickLabs = string(num2str(newYTicks'));
            ax.YTickLabel = newYtickLabs;
            title(sprintf('Ch#%d CWT - %s, %s, %s',chan,decompType,flagType,bssType))
            xlabel('Time [s]')
            ylabel('Frequency [Hz]')
            c = colorbar;
            c.Label.String = 'CWT coeff. magnitude';
            clear cfs f ax
            
            linkaxes(findobj(gcf,'Type','axes'),'x')
        end
    end
end