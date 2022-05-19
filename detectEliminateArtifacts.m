function data_denoized = detectEliminateArtifacts(data,fs,method,slideWinSize,autoCorrThr)

if nargin == 0 || isempty(data) || isempty(fs)
    rhdStruct = read_Intan_RHD2000_file_szb;
    data = rhdStruct.amplifier_data;
    goodChans = find([rhdStruct.amplifier_channels.electrode_impedance_magnitude] < 10^6);
    if isempty(goodChans)
        errordlg('all chans'' impedance > 10^6')
        return
    elseif length(goodChans) > 5
        goodChans = goodChans(1:5);
    end
    data = data(goodChans,:);
    fs = rhdStruct.fs;
    assignin('base','data',data)
    assignin('base','fs',fs)
    assignin('base','rhdStruct',rhdStruct)
    assignin('base','goodChans',goodChans)
end

% make sure channels are in rows
if size(data,1) > size(data,2)
    data = data';
end
% save number of channels and samples for easier access downstream
numChans = size(data,1);
numSamples = size(data,2);

% ask user whether they want to run periodic noise filter
answer = questdlg('Do you want to run the periodic noise filter first?','Periodic Filter');
if strcmp(answer,'Yes')
    doPeriodicFilt = true;
else
    doPeriodicFilt = false;
end

% running the periodic noise filter if requested
if doPeriodicFilt
    data = periodicNoise(data,fs);
    assignin('base','perio',data)
end

% set a window size to be used later for going through channels
if nargin < 4
    slideWinSize = 2000;
end

if nargin < 5
    autoCorrThr = 0.999;
end

% set frequency ranges where no information of interest is expected
% badFreqs = badFreqInput(fs);

% compute EEMD 
% % [imfsCell, residsCell, numImfs] = runEEMD(data,fs,20,20);
% % assignin('base','imfsCell',imfsCell)
% % assignin('base','residsCell',residsCell)
% % assignin('base','numImfs',numImfs)
% imfsCell = evalin('base','imfsCell');
% residsCell = evalin('base','residsCell');
% numImfs = evalin('base','numImfs');
% % numImfs = 10;
% imfsCleaned = imfsCell;
% 
% % determine up to which IMF level you want to use
% upToImfNum = min(min(numImfs),20);

% compute DWT decomposition
% wname = 'db4';
% freqLow = 4;
% decLvl = round(log2(fs/freqLow));
% for i = 1:numChans
%     [c(i,:),l(i,:)] = wavedec(data(i,:),decLvl,wname);
%     details(i,:) = detcoef(c(i,:),l(i,:),1:decLvl);
%     approx(i,:) = appcoef(c(i,:),l(i,:),wname);
% end
% detailsClean = details;
% approxClean = approx;

% compute SWT decomposition
wname = 'db4';
freqLow = 4;
decLvl = round(log2(fs/freqLow));
padLen = numSamples/(2^decLvl);
padLen = (2^decLvl)*round(padLen) - numSamples;
dataPad = [data, zeros(numChans,padLen)];
swc = cell(numChans,1);
for chan = 1:numChans
    swc{chan} = swt(dataPad(chan,:),decLvl,wname);
end
swcClean = swc;

if nargin < 3
    method = 7;
end
switch method
    case 1 % simple version just based on IMF autocorrelations, when a IMF has higher autocorrelation than a
           % given threshold then its just discarded completely
           
        for chan = 1:numChans
            for imfNum = 1:numImfs(chan)
%                 for i = 1:slideWinSize:numSamples
%                     if (i + slideWinSize -1) > numSamples
%                         win = i:numSamples;
%                     else
%                         win = i:(i + slideWinSize - 1);
%                     end
%                     
                    r = xcorr(imfsCell{chan}(imfNum,:),3,'coeff');
                    if r(end) < 0.999
                        imfsCleaned{chan}(imfNum,:) = 0;
                    end
%                 end
            end
        end
        reconstr = cell2mat(cellfun(@sum, imfsCleaned, 'UniformOutput',false));
        reconstr = reconstr + cell2mat(residsCell);
        
    case 2 % throw away IMF segments based on IMF-wise IEC
        
        IECs = IMFwiseIEC(imfsCell,upToImfNum,numSamples,slideWinSize);
        figure;
        imagesc(IECs)
        title('IMF wise IECs')
        highIECwins = IECs > 0.85;
        
         % variable to store the high IEC segments for each IMF level
        segments = extractSegments(highIECwins,0);
        for imfNum = 1:upToImfNum

            currIMFs = cell2mat(cellfun(@(x)x(imfNum,:), imfsCell, 'UniformOutput',false));
            currIMFs_cleaned = currIMFs;
            
            for i = 1:size(segments{imfNum},1)
                currSegInds = segments{imfNum}(i,1):segments{imfNum}(i,2);
                currIMFs_cleaned(:,currSegInds) = 0;
            end
            
            % storing the IMFs containing cleaned flagged segments
            for i = 1:numChans
                imfsCleaned{i}(imfNum,:) = currIMFs_cleaned(i,:);
            end
            
        end
        
        reconstr = cell2mat(cellfun(@sum, imfsCleaned, 'UniformOutput', false));
        reconstr = reconstr + cell2mat(residsCell);
        
    case 3 % IMF thresholding
        
        % compute IMF-wise inter-electrode-correlation
        IECs = IMFwiseIEC(imfsCell,upToImfNum,numSamples,slideWinSize);
        % finding high IEC windows
        highIECwins = IECs > (mean(IECs,2) + 2*std(IECs,0,2));
        
        for chan = 1:numChans
            for imfNum = 1:upToImfNum
                % get the peaks of the IMF
                [~,locs] = findpeaks(abs(imfsCell{chan}(imfNum,:)));
                % check which peaks are in a high IEC region
                highIECpeakInds = ismember(locs,find(highIECwins(imfNum,:)));
                if isempty(highIECpeakInds)
                    display('sad')
                    continue
                end
                avg = mean(abs(imfsCleaned{chan}(imfNum,highIECpeakInds)));
                sd = std(abs(imfsCleaned{chan}(imfNum,highIECpeakInds)));
                imfsCleaned{chan}(imfNum,highIECpeakInds) = (abs(imfsCleaned{chan}(imfNum,highIECpeakInds)) - (avg + 2*sd))...
                    .* sign(imfsCleaned{chan}(imfNum,highIECpeakInds));
            end
        end
        
        reconstr = cell2mat(cellfun(@sum, imfsCleaned, 'UniformOutput',false));
        reconstr = reconstr + cell2mat(residsCell);
        
    case 4 % check IEC and threshold on raw time domain
        
        IECs = zeros(numChans,numSamples);
        for i = 1:slideWinSize:numSamples
            win = i:i+slideWinSize-1;
            if win(end) > numSamples
                win = i:numSamples;
            end
            IECs(:,win) = mean(corrcoef(data(:, win)'), 'all');
        end
        
        slideAvg = movmean(abs(data),slideWinSize,2);
        slideSD = movstd(abs(data),slideWinSize,0,2);
        
        suspInds = (IECs > 0.8) & (abs(data) > (slideAvg + 1*slideSD));
        
        reconstr = abs(data);
        reconstr(suspInds) = reconstr(suspInds) - (slideAvg(suspInds) + 2*slideSD(suspInds));
        
        reconstr = reconstr.*sign(data);
    case 5 
        
        % compute IMF-wise inter-electrode-correlation
        IECs = IMFwiseIEC(imfsCell,upToImfNum,numSamples,slideWinSize);

        % finding high IEC windows
%         highIECwins = IECs > (mean(IECs,2) + 2*std(IECs,0,2));
        highIECwins = IECs > 0.8;
        figure;
        imagesc(highIECwins)
        dewIT = msgbox('go forth?');
        waitfor(dewIT)

        % variable to store the high IEC segments for each IMF level
        segments = extractSegments(highIECwins,slideWinSize);
        flagSegsFig = figure('WindowState','maximized');
        for imfNum = 1:upToImfNum

            % run ICA or CCA on the flagged segments
            currIMFs = cell2mat(cellfun(@(x)x(imfNum,:), imfsCell, 'UniformOutput',false));
            currIMFs_cleaned = currIMFs;
            
            ICAorCCA = 1;
            for i = 1:size(segments{imfNum},1)
                currSegInds = segments{imfNum}(i,1):segments{imfNum}(i,2);
                
                switch ICAorCCA
                    case 1
                        % computing ICA and then asking the user to choose which ICs they want to discard
                        [badSegICs,A,W] = fastica(currIMFs(:,currSegInds));
                        flagSegsFig.Name = sprintf('IMF#%d, flagged segment#%d, ICA',imfNum,i);
                        for j = 1:size(badSegICs,1)
                            subplot(size(badSegICs,1),1,j)
                            plot(badSegICs(j,:))
                        end
                        [IC2discard,tf] = listdlg('ListString',num2str((1:size(badSegICs,1))'));
                        if ~tf
                            errordlg('Algorithm stopped')
                            return
                        end
                        badSegICs(IC2discard,:) = 0;
                        reconstr = A*badSegICs;
                    case 2
                        % compute CCA, discard sources with the autocorrelations lower than threshold
                        reconstr = doCCAdiscard(currIMFs(:,currSegInds),autoCorrThr);
                end
                currIMFs_cleaned(:,currSegInds) = reconstr;
            end

            % storing the IMFs containing cleaned flagged segments
            for i = 1:numChans
                imfsCleaned{i}(imfNum,:) = currIMFs_cleaned(i,:);
            end
        end
        close(flagSegsFig)

        % reconstructing channels from the cleaned IMFs (residuals not included
        % yet!!!)
        reconstr = cell2mat(cellfun(@sum, imfsCleaned, 'UniformOutput', false));
        reconstr = reconstr + cell2mat(residsCell);
        
    case 6 % check autocorrelation of raw signal in windows, flag suspicious windows, decompose them and get
        % rid of bad IMFs
        
        slideAutoCorr = zeros(size(data));
        reconstr = data;
        for chan = 1:numChans 
            for i = 1:slideWinSize:numSamples
                if (i + slideWinSize - 1) > numSamples
                    win = i:numSamples;
                else
                    win = i:(i + slideWinSize - 1);
                end
                r = xcorr(data(chan,win),10,'coeff');
                slideAutoCorr(chan,win) = r(end);
                if r(end) < 0.9
                    temp = doCCAdiscard(imfsCell{chan}(:,win),0.99);
                    reconstr(chan,win) = sum(temp) + residsCell{chan}(win);
                end
            end
            figure;
            yyaxis left
            plot(data(chan,:))
            yyaxis right
            plot(slideAutoCorr(chan,:))
            title(sprintf('Chan#%d - autocorrs',chan))
        end
        
    case 7 % check IEC on dwt components, here the identity of the component levels is assured across channels
        % unlike in non multivariate EMD
        
        % computing dwt component wise IEC
        IECs = dwtIEC(details,approx,numSamples,slideWinSize);
        
        % looping through dwt components, indices which had high IEC get set to 0
        for dwtLvl = 1:decLvl+1
            highIECinds = IECs{dwtLvl} > 0.75;
            
            if dwtLvl <= decLvl
                temp = cell2mat(detailsClean(:,dwtLvl));
                temp(:,highIECinds) = 0;
                temp = mat2cell(temp,ones(1,numChans));
                detailsClean(:,dwtLvl) = temp;
            else
                approxClean(:,highIECinds) = 0;
            end
        end
        
        % reconstructing vector c for dwt function
        for chan = 1:numChans
            cClean = zeros(size(c(chan,:)));
            for i = 1:length(l(chan,:))-1
                if i == 1
                    cClean(i:l(chan,1)) = approxClean(chan,:);
                else
                    cClean(sum(l(chan,1:i-1))+1:sum(l(chan,1:i))) = detailsClean{chan,length(l(chan,:))-i};
                end
            end
            reconstr(chan,:) = waverec(cClean,l(chan,:),wname);
        end
        
    case 8 % using SWT discard segments with high component-wise IEC
        IECs = swtIEC(swc,numSamples+padLen,slideWinSize);
        highIECinds = IECs > 0.9;
        figure;
        imagesc(highIECinds)
        title('high SWT IEC inds')
        
        reconstr = zeros(size(data));
        for chan = 1:numChans 
            temp = swcClean{chan};
            temp(highIECinds) = 0;
            
            reconstrPad = iswt(temp,wname);
            reconstr(chan,:) = reconstrPad(1:numSamples);
        end
        
        
end

figure('Name','Comparing OG to cleaned','WindowState','maximized');
for i = 1:numChans
    subplot(211)
    plot(data(i,:))
    subplot(212)
    plot(reconstr(i,:))
%     hold on
%     transfSuspInds = 1:numSamples;
%     transfSuspInds = transfSuspInds(suspInds(i,:));
%     plot(transfSuspInds,zeros(size(transfSuspInds)),'r*')
%     hold off
    linkaxes(findobj(gcf,'Type','axes'),'xy')
    next = msgbox('Press ok to go to next channel');
    waitfor(next)
end

% assign output value
if nargout > 0
    data_denoized = reconstr;
end

% % start looping over the channels
% for chan = 1:numChans
%     
%     % step through channel in computationally efficient windows and check
%     % characteristics which can help decide whether the window is
%     % contaminated with artifacts or not
%     for winStart = 1:winSize:numSamples
%         
% %         % compute the CWT in the current window
% %         [cwtCoeffs, cwtFreqs] = cwt(data(chan, winStart:(winStart+winSize-1)), 'amor', fs,...
% %             'FrequencyLimits',[min(badFreqs,[],'all'), max(badFreqs,[],'all')]);
% %         
% %         badCwtRows = false(size(cwtFreqs));
% %         for i = 1:size(badFreqs,1)
% %             badCwtRows(cwtFreqs <= badFreqs(i,1) & cwtFreqs >= badFreqs(i,2)) = true;
% %         end
%         
%         
%     end
%     
% end
% 
% end
end

function reconstr = doCCAdiscard(inputMat,discardThr)
% compute CCA, discard sources with the autocorrelations lower than threshold
xDelay = [zeros(size(inputMat,1),3), inputMat(:,1:end-3)];
[A,~,r,U,~,~] = canoncorr(inputMat',xDelay');
sources2discard = r < discardThr;

Amod = A;
Amod(:,sources2discard) = 0;
reconstr = (U*pinv(Amod))';

end

function segments = extractSegments(inputMat,minSegLen)
    
    segments = cell(size(inputMat,1),1);
    
    for i = 1:size(inputMat,1)
        % find indexes above the threshold for a given amount of time
        if ~iscell(inputMat)
            aboveThrInds = find(inputMat(i,:));
        else
            aboveThrInds = find(inputMat{i});
        end
        
        if isempty(aboveThrInds)
            continue
        end
        indDiffs = diff(aboveThrInds);
        disconts = find(indDiffs ~= 1);

        currSegments = [aboveThrInds([1, disconts+1])', aboveThrInds([disconts, length(aboveThrInds)])'];

        currSegments(currSegments(:,2)-currSegments(:,1) < minSegLen,:) = [];

        segments{i} = currSegments;
    end
end

function IECs = swtIEC(swcCell,numSamplesPad,slideWinSize)
    numComps = size(swcCell{1},1);
    % create IEC placeholder
    IECs = zeros(numComps,numSamplesPad);
    
    for compNum = 1:numComps
        currComps = cell2mat(cellfun(@(x) x(compNum,:), swcCell,'UniformOutput',false));
        for i = 1:slideWinSize:numSamplesPad
            if (i + slideWinSize -1) > numSamplesPad
                win = i:numSamplesPad;
            else
                win = i:(i + slideWinSize -1);
            end
            corrVals = corrcoef(currComps(:,win)');
            corrVals = tril(corrVals,-1);
            corrVals(corrVals == 0) = [];
            IECs(compNum,win) = mean(corrVals);
        end
    end
end

function IECs = dwtIEC(detailsCell,approx,numSamples,slideWinSize)
    % create the IECs placeholder, here a cell is needed because of the unequal lengths of components
    IECs = cell(size(detailsCell,2)+1,1);
    
    for detCoefNum = 1:(size(detailsCell,2)+1)
        if detCoefNum <= size(detailsCell,2)
            currComps = cell2mat(detailsCell(:,detCoefNum));
        else
            currComps = approx;
        end
        currCompsLen = size(currComps,2);
        slideWinSizeComp = round(slideWinSize / (2^detCoefNum));
        for i = 1:slideWinSizeComp:currCompsLen
            if (i + slideWinSizeComp -1) > currCompsLen
                win = i:currCompsLen;
            else
                win = i:(i + slideWinSizeComp -1);
            end
            corrVals = corrcoef(currComps(:,win)');
            corrVals = tril(corrVals,-1);
            corrVals(corrVals == 0) = [];
            IECs{detCoefNum}(win) = mean(corrVals);
        end
    end
end

function IECs = IMFwiseIEC(imfsCell,upToImfNum,numSamples,slideWinSize)
    % compute IMF-wise inter-electrode-correlation
    IECs = zeros(upToImfNum,numSamples);
    for imfNum = 1:upToImfNum
        currIMFs = cell2mat(cellfun(@(x)x(imfNum,:), imfsCell, 'UniformOutput',false));
        for i = 1:slideWinSize:numSamples
            win = i:i+slideWinSize-1;
            if win(end) > numSamples
                win = i:numSamples;
            end
            corrVals = corrcoef(currIMFs(:,win)');
            corrVals = tril(corrVals,-1);
            corrVals(corrVals == 0) = [];
            IECs(imfNum, win) = mean(corrVals);
        end
    end
end

function badFreqs = badFreqInput(fs)
    inputBadFreqsFig = figure('NumberTitle','off',...
        'Name','Frequency ranges of no interest',...
        'CloseRequestFcn',@(h,e) saveRanges);
    
    colNames = {'Start [Hz]','End [Hz]'};
    badFreqTable = uitable(inputBadFreqsFig,...
        'Units','normalized',...
        'Position',[0.01, 0.1, 0.9, 0.8],...
        'ColumnName',colNames,...
        'ColumnEditable',true,...
        'CellEditCallback',@ freqInputControll);
    
    uicontrol(inputBadFreqsFig,...
        'Style','pushbutton',...
        'Units','normalized',...
        'Position',[0.01, 0.01, 0.1, 0.05],...
        'String','Add range',...
        'Callback',@(h,e) addRow);
    
    uiwait(inputBadFreqsFig)
    
    function freqInputControll(~,event)
        ind1 = event.Indices(1);
        ind2 = event.Indices(2);
        if event.NewData > (fs/2)
            badFreqTable.Data(ind1,ind2) = fs/2;
        elseif event.NewData < 1
            badFreqTable.Data(ind1,ind2) = 1;
        end
    end
    
    function addRow
        badFreqTable.Data = [badFreqTable.Data; [0, 0]];
    end

    function saveRanges
        if ~isempty(find(badFreqTable.Data(:,1) >= badFreqTable.Data(:,2), 1))
            errordlg('Start shouldnt be higher then end frequency!')
            return
        end
        
        badFreqs = badFreqTable.Data;
        delete(inputBadFreqsFig)
    end
end