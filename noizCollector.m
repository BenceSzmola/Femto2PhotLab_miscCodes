if ~exist('fn','var')
    fn = 'no prev file';
end

[fn,path] = uigetfile('*.rhd',['Prev: ',fn]);
if ~fn
    return
end

rhdStruct = read_Intan_RHD2000_file_szb([path,fn]);
data = rhdStruct.amplifier_data;
numCh = size(data,1);

dataFig = figure('NumberTitle','off','Name',fn,'WindowState','maximized');
for i = 1:numCh
    subplot(numCh,1,i)
    plot(data(i,:))
    title(sprintf('Ch#%d',i))
    linkaxes(findobj(dataFig,'Type','axes'),'x')
end

inpGud = false;
while ~inpGud
%     opts.WindowStyle = 'normal';
%     answer = inputdlg('Specifiy noise segments! (startInd1,stopInd1,...)','Bords',[1 50],{''},opts);
%     if isempty(answer)
%         return
%     end
%     
%     answer = str2num(answer{:});

    [answer,~] = ginput;
    answer = round(answer);
    answer(answer <= 0) = 1;
    answer(answer > size(data,2)) = size(data,2);

    if mod(length(answer),2) == 0
        bordInds = zeros(length(answer)/2,2);
        inpGud = true;
    else
        eD = errordlg('Give pairs of indexes (start of segment, end of segment)');
        waitfor(eD)
        clear eD
    end
end

if ~isempty(answer)
    bordInds(:,1) = answer(1:2:end);
    bordInds(:,2) = answer(2:2:end);

    temp = cell(size(bordInds,1),1);
    for i = 1:size(bordInds,1)
        temp{i} = data(:,bordInds(i,1):bordInds(i,2));
    end

    if exist('spwChunks','var')
        spwChunks = [spwChunks; temp];
    else
        spwChunks = temp;
    end
    clear temp
end

if exist('dataFig','var')
    close(dataFig)
end

cont = questdlg(['Load next RHD? (prev was ',fn,')']);
if strcmp(cont,'Yes')
    noizCollector
end
clearvars -except spwChunks fn