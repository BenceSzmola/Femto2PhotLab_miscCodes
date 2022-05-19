function [meanTable,stdTable,medianTable] = quickParamExport

[fnames,path] = uigetfile('MultiSelect','on');

for i = 1:length(fnames)
    load([path,fnames{i}],'saveStruct')
    for j = 1:length(saveStruct)
        if isempty(saveStruct(j).imagingEvents.Params.RiseTime2080)
            saveStruct(j).imagingEvents.Params.RiseTime2080 = saveStruct(j).imagingEvents.Params.RiseTime;
        end
        if isempty(saveStruct(j).imagingEvents.Params.DecayTime8020)
            saveStruct(j).imagingEvents.Params.DecayTime8020 = saveStruct(j).imagingEvents.Params.DecayTime;
        end
    end
    iEvs = [saveStruct.imagingEvents];
    iAvg = mean(cell2mat(squeeze(struct2cell([iEvs.Params]))),2);
    iSd = std(cell2mat(squeeze(struct2cell([iEvs.Params]))),[],2);
    iMed = median(cell2mat(squeeze(struct2cell([iEvs.Params]))),2);
    means(i,:) = iAvg;
    stds(i,:) = iSd;
    medians(i,:) = iMed;
    rowNames{i} = ['Block ',num2str(i)];
end

meanTable = array2table(means,'VariableNames',fieldnames([iEvs(1).Params]),'RowNames',rowNames);
stdTable = array2table(stds,'VariableNames',fieldnames([iEvs(1).Params]),'RowNames',rowNames);
medianTable = array2table(medians,'VariableNames',fieldnames([iEvs(1).Params]),'RowNames',rowNames);

[saveFn,saveP] = uiputfile('*.xlsx');
writetable(meanTable,[saveP,saveFn],'Sheet','Means','WriteRowNames',1)
writetable(stdTable,[saveP,saveFn],'Sheet','SD','WriteRowNames',1)
writetable(medianTable,[saveP,saveFn],'Sheet','Medians','WriteRowNames',1)