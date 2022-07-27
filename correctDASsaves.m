path = uigetdir(cd, 'Select directory to correct!');
path = [path,'\'];
saveFnames = dir([path,'DASsave*.mat']);
saveFnames = {saveFnames.name};

for i = 1:length(saveFnames)
    load([path, saveFnames{i}], 'ephysSaveData', 'ephysSaveInfo', 'imagingSaveData', 'imagingSaveInfo', 'runData',...
        'simultSaveData', 'simultSaveInfo', 'comments')
    
    if isempty(ephysSaveData) || isempty(ephysSaveInfo)
        continue
    end
    
    if isfield(ephysSaveInfo.DetSettings, 'RefCh')
        refch = ephysSaveInfo.DetSettings.RefCh;
        if any(ismember(refch, ephysSaveInfo.DetChannel))
            row2del = ismember(ephysSaveInfo.DetChannel, refch);
            ephysSaveInfo.DetChannel(row2del) = [];
            ephysSaveInfo.AllChannel(ismember(ephysSaveInfo.AllChannel, refch)) = [];
            ephysSaveData.RawData(ismember(ephysSaveInfo.AllChannel, refch),:) = [];
            ephysSaveData.Dets(row2del) = [];
            ephysSaveData.DetBorders(row2del) = [];
            ephysSaveData.DetParams(row2del) = [];
            ephysSaveData.EventComplexes(row2del) = [];
            
        end
    else
        continue
    end
    
    save([path, saveFnames{i}], 'ephysSaveData', 'ephysSaveInfo', 'imagingSaveData', 'imagingSaveInfo', 'runData',...
        'simultSaveData', 'simultSaveInfo', 'comments')
end

clear all