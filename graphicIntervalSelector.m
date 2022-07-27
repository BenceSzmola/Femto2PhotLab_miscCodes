function inds2use = graphicIntervalSelector(selFig,tAxis)


inpGud = false;
figure(selFig)
while ~inpGud
    try
        [selPoints, ~] = ginput;
    catch
        selPoints = [];
    end
    selPoints(selPoints <= tAxis(1)) = tAxis(1);
    selPoints(selPoints > tAxis(end)) = tAxis(end);

    if mod(length(selPoints), 2) == 0
        inpGud = true;
    else
        eD = errordlg('Repeat selection and give pairs of indexes (start of segment, end of segment)!');
        pause(1)
        if ishandle(eD)
            close(eD)
        end
    end
end


if (length(selPoints) / 2) == 1
    msgTxt = '%d interval selected!';
elseif (length(selPoints) / 2) == 0
    msgTxt = 'No intervals selected!';
else
    msgTxt = '%d intervals selected!';
end
mb = msgbox(sprintf(msgTxt, length(selPoints) / 2));
pause(1)
if ishandle(mb)
    close(mb)
end

inds2use = 1:length(tAxis);
for i = 1:2:length(selPoints)
    [~, intervalStart] = min(abs(selPoints(i) - tAxis));
    [~, intervalEnd] = min(abs(selPoints(i+1) - tAxis));
    inds2use(intervalStart:intervalEnd) = nan;
end
inds2use(isnan(inds2use)) = [];