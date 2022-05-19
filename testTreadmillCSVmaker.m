function testTreadmillCSVmaker(fname,tStart,tEnd,movementStartStops)

varNames = {'time','velocity','AbsolutePosition','Lap','RelativePosition','Lick'};
fs = 1/0.005;
duration = tEnd - tStart;
table2write = table('Size',[duration*fs,6],'VariableTypes',repmat("double",1,6),'VariableNames',varNames);
assignin('base','tab2wr',table2write)

% creating time variable, it is in ms
table2write.time = linspace(tStart*1000,tEnd*1000,duration*fs)';

movementInds = cell(size(movementStartStops,1),1);
for i = 1:size(movementStartStops,1)
    movementInds{i} = round(movementStartStops(i,1)*fs):round(movementStartStops(i,2)*fs);
    table2write.velocity(movementInds{i}) = 5 + (10-5)*rand(length(movementInds{i}),1);
end

writetable(table2write,[fname,'.csv'])

end