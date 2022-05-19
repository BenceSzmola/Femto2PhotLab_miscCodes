
inSamps = find(rhdStruct.tAxis > 12, 1) : find(rhdStruct.tAxis > 14, 1);

figure;
for i = inSamps(1):1000:(inSamps(end)-2000)
corrMat = corr(data(:,i:i+2000)');
t = num2cell(corrMat); % extact values into cells
t = cellfun(@num2str, t, 'UniformOutput', false);
x = repmat(1:size(corrMat,1),size(corrMat,1),1); % generate x-coordinates
y = x'; % generate y-coordinates
subplot(211)
plot(rhdStruct.tAxis(i:i+2000), data(1,i:i+2000))
subplot(212)
imagesc(corrMat)
text(x(:), y(:), t, 'HorizontalAlignment', 'Center')
colorbar
waitforbuttonpress
end