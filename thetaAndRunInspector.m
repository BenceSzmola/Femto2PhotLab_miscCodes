rhdStruct = read_Intan_RHD2000_file_szb;

data = rhdStruct.amplifier_data;
fs = rhdStruct.fs;
tAxis = rhdStruct.tAxis;
goodChs = find([rhdStruct.amplifier_channels.electrode_impedance_magnitude] < 10^6)
data = data(goodChs(1:5),:);

clearvars -except data fs tAxis


fsDs = 1250;
[dataDs,tAxisDs] = downSampData(fs,fsDs,tAxis,data);

data_cl = periodicNoise(dataDs,[],fsDs,300,27.78);
data_cl = periodicNoise(data_cl,[],fsDs,300,150);

[csvName, path] = uigetfile('*.csv','Select run csv file!');
gramoData = csvread([path,csvName],1,0);

selFig = figure;
plot(tAxisDs,dataDs(1,:))
inds2use = graphicIntervalSelector(selFig,tAxisDs);

for i = 1:5
    [coeffs(:,:,i),f] = cwt(data_cl(i,:), 'amor', fsDs, 'FrequencyLimits', [1, 300]);
    coeffs = abs(coeffs).^2;
    instE(i,:) = trapz(squeeze(coeffs(:,:,i)));
    thetaFreqs = find(round(f,2)>6 & round(f,2)<12);
    thetaInstE(i,:) = trapz(squeeze(coeffs(thetaFreqs,:,i)));
    thetaFract(i,:) = thetaInstE(i,:) ./ instE(i,:);
    thetaInstE(i,setxor(1:length(tAxisDs), inds2use)) = nan;
    thetaFract(i,setxor(1:length(tAxisDs), inds2use)) = nan;
end

figure;
for i = 1:5
    subplot(5,1,i)
    plot(tAxisDs,thetaFract(i,:))
end
subplot(515)
plot(gramoData(:,1)/1000,gramoData(:,2))
linkaxes(findobj(gcf,'Type','axes'),'x')
sgtitle('Theta band (6-12 Hz) instantaneous energy ratio to 1-300 Hz instantaneous energy')

figure;
for i = 1:5
    subplot(5,1,i)
    plot(tAxisDs,thetaInstE(i,:))
end
subplot(515)
plot(gramoData(:,1)/1000,gramoData(:,2))
linkaxes(findobj(gcf,'Type','axes'),'x')
sgtitle('Theta band (6-12 Hz) instantaneous energy')