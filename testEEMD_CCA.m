testEEMDscript
clear A B r U V stats imfs4cca imfs4ccaDelay ccaReconstrIMFs ccaIMFreconstr

for i = 1:length(imfs)
    imfs4cca = [imfs{i}; resids{i}];
    imfs4ccaDelay = imfs4cca(:,2:end);
    imfs4ccaDelay = [imfs4ccaDelay,zeros(size(imfs4ccaDelay,1),1)];
    [A{i},B{i},r{i},U{i},V{i},stats{i}] = canoncorr(imfs4cca',imfs4ccaDelay');
    
    source2discard{i} = false(size(U{i},2),1);
%     [faxis,ccaSourcePSDs] = freqspec(U{i},fs,0,0,1000);
%     wgn = 0.2 * std(U{i}) .* randn(size(U{i}));
%     [wgn_faxis,wgn_psd] = freqspec(wgn,fs,0,0,1000);
%     ccaSourceFig = figure('Name',['CCA sources'],'WindowState','maximized');
%     for j = 1:size(U{i},2)
%         subplot(211)
%         plot(U{i}(:,j))
%         subplot(212)
%         plot(faxis(j,1:find(faxis(j,:) > 1000,1)),ccaSourcePSDs(j,1:find(faxis(j,:) > 1000,1)))
%         sim2wgn = corr(ccaSourcePSDs(j,:)',wgn_psd(j,:)');
%         title(sprintf('Correlation with WGN: %.2f',sim2wgn))
% %         delORno = questdlg('Delete this CCA source?');
% %         if strcmp(delORno,'Yes')
% %             source2discard{i}(j) = true;
% %         end
%         waitforbuttonpress
%     end
%     delete(ccaSourceFig)
    lowAutoCorrSources = (1 - r{i}) > 0.001;
    source2discard{i}(lowAutoCorrSources) = true;
    
    Amod = A{i};
    Amod(:,source2discard{i}) = 0;
    ccaReconstrIMFs = (U{i}*pinv(Amod))';
    ccaIMFreconstr(i,:) = sum(ccaReconstrIMFs);
    
    figure('Name',['EEMD-CCA - ch#',num2str(i)]);
    subplot(211)
    plot(data(i,:))
    title(['Raw channel #',num2str(i)])
    subplot(212)
    plot(ccaIMFreconstr(i,:))
    title(sprintf('Reconstructed after EEMD-CCA denoising, last %d discarded',length(find(lowAutoCorrSources))))
    linkaxes(findobj(gcf,'Type','axes'),'x')
end