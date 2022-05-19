% testEEMDscript
clear imfs4ICA imfICs numICs

oneOverFun = @(theta,xdata) theta(1)./xdata;

for i = 1:length(imfs)
    imfs4ICA = [imfs{i}; resids{i}];
    [imfICs{i},A{i},W{i}] = fastica(imfs4ICA,'maxNumIterations',2000);
    numICs = size(imfICs{i},1);
    
    fftFitParam{i} = zeros(numICs,1);
    [faxis,icaPSDs] = freqspec(imfICs{i},fs,0);

    icsFig = figure('Name',['Independent components'],'WindowState','maximized');
    ICs2discard{i} = false(numICs,1);
    for j = 1:numICs
        subplot(211)
        plot(imfICs{i}(j,:))
        subplot(212)
        plot(faxis(j,1:find(faxis(j,:) > 1000,1)),icaPSDs(j,1:find(faxis(j,:) > 1000,1)))
        discardIC = questdlg('Discard IC?','IC checking');
        if strcmp(discardIC,'Yes')
            ICs2discard{i}(j) = true;
        end
    end
    delete(icsFig)
    
%     figure('Name',['EEMD-ICA - ch#',num2str(i),' - ICs from IMFs']);
%     subplot(round(numICs/2)+1,2,1)
%     plot(data(i,:))
%     title('Raw data')
%     
%     for j = 1:numICs
%         faxis1000 = find(faxis(j,:) > 1000, 1);
%         fftFitParam{i}(j) = lsqcurvefit(oneOverFun,0.9,faxis(j,2:faxis1000),icaPSDs(j,2:faxis1000));
%         subplot(ceil(numICs/2)+1,2,j+1)
%         plot(imfICs{i}(j,:))
%         title(sprintf('IC #%d - fit param = %.3f',j,fftFitParam{i}(j)))
%     end
%     linkaxes(findobj(gcf,'Type','axes'),'x')

    modA = A{i};
    modA(:,ICs2discard{i}) = 0;
    reconstrImfs = modA*imfICs{i};
    
%     reconstrImfsFig = figure('Name',sprintf('Reconstructed Imfs, chan #%d',i),'WindowState','maximized');
%     for j = 1:size(imfs4ICA,1)
%         subplot(211)
%         plot(imfs4ICA(j,:))
%         title('OG IMF')
%         subplot(212)
%         plot(reconstrImfs(j,:))
%         title('Reconstr')
%         waitforbuttonpress
%     end
%     delete(reconstrImfsFig)

    emdIcaReconstr(i,:) = sum(reconstrImfs);
    
    figure('Name',['EMD-ICA - ch#',num2str(i)]);
    subplot(211)
    plot(data(i,:))
    title(['Raw channel #',num2str(i)])
    subplot(212)
    plot(emdIcaReconstr(i,:))
    title('Reconstructed after EMD-ICA denoising')
    linkaxes(findobj(gcf,'Type','axes'),'x')
end