function [gmmerravg, dataerravg] = calcGMMerr(gmmdata, Hdata, sampleInd, sampsize, binwidths, plotBool)
% function calcGMMerr(gmmdata, Hdata, sampleInd)
% This function calculates the GMM error as well as the data resampling
% error. Error is P(orig) - P(gmm) or P(resample) across all bins of the
% probability space. binsize denotes the binsize along each axis of the
% feature space. 
% 
% good defaults to use: 
% sampsize = 2000;
% binsize = 5;

maroon = [162,20,47]/255;

gmmlist = gmmdata.gmmlist;%[{gmmdata.refmodel}, gmmdata.testlist];
Hlist = Hdata.Hlist;%;[Hdata.reflist, Hdata.testlist];
nameslist = Hdata.dfnames;%[Hdata.refnames, Hdata.testnames];
currgmm = gmmlist{sampleInd};
currH = Hlist{sampleInd};
currname = nameslist{sampleInd};

n = size(currH, 2);
M = random(currgmm, n)';

errv = calc2DErr(M,currH, binwidths);

%% recalculate the error using random sampling from the original data: 

if (sampsize < (n/2))
    inds = randsample(n,sampsize);
    temp = setdiff(1:n,inds);
    altinds = randsample(temp, sampsize);
else
    inds = sort(randsample(1:n,n/2))';
    altinds = setdiff(1:n,inds);
    if length(altinds)>length(inds)
        altinds = altinds(1:end-1);
    end
end

H1 = currH(:,inds);
H2 = currH(:,altinds);

dataerrv = calc2DErr(H1, H2, binwidths);

%% Plot error together

gmmerravg = mean(errv);
dataerravg = mean(dataerrv);

if plotBool
    
    figure; plot(dataerrv, 'k-','LineWidth',1.5); hold on;
    plot(errv, '-', 'LineWidth', 1.5, 'Color', maroon); hold on;
    line([0, length(dataerrv)],[gmmerravg, gmmerravg], 'LineWidth', 4,'Color', maroon);
    line([0, length(dataerrv)],[dataerravg, dataerravg], 'LineWidth', 4, 'Color', 'k');
    set(gca,'FontSize',18);
    set(gca, 'xtick',1:length(dataerrv),'xticklabel',[]);
    xlim([0,length(dataerrv)])
    ylabel('% error')
    xlabel('2d projection index')
    legend({'data','model'})

    set(gcf, 'Position', [0,0,300,400])
end

