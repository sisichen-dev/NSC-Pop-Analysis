
function Hdata = checkGMMs(gmmdata, Hdata, logdata,feats)
% function checkGMMs(gmmdata, Hdata, logdata,feats)
%  This function checks gmms by displaying plots showing their quality: 
% 1. tsne colored by subpouplation
% 2. pairwise distance distributions by subpopulation
% 3. Hmatrix clustered
% 4. gene expression heatmap

numtoplot = length(gmmdata.bestlist);

for i = 3 %1:numtoplot
    modelnums = gmmdata.bestlist; % 1:length(gmmdata.sampnums);
    modelind = modelnums(i);
    sampind = gmmdata.sampnums(i);
    
    currdfname = Hdata.dfnames{sampind};
    currH = Hdata.Hlist{sampind};
    clusters = gmmdata.clusterlist{modelind};
    Y = Hdata.tsnelist{sampind}; % get tsne data    

    % make names for figures;
    fname1= ['figures/check_tsne_',num2str(sampind),'_', currdfname, '_model-', num2str(modelind), '.pdf'];
    fname2= ['figures/check_pdist_',num2str(sampind),'_', currdfname, '_model-', num2str(modelind), '.pdf'];
    fname3= ['figures/check_Hmat_',num2str(sampind),'_', currdfname, '_model-', num2str(modelind),'.pdf'];
    fname4= ['figures/check_gex_',num2str(sampind),'_', currdfname, '_model-', num2str(modelind),'.pdf'];
    allname = ['figures/check_',num2str(sampind),'_', currdfname, '_model-', num2str(modelind),'.pdf'];
    
    plottsne(clusters,Y,currdfname)
    print(fname1, '-dpdf', '-r300')
    
    plotpdist(clusters,currH)
    print(fname2, '-dpdf', '-r300')
    
    plotHmat(clusters,currH,feats)
    print(fname3, '-dpdf', '-r300')
    
    plotGMMheatmaplow(gmmdata.clusterlist{sampind}, logdata.logdflist{sampind}, logdata.genes, false)
    print(fname4, '-dpdf', '-r300')
    
    append_pdfs(allname, fname1, fname2, fname3,fname4);
    delete(fname1)
    delete(fname2)
    delete(fname3)
    delete(fname4)
end

close all

function plottsne(clusters, Y, currdfname)
% function plottsne(clusters, Y)
% plot tsne

figure; 
colormap(brewermap(length(unique(clusters)),'Set2'))
scatter(Y(:,1), Y(:,2), 20, clusters,'filled');
set(gca,'XTick', [], 'YTick', []);
title(currdfname)
xlabel('t-sne 1')
ylabel('t-sne 2')
colorbar


    