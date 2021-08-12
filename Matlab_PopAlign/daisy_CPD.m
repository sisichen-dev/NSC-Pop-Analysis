function[newdflist, transforms] = daisy_CPD(dflist, sampsize, niters, dfnames)
% function [newdflist,transforms] = daisy_CPD(dflist, sampsize, niters, )
% This function iteratively determines transformations 

numX = length(dflist);
numT = numX - 1;

% Discover transform parameters: 
for i=1:numT
    df1 = dflist{i} ; % target dataset
    df2 = dflist{i+1} ; % to be transformed dataset

    [newdf2, R, t, ax1] = CPD_correct(df1,df2,sampsize,niters);
    set(gcf,'pos',[ 440   553   560   245]);
    print(['figures/CPD_',dfnames{i+1},'_correct.pdf'],'-dpdf','-r300');
    Rlist{i} = R;
    tlist{i} = t;
end

% Apply transform in daisy chain
newdflist = {}
for i=2:numX
    currX = dflist{i};
    [M,D2]=size(currX);
    
    % get list of transforms to apply:
    Torder = fliplr(1:(i-1));
    
    % Apply transformations backwards: 
    for j=1:length(Torder)
        Tnum = Torder(j);
        currR = Rlist{Tnum};
        currt = tlist{Tnum};
        newX = currX*currR' + ones(M,1)*currt';
        currX = newX;

    end
    newdflist{i} = currX;
end

newdflist{1} = dflist{1};

transforms.Rlist = Rlist;
transforms.tlist = tlist;

%% Temp space: 

% axlist = {};
% totaln = length(newdflist)
% newdfnames = dfnames;
% [r,c]=findsubplotsize(totaln);
% figure;
% i1=1;
% i2=2;
% i3=3;
% svz0pca = newdflist{1};
% for i=1:totaln
%     axlist{i} = subplot(r,c,i); 
%     currpca = newdflist{i};
%     scatter3(currpca(:,i1), currpca(:,i2), currpca(:,i3), 20, 'filled', 'k', 'MarkerFaceAlpha', 0.3); hold on; 
%     %scatter3(svz0pca(:,i1), svz0pca(:,i2), svz0pca(:,i3), 20, 'filled', 'r', 'MarkerFaceAlpha', 0.3); 
%     title(newdfnames{i}, 'interpreter', 'none')
%     xlim(XL);
%     ylim(YL);
%     zlim(ZL);
%    view([52,16]);
%   
% end
% linkprop([axlist{:}],{'Xlim', 'Ylim', 'Zlim', 'CameraPosition'})
% 


