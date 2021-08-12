function makegmms(trainH, validH, name, gmmfolder, rlist, trainIdx, validIdx)
% function makegmms(trainH, validH, name, gmmfolder, rlist, options)
% If sample specific folder exists, exit, otherwise make folder 
% and build gmms. Use cross-validation data to calculate the BIC for the 
% data 

% calculate a regularization value based on the eigenvalue of the
% covariance of the entire matrix:
options = statset('MaxIter',300);
result = cov(trainH');
eigenvalues = eig(result);
regvalue = max(eigenvalues)/25;

for i=1:length(rlist)
    r = rlist(i);
    %%%%%%%%%%%%%%%%%%%%%
    % OPTIONS
    %%%%%%%%%%%%%%%%%%%%%
    
    % [1] directly use EM
    gmfitall = fitgmdist(trainH', r, 'Replicates', 5, 'regularization', regvalue,'Options',options);
    
%     % [2] use kmeans to preseed 
%     S = kmeansSeed(trainH,r);
%     gmfitall = fitgmdist(trainH', r, 'Replicates', 1, 'regularization', regvalue, 'Options', options, 'Start', S);
% 
%     % [3] use kmeans to make model
%    mindet = 0;
%    iteration = 1;
%    
%    while (mindet < 1e-10)
%     disp(['Currently learning iteration  ', num2str(iteration)])
%     S = kmeansSeed(trainH,r);
%     % Calculate covariance matrix determinants
%     dets = [];
%     for i=1:r
%         dets = [dets, det(S.Sigma(:,:,i))];
%     end
%     mindet = min(dets);
%     iteration = iteration + 1;
%     if iteration > 10
%         return
%     end
%    end
%    
%    S = kmeansSeed(trainH,r);
%     gmfitall = gmdistribution(S.mu,S.Sigma,S.ComponentProportion);
%     
    %%%%%%%%%%%%%%%%%%%%%%
    % save gmms
    save([gmmfolder, name,'_', num2str(rlist(i)), '_gmm.mat'], 'gmfitall', 'trainIdx', 'validIdx')
end









