function feats = makefeat(filtdata, klist, sampsize, featfolder)
% function feats = makefeat(filtdata, klist, sampsize,featfolder)
% Uses oNMF to define featuresets using data sampled from all datasets. To
% choose the best featureset, we compute oNMF features across a range of k
% and then pick the featureset that minimizes the cost function
% (f(residual,k)). oNMF is compute-intensive so we use parfor loops to run
% oNMF across many k at once. 
% filtdata = data structure output from filterdata
% klist = list of the number 
 

fprintf('\n')
disp('************************************')
disp('Making features....')
disp('************************************')
fprintf('\n')

% sample the data
filtdatalist = filtdata.filtdflist;
featfigfolder = [featfolder,'figs/'];
expname = filtdata.expname;
featnames = strcat(featfolder, expname,'_W',strsplit(num2str(klist)), '_features.mat');
nvalid = 100; % number of samples for validation data

opts=struct();
opts.orthogonal=[1,0];
opts.iter = 200;

klistforNMF = klist;
tally =[];

% Check if filtfolder already exists, if not make it.
if (exist(featfolder,'dir')==7)
    
    if (exist(featfigfolder,'dir')==7)
        % tally folders that already exist
        for i=1:length(klist)
            if (exist(featnames{i},'file'))==2
                tally = [tally, i];
            end
        end
    else
        mkdir(pwd, featfigfolder)
    end
    
else
    mkdir(pwd,featfolder)
    mkdir(pwd, featfigfolder)
end

% remove the already computed k from the NMF list: 
if length(tally)>0
    klistforNMF(tally)=[];
end

% set parameters for parfor
myCluster = parcluster('local');
totworkers = myCluster.NumWorkers;
useworkers = totworkers - 1;
initVal = 1;
endVal = length(klistforNMF);

% Sample the imput data files into a sample df for orthogonal nmf
rng(1)
traindf = [];
validdf = [];
ctrainlist = {};
cvalidlist = {};
    
for i=1:length(filtdatalist)
    currdata = filtdatalist{i};
    currn = size(currdata,2);
    if currn<sampsize
        s = currn;
    else
        s = sampsize;
    end
    traininds = randsample(currn,s);
    otherinds = setdiff(1:currn, traininds);
    validinds = randsample(otherinds,nvalid,true); % with replacement for test indices
    
    traindf = [traindf, currdata(:,traininds)];
    validdf = [validdf, currdata(:,validinds)];
    
    ctrainlist{i} = traininds;
    cvalidlist{i} = validinds;

end
traindf = full(traindf);
validdf = full(validdf);
   
% run orthogonal NMF
%parfor(loopvar = initVal:endVal,useworkers);
%    k = klist(loopvar);

% First run NMF on k that needs computing
parfor (loopvar = initVal:endVal) 
    
    k = klistforNMF(loopvar);
    currname = featnames{loopvar};
    
    % Run orthogonal nmf
    [W, S, H, numIter, tElapsed, finalResidual] = orthnmfrule(traindf,k,opts);
   
    % Compute rescaling matrices for W and H
    a = []
    for i=1:k
        currw = W(:,i);
        ai = norm(currw);
        a = [a, ai];
    end
    g = size(W,1);
    A = repmat(1./a, g, 1)

    c = size(H,2)
    B = repmat(a',1,c)
    
    % rescale W's so that unit norm is 1
    % rescale H's 
    What = A.*W;
    Hhat = H.*B
    
    % solve for H for the validation data 
    validH = nnls(What,validdf');
    
    parsave(currname,What,Hhat,numIter,tElapsed,finalResidual, validH, validdf, filtdata, ctrainlist, cvalidlist)
end

% Now iterate through all learned feature sets and make plots of the
% features
resids = [];
for i=1:length(klist)
    %k = klistforNMF(i);
    
    k = klist(i);
    currname = featnames{i};
    
    % Load in saved features from parfor loop
    disp('Loading in saved features from parfor loop... ');
    load(currname)
    W = feats.W;
    
    % Run GSEA if not already run: 

    test = sum(strcmp(fieldnames(feats),'setnames'));
    if test==0

        genesets = gsea(feats);
        feats.setnames = genesets.setnames;
        feats.allsetnames = genesets.allsetnames;
        feats.sets_genes = genesets.sets_genes;
        feats.sets_p = genesets.sets_p;

        save(currname, 'feats');
    end
    
    % need to calculate gorder and then re-save 
    [gorder,featlist] = Worder(W, 4);
    feats.geneorder = gorder;
    save(currname, 'feats');
    
    
    % load in validation data and calculate residuals
    filtdata.genes = feats.genes;
    validH = feats.validH;
    validdf = feats.validdf;
    testresid = sqrt(sum(sum((validdf - W*validH).^(2))));    
    resids=[resids,testresid];
    
    
    % Plot error for each dataset individually
    plotErr(W*validH,validdf);
    print([featfigfolder, expname, '_W', num2str(k), '_err.pdf'],'-dpdf','-r0')
    close;
    
    % Plot the W and H separately 
    figure; 
    subplot(2,1,1); imagesc(W(gorder,:)); title('W - sorted');
    subplot(2,1,2); imagesc(validH); title('H - sample order');
    set(gca,'YTick', 1:size(validH,1), 'YTickLabel',feats.setnames,'TickLabelInterpreter','none');
    set(gca,'FontSize',8)
    currpos = get(gcf,'Position');
    newpos = currpos;
    newpos(3)=800;
    set(gcf,'Position',newpos);
    
    print([featfigfolder, expname, '_W', num2str(k), '_WHplot.pdf'],'-dpdf','-r0')

end
    

% plot the cost function for all these datasets
errcost = resids + (klist).^(1);
figure;plot(klist,errcost,'k*-');
xlabel('number of features');
ylabel('residual - cost(k)');
set(gca, 'FontSize', 16);
print([featfigfolder,expname,'_cost_vs_features.pdf'],'-dpdf','-r300')
close;

% Load feature set that minimizes the cost:
ind = find(errcost==min(errcost));
disp(['The best feature set has ', num2str(klist(ind)), ' features'])
name=[featfolder, expname, '_W', num2str(klist(ind)), '_features.mat'];
load(name)

function parsave(currname, What, Hhat, numIter,tElapsed,finalResidual, validH, validdf, filtdata, ctrainlist, cvalidlist)

% Find H for validdf
disp('Saving features from parfor loop... ');

feats = struct();
feats.expname = filtdata.expname;
feats.W = What;
feats.H = Hhat;
feats.numIter = numIter;
feats.tElapsed = tElapsed;
feats.finalResidual = finalResidual;
feats.genes = filtdata.genes;
feats.validH = validH;
feats.validdf = validdf;
feats.k = size(feats.W,2);
feats.ctrainlist = ctrainlist;
feats.cvalidlist = cvalidlist;
feats.feattype = 'ONMF';

save(currname, 'feats');


