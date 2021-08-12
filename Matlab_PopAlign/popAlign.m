function alignedPops = popAlign(gmmdata, feats, strategy, refindex)
% function alignedPops = popAlign(gmmdata, feats, strategy, refindex)
%  This function takes in a gmmdata structure, in which the reference model
%  is specified, and aligns all other models to the referencemodel. If a
%  newrefmodel is specific, then we use the newrefmodel instead of the one
%  specified within the gmmdata structure. strategy is a string that
%  specifies which strategy we need to use:  'ref2test' (merge), 'test2ref' (bifurcate),
%  'conservative'(default)
% 
%  ref: model with k components
%  testlist: list of j models
%  each test model has i components
%
%  Assume n is the number of total aligned subpopulations
% 
%  The returned structure align has the following fields:
%  alignedPops.JDmatlist [ 1 x m cell array ] : pairwise JD matrix for each
%                                               test sample with ref sample
%
%  alignedPops.alignidx [ n x 3 table ] : has the following structure
%                                           [ref component number , 
%                                            test sample # , 
%                                            test component # ]
%  alignedPops.JDlist [ 1 x n array ] : JD values for aligned pairs
%  alignedPops.delta_mu [ 1 x n array ] : delta_mu for aligned pairs
%  alignedPops.delta_cov[ 1 x n array ] : delta_cov for aligned pairs
%  alignedPops.delta_w [ 1 x n array ] : delta_w for aligned pairs
%  alignedPops.subpopEnts [ 1 x n array ] : subpopEnts for aligned pairs
%  alignedPops.testnames [ 1 x n cell array ] : test sample name
%  alignedPops.refcompnames [ 1 x n cell array ] : ref component name
%  alignedPops.testcompnames [ 1 x n cell array ] : test component name
 

if nargin==2
    strategy = 'conservative';
    refindex = gmmdata.refmodel;
elseif nargin==3
    % check to make sure strategy is one of allowed strings: 
    test = startsWith({'ref2test','test2ref', 'conservative'}, strategy);
    if sum(test)~=1
        error('Strategy must be one of three options: test2ref, ref2test, or conservative')
    end
    refindex = gmmdata.refmodel;
end

% set the dimension along which minmask will operate
if (startsWith('ref2test',strategy))
    dim = 1;
elseif (startsWith('test2ref',strategy))
    dim = 2;
elseif (startsWith('conservative',strategy ))
    dim = 3;
end

% Get reference index

ref = gmmdata.gmmlist{refindex};
refname = gmmdata.refname;

testlist = gmmdata.gmmlist;
%testlist(refindex)=[];
testnames = gmmdata.modelnames;
%testnames(refindex)=[];
sampnumlist = gmmdata.sampnums;
%sampnumlist(refindex)=[];

W = feats.W;

% Construct JD matlist first: 
JDmatlist = cell(1,length(testlist));

for j=1:length(testlist)
    test = testlist{j};
    
    currJDmat=[];
    for k=1:ref.NumComponents
        mu0 = ref.mu(k,:)';
        sigma0 = ref.Sigma(:,:,k);
        
        %initialize JD vector 
        JD = zeros(1,test.NumComponents);
        
        for i=1:test.NumComponents

            % calculate alignment score and covariance distance
            AS = zeros(1,ref.NumComponents);

            mu1 = test.mu(i,:)';
            sigma1 = test.Sigma(:,:,i);
            w1=test.ComponentProportion(i);
    
            curr_jd = gkld(mu0,sigma0,mu1,sigma1);
        	JD(i) = curr_jd;
       
        end
        currJDmat = [currJDmat;JD];
    end
        
    JDmatlist{j} = currJDmat;
end

% find best alignment based on the strategy: 

% compute alignments matrix: 
alignidx = [];
testnames = {};
for j = 1:length(testlist)
    currJDmat = JDmatlist{j};
    % get masking matrix from JD matrix. All 1 entries are alignments to keep
    curr_mask = minmask(JDmatlist{j},dim);
    % find indexes of aligned components
    [ref_compidx,test_compidx] = find(curr_mask==1);
    
    npairs = length(ref_compidx);
    
    % store the sample id number
    test_idx = repmat(j,1, npairs);
    
    % keep the sample names also
    testnames = [testnames;cellstr(repmat(gmmdata.modelnames(j),npairs,1))];
    ref_idx = repmat(refindex,1, npairs);
    
    curr_Aidx = [ref_idx; ref_compidx';test_idx; test_compidx'];
    alignidx = [alignidx, curr_Aidx];
end

% compute deltas: 

JDlist = [];
delta_mu_list = [];
testmu_list = [];

delta_cov_list = [];
delta_w_list = [];

w_list = [];

refcompnames_list = {};
testcompnames_list = {};
testnames_list = {};


for y = 1:size(alignidx,2)
    % find positions of current population pair
    curr_refcompidx = alignidx(2,y);
    curr_testidx = alignidx(3,y);
    curr_testcompidx = alignidx(4,y);
    
    % transform curr_testidx to the position within the list
    curr_j = curr_testidx;
    
    test = testlist{curr_j};
    testname = testnames{y};
    
    % get alignment score
    currJD = JDmatlist{curr_j}(curr_refcompidx,curr_testcompidx);
    JDlist = [JDlist, currJD];
    
    % calculate delta mu
    refmu = ref.mu(curr_refcompidx,:)';
    testmu = test.mu(curr_testcompidx,:)';
    delta_mu = mumetric(refmu, testmu, W, feats.feattype);
    delta_mu_list = [delta_mu_list, delta_mu];
    testmu_list = [testmu_list, testmu];
    
    % calculate delta cov
    refSigma = ref.Sigma(:,:,curr_refcompidx);
    testSigma = test.Sigma(:,:,curr_testcompidx);
    delta_cov = covmetric(refSigma, testSigma);
    delta_cov_list = [delta_cov_list,delta_cov];
    
    % calculate delta w
    ref_w = ref.ComponentProportion(curr_refcompidx);
    test_w = test.ComponentProportion(curr_testcompidx);
    delta_w = test_w - ref_w;
    delta_w_list = [delta_w_list, delta_w];
    
    % keep original w
    w_list = [w_list, test_w];
    
    % refcomponentname 
    refcompnames_list = [refcompnames_list, {[refname, '-', num2str(curr_refcompidx)]}];

    % test component name 
    testcompnames_list = [testcompnames_list, {[testname, '-', num2str(curr_testcompidx)]}];
    
end

%% set all variables into return structure

 alignedPops.JDmatlist = JDmatlist;

 % add reference idx to alignidx
 temptable = array2table(alignidx', 'VariableNames',{'ref_idx','ref_comp_idx', 'test_idx', 'test_comp_idx'});
 temptable = [testnames,temptable];
 temptable.JD = JDlist';
 alignedPops.alignidx = temptable;
 alignedPops.JDlist = JDlist;
 alignedPops.delta_mu = delta_mu_list;
 alignedPops.testmu_list = testmu_list;
 alignedPops.delta_cov = delta_cov_list;
 alignedPops.delta_w = delta_w_list;
 alignedPops.w = w_list;
 alignedPops.refcompnames = refcompnames_list;
 alignedPops.testcompnames = testcompnames_list;
 alignedPops.strategy = strategy;
 alignedPops.refindex = refindex;
 alignedPops.dfnames = gmmdata.dfnames;
    