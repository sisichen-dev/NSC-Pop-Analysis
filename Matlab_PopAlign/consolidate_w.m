function abmat = consolidate_w(alignedPops, gmmdata)
% function fig1 = plot_w(alignedPops, refpop, Hdata, rawdata, labeltype)
% This function collects abundances aligned to the reference subpopulations 

numrefcomp = length(unique(alignedPops.alignidx.ref_comp_idx));
nummodels = length(unique(alignedPops.alignidx.test_idx));
% make abundance matrix 
abmat = zeros(numrefcomp, nummodels);

for i=1:numrefcomp
    refpop = i;
    % extracts the original abundances
    origw = alignedPops.w;

    %% Find indexes of all data aligned to this subpopulation
    idx = find(alignedPops.alignidx.ref_comp_idx==refpop);
    modelidx = alignedPops.alignidx.test_idx(idx);
    modelnames = alignedPops.alignidx.Var1(idx);

    %%
    % Now consolidate the w's for each model (i.e. sum the components that
    % align to the same population) 

    w_subset = origw(idx);
    modelvals = unique(modelidx);
    summed_w = [];
    nameslist = {};
    for i=1:length(modelvals)
        currval = modelvals(i);
        totw = sum(w_subset(find(modelidx==currval)));
        % populate the abundance matrix
        abmat(refpop, currval) = totw;
        %summed_w = [summed_w, totw];
        currname = modelnames(find(modelidx==currval));
        nameslist{i} = currname{1};
    end

    orignames = alignedPops.dfnames;
    sampidx = cellfun(@(x) find(strcmp(orignames,x)), nameslist);
    
end

% try plotting this for one of the subpopulations

sampidx  = cellfun(@(x) find(strcmp(orignames,x)), gmmdata.modelnames);

% only keep sampidx==3:8
keep = find(sampidx>0)

figure;
for refpop = 1:numrefcomp
    subplot(1,numrefcomp,refpop)
    errorplot(sampidx(keep),abmat(refpop,keep),unique(keep),'LineWidth',3); hold on;
    scatter(sampidx(keep),abmat(refpop,keep),10, 'filled', 'Jitter', 'on', 'JitterAmount', 0.1)
    set(gca,'XTIck', unique(sampidx),'XTickLabel', orignames(unique(sampidx)), 'TickLabelInterpreter', 'none');
    set(gca,'FontSize',16)
    xtickangle(45)
    title(['Population ', num2str(refpop)])
    ylim([0,inf]);
end

% Get results and 
result.sampidx = sampidx;
result.modelvals = modelvals';
result.summed_w = summed_w;
result.names = nameslist;

