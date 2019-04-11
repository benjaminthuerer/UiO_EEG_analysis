% perform cluster-based statistics on results from UiO_results.m
% 
% UiO_cluster(dataset1, dataset2, chanlocs, perm_s)
% 
% dataset1: 2d matrix with subjects / sessions in the first and channels 
%           (or what clusters should be computed on) in the second 
%           dimension
% dataset2: like dataset1 but of a different condition (contrast which
%           should be computed)
% chanlocs: EEG channellocations structure according to EEGlab
% perm_s: number of permutations (default: 5000)

% for questions:
% benjamin.thuerer@gmail.com

function [] = UiO_cluster_statistics(mData, groupA, groupB, chanlocs, perm_s)

if nargin < 4
    perm_s = 5000;
end

n = 5000;
[P_map, O_map] = permutationTest(mData, groupA, groupB, n, 0, 1);

% topoplot grids of P-values
[~,topo_P] = topoplot(P_map, chanlocs, 'noplot', 'on', 'conv', 'on');

% topoplot grids of observed differences
[~,topo_O,plotrad,xmesh,ymesh] = topoplot(O_map, chanlocs, 'noplot', 'on');

% building clusters
topo_P = topo_P <= .05;
[obs_P,nblobs] = bwlabeln(topo_P);

clustsum_obs = zeros(1,nblobs);

abs_topo_O = abs(topo_O);
for i=1:nblobs
    clustsum_obs(i) = sum(abs_topo_O(obs_P(:)==i));
end

%% Interaction effect:
if ~isempty(clustsum_obs)
    permP_map = ones(length(groupA),1);
    permF_map = ones(length(groupB),1);

    perm_set = mData;

    n_set = size(perm_set,1);
    percentages = perm_s/10:perm_s/10:perm_s;
    percentages2 = 10:10:100;

    perm_i = 1;
    max_perm = zeros(1,perm_s);
    disp('start permutation testing. This may take a while');

    while perm_i <= perm_s
        rpt = randperm(n_set);
        permData = perm_set(rpt,:);
        [permP_map, permO_map] = permutationTest(permData, groupA, groupB, n, 0, 1);

           % topoplot grids of P-values
        [~,perm_P] = topoplot(permP_map,chanlocs,'noplot','on', 'conv', 'on');

        % topoplot grids of F-values
        [~,perm_O] = topoplot(permO_map,chanlocs,'noplot','on');

        % building clusters
        perm_P = perm_P <= .05;
        [perm_P, perm_nblobs] = bwlabeln(perm_P);

        clustsum = zeros(1, perm_nblobs);
        perm_O = abs(perm_O);
        if find(perm_nblobs)
            for iii = 1 : perm_nblobs
                clustsum(iii) = sum(perm_O(perm_P(:) == iii));
            end
            max_perm(1, perm_i) = max(clustsum);
        end

        if ~isempty(find(percentages == perm_i))
            disp([num2str(percentages2(percentages == perm_i)) ' % done']);
        end
        perm_i = perm_i +1;

        clear rpt permData permP_map permO_map perm_map nblobs clustsum perm_nblobs
    end

    stat_cluster = sum(repmat(abs(max_perm)',1,length(clustsum_obs)) > repmat(clustsum_obs,length(max_perm),1), 1) / perm_s;
else
    stat_cluster = [];
end

% start plotting
if find(stat_cluster)
    combi = [];
    [~,idx] = find(stat_cluster < 0.05);
    sum_map = zeros(length(idx),size(topo_O,1),size(topo_O,2));
    for i = 1:length(idx)
        sum_map(i,:,:) = obs_P == idx(i);
    end
    sum_map = squeeze(sum(sum_map,1));
    sum_map( sum_map >= 1 ) = 100;
    combi(1,:,:) = sum_map;
    combi(2,:,:) = topo_O;
    combi = squeeze(sum(combi,1));
    h = figure;
    toporeplot(combi,'plotrad',plotrad,'xsurface',xmesh,'ysurface',ymesh,'chanlocs',chanlocs,'electrodes','on','headrad',0.54,'style','map');
    set(gca,'clim',[0 12]);
    if isempty(idx)
        title('observed difference: no sign clusters')
    else
        title(['observed difference: significant cluster']);
    end
else
    disp('no significant clusters');
end

end