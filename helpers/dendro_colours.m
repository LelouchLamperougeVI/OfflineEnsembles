function I = dendro_colours(Z, clust)
% find the indices of the cluster in the dendrogram

I = [];
if(length(clust)>1)
    [~, i] = intersect(Z(:,1), clust);
    [~, j] = intersect(Z(:,2), clust);
    idx = intersect(i,j);
    
    mask = nan(size(Z,1),1);
    mask(idx) = 1;
    
    [~,idx] = min(Z(:,3) .* mask);
    
    clust(clust == Z(idx,1)) = [];
    clust(clust == Z(idx,2)) = [];
    clust = [clust; idx + size(Z,1) + 1];
    
    I = [I idx];
    I = [I dendro_colours(Z, clust)];
end