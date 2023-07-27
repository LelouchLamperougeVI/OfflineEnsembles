function hclust(obj)
% perform agglomerative clustering

if isempty(obj.ensembles.R)
    obj.corr;
end

Dm = 1 - obj.ensembles.R;
Dm(1:(length(Dm) + 1):numel(Dm))=0;
D = squareform(Dm);

obj.ensembles.tree = linkage(D,'average');
try
    obj.ensembles.order = optimalleaforder(obj.ensembles.tree, D);
catch
    tmp_d = D;
    tmp_d(isnan(tmp_d)) = 1;
    obj.ensembles.order = optimalleaforder(obj.ensembles.tree, tmp_d);
end

c = cluster(obj.ensembles.tree,'cutoff',obj.ops.thres,'criterion','distance');

clust = c(:) == (1:max(c));
clust(:, sum(clust) < obj.ops.e_size) = [];
[r, c] = find(clust);
clust = mat2cell(r, histcounts(c));

obj.ensembles.clust = clust;