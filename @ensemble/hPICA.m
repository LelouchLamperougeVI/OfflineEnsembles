function hiepi = hPICA(obj)
% Core hierarchical PCA/ICA (hPICA) method

X = obj.twop.deconv;
if ~any(isnan(X), 'all')
    warning('Deconv does not contain NANs; possible omission of movement rejection.');
end

if isempty(obj.ensembles.clust)
    error('No ensembles exist. Possible omission of ensemble detection.');
end

% step 0. get rid of NANs and normalize to null mean, unitary variance
nan_idx = any(isnan(X), 2);
X(nan_idx, :) = []; % get rid of movement epochs
X = (X - mean(X)) ./ std(X); % z-score X to null mean and unit variance

% step 1. find initial estimate of basis vectors
init = zeros(size(X, 2), length(obj.ensembles.clust)); % construct initial estimate of mixing matrix W tilde
for ii = 1:length(obj.ensembles.clust)
    init(obj.ensembles.clust{ii}, ii) = 1;
end
init = init ./ vecnorm(init);

% step 2. estimate fraction of significant variance that need to be
% accounted for using Marchenko-Pastur theorem
[~, L] = eig(cov(X));
L = diag(L);
thres = march_pastur_lambda(X);
r_thres = size(X, 2) - sum(L(L < thres)); % the fraction of significant variance

% step 3. find eigenvectors of residual covariance matrix, get rid of the
% zero eigenvalue associated vectors and concatenate to basis estimate from
% hierarchical clustering to account for significant variance
resid = X - X * (init * init'); % residuals
[S, L] = eig(cov(resid));
[~, idx] = sort(diag(L));
S(:, idx(1:size(init, 2))) = [];
pc_comps = cat(2, init, S);
r = var(X * pc_comps);
[r((size(init, 2) + 1):end), idx] = sort(r((size(init, 2) + 1):end), 'descend');
pc_comps(:, (size(init, 2) + 1):end) = pc_comps(:, idx + size(init, 2));
r = cumsum(r);
W_ = pc_comps(:, r < r_thres);
if size(W_, 2) < size(init, 2)
    W_ = init;
end

X_ = X * W_;
md = rica(X_, size(W_, 2), 'InitialTransformWeights', eye(size(W_, 2)), 'ContrastFcn', 'logcosh');
W = md.TransformWeights' * W_';
W = W';
pc = W(:, 1:size(init, 2));

hiepi.pc = pc;
hiepi.W = W;
hiepi.init = init;

hiepi.X = nan(size(obj.twop.deconv));
hiepi.X(~nan_idx, :) = X;
hiepi.z = nan(size(obj.twop.deconv, 1), size(pc, 2));
hiepi.z(~nan_idx, :) = X * pc;

obj.hiepi = hiepi;


function l = march_pastur_lambda(X)
% Upper bound of eigenvalue defined by Marchenko-Pastur Law
dim = size(X);
l = (1 + sqrt(dim(2) / dim(1))) ^ 2;

