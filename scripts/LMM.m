% Laplacian mixture model
%   clusters = LMM(data, k)
% Clusters is an array, with each row being a cluster with columns:
% First:  the mean of the cluster
% Second: the variance of the Laplacian
% Third:  the weight of the cluster
%
% Currently a stub, using k-means instead
function clusters = LMM(data, k)
    [assign, centroids, dists] = kmeans(data, k);
    clusters(:,3) = hist(assign,1:k) / length(data,1);
    clusters(:,2) = dists ./ clusters(:,3);
    clusters(:,1) = centroids;
end
