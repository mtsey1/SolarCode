function GMM = GMM_AIC(data, weights, fuzz)
%  mean_var_weight = GMM_AIC(data)
% [mean_var_weight, associations] = GMM_AIC(data)
%
% Find a Gaussian mixture model describing the rows of  data,  using the
% Akaike Information Criterion (AIC) to determine the number of clusters

    if nargin < 2 || isempty(weights)
        weights = ones(size(data,1), 1);
    end
    if nargin == 3
        weights(:,2) = fuzz;
    else
        weights (:,2) = sum (var (data)) / 5;
    end

    GMMs = cell(1,2);           % keep best GMM without copying
    min_AIC = Inf;
    argmin_AIC = 1;
    k = 0;
    best = 2; working = 1;
    nCols = size(data,2);
    nRows = size(data,1);

    while k < 3 || (argmin_AIC > k/2 && k < size(data,1)/2)
        k = k+1;
        p00 = params(k, nCols, false, false);
        p01 = params(k, nCols, false, true);
        p10 = params(k, nCols, true, false);
        p11 = params(k, nCols, true, true);
        if (p00 < nRows)
          sharedCovar = false;
          diagonalCovar = 'full';
        elseif (p01 < nRows)
          sharedCovar = false;
          diagonalCovar = 'diagonal';
        elseif (p10 < nRows)
          sharedCovar = true;
          diagonalCovar = 'full';
        elseif (p11 < nRows)
          sharedCovar = true;
          diagonalCovar = 'diagonal';
        else
          break
        end
        %[GMMs{working}, extra] = w_fitgmdist(data,k, 'start', 'plus', 'RegularizationValue', 1e-6, 'Replicates', 5, 'weights', weights);
        [GMMs{working}, extra] = w_fitgmdist(data,k, 'start', 'plus', 'RegularizationValue', 1e-6, 'weights', weights, 'sharedCovariance', sharedCovar, 'covarianceType', diagonalCovar);
        [k extra.AIC]
        if extra.AIC < min_AIC
            min_AIC = extra.AIC;
            argmin_AIC = k;
            best = working;
            working = 3-best;
        end
    end

    GMM = GMMs{best};
end

function p = params(k, nCols, sharedCovar, diagonalCovar)
  if (diagonalCovar)
    p = nCols;
  else
    p = nCols * (nCols+1) / 2;
  end
  if (sharedCovar)
    p = p*k;
  end
  p = p + 2*nCols - 1;
  
end

%
%    if isempty(weights)
%        weights = ones(size(data,1), 1);
%    end
%
%    i = 1;
%    mvw(i).means = w_mean(data, weights)
%    mvw(i).vars  = w_cov(data, weights, mvw(i).means);
%    mvw(i).weights = 1;
%    max_AIC = AIC(data, weights, mvw(i));
%    argmax_AIC = i;
%    
%    while argmax_AIC > i/2
%       new_p = addClust(mvw(i), data, weights);
%       i++;
%       mvw(i) = GMM(data, weights, new_p);
%       a = AIC(data, weights, mvw(i));
%       if a > max_AIC
%           max_AIC = a;
%           argmax_AIC = i;
%       end
%    end
%end
%
%% Weighted mean
%function wm = w_mean(data, weights)
%    wm = sum(bsxfun(@times, data, weights)) / sum(weights);
%end
%
%% Weighed variance
%function wv = w_cov(data, weights, means)
%    d = bsxfun(@minus, data, means);
%    wv = d' * bsxfun(@times, d, weights);
%end
%
%% form a new set of  n+1  clusters from a set of  n  clusters.
%% Currently, this is by splitting the cluster with the largest
%% (weighted) covariance matrix along its principal component
%function new_p = addClust(old_par, ~, ~)
%    new_p = old_par;
%    [~, toSplit] = max(trace(old_par.var) .* old_par.weight);
%    [splitDir, step] = pcacov(old_par.var(toSplit,:,:));
%                               % 0.798 = mean of positive half of std normal
%    new_p.means(toSplit)=old_par.means(toSplit,:)-0.798*step(1)*split_Dir(:,1)];
%    new_p.means(end+1)  =old_par.means(toSplit,:)+0.798*step(1)*split_Dir(:,1)];
%    step(1) = 0.364 * step(1); % 0.364 = var of one side of a standard normal
%    new_p.vars(toSplit) = split_Dir * step * splitDir';
%    new_p.weights(toSplit) = old_parweights(toSplit)/2;
%
%    new_p.vars   (end+1) = new_par.vars(toSplit);
%    new_p.weights(end+1) = new_par.weights(toSplit);
%end
%
%
