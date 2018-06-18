function [guessScalar, guessCat] = knnaverage(idx, D, scalar, categorical)
weights = repmat(size(idx,2):-1:1, [size(idx,1), 1]);

guessScalar = weightedMean(idx, weights, scalar);
guessCat = weightedMode(idx, weights, categorical);

%guess = [guessScalar, guessCat];

end

% jth element of row i of guess is a mean, weighted by weights(i,:),
% of the elements of scalar(idx(i),j).
% (Should remove  for  loops, but I just want to get it working...)
function guess = weightedMean(idx, weights, scalar)
guess = zeros(size(idx,1), size(scalar,2));
for i = 1:size(idx,1)
    for j = 1:size(scalar,2)
        guess(i,j) = weights(i,:) * scalar(idx(i,:),j) ;
    end
    guess(i,:) = guess(i,:) / sum(weights(i,:));
end
end

function guess = weightedMode(idx, weights, categorical)
guess = zeros(size(idx,1), size(categorical,2));
for i = 1:size(idx,1)
    for j = 1:size(categorical,2)
        c = categorical (idx(i,:),j);
        values = unique(c);
        count = zeros(size(values));
        for k = 1:length(values)
            count(k) = sum(weights(i, c == values(k)));
        end
        [~, pos] = max(count);
        guess(i,j) = values(pos);
    end
end
end
