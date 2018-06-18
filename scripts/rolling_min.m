function out = rolling_min(data, window)
% out = rolling_min (data, window)
%
% out(i,:) = min(data(i + (-window/2 : window/2),:))
% If window is not specified, take  window=3.
    if nargin < 2
        window = 3;
    end
    out = zeros(size(data));
    window = min(window, size(data,1));
    if window == 0
        return
    end
    w2 = floor(window/2);
    ww = window - w2 - 1;

    % take min over sliding window of size  w
    for i = 1:w2
       out(i,:) = min(data(1:i+ww,:),[],1);
    end

    % The loop below is semantically the same as the following, but gathers
    % non-overlapping blocks and performs a single "min" on them
    % That should be more efficient on long vectors and short windows
    %           for i = 1+w2:size(data,1)-ww
    %               out(i,:) = min(data(i-w2:i+ww,:),[],1);
    %           end
    middle = size(data,1) - w2 - ww;
    count = floor(size(data,1) / window);
    
    batch = 1:(1+size(data,1)-count*window);
    per_batch = ceil(middle / window);

    batch_2 = (2+size(data,1)-count*window):window;
    per_batch_2 = per_batch - 1;

    for j = 1:2
	% O(window) loops, which is usually better than O(length(data)).
        for i = batch
            d = reshape(data(i:i-1+per_batch*window,:), ...
                        [window, per_batch*size(data,2)]);
            out(i+w2:window:size(data,1)-ww,:) = ...
                reshape(min(d,[],1), [per_batch, size(data,2)]);
        end
        batch = batch_2;
        per_batch = per_batch_2;
    end
    for i = size(data,1)+1-(ww:-1:1)
       out(i,:) = min(data(i-w2:end,:),[],1);
    end

    %% max over same sliding window  w  gives
    %% "out(x) >= y if x is in a neighbourhood of size  w  that is at least y"
    %for i = w2:-1:1
    %   tmp = max(out(1:i+ww,:),[],1);
    %   data(i,:) = max(tmp, data(1,:));
    %end
    %for i = 1+w2:size(out,1)-ww
    %    data(i,:) = max(out(i-w2:i+ww,:),[],1);
    %end
    %for i = size(out,1)+1-(ww:-1:1)
    %   tmp = max(out(i-w2:end,:),[],1);
    %   data(i,:) = max(tmp, data(end,:));
    %end

end
