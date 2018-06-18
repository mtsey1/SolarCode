% p = secondPeak(data)
% p = secondPeak(data, num)
%
% [p, conf, sep, rect] = secondPeak (_)
%
% For each row of data, find the "secondary peak".  It is assumed that there
% is a "primary peak" near the start of the row, and a general downward trend.
% The second peak is the point after the primary peak that maximizes
%    d(i) - min_{j<i} d(j)
%
% If num is supplied, find the  n  largest secondary peaks.
%
% Return value conf is the "confidence" in each peak, currently measured as
% the objective above.
%
% Return value sep is the "separation" of the peak from the peak preceding it
% (sorted by value, not index in p).
% Currently the sum of the fraction of  intervening values less than 10% of
% the height of this peak and
% (1 - fraction of this peak of 1/4  of intervening values)
%
% Return value "rect" is the area of the largest rectangle that
% can fit between pairs of peaks such that it is entirely below the peaks,
% and entirely above the histogram.
function [peaks, conf, sep, rect] = secondPeak(data, num)
  if nargin < 2
    num = 1;
  end
  cm = zeros(size(data));
  r(1) = 0; r(2) = length(data)+1;	% region boundaries
  [c_pool(1), p_pool(1)] = max(data);

  for i = 1:num
    [conf(i), p] = max(c_pool);
    peaks(i) = p_pool(p);

    p_pool(p+1:end+1) = p_pool(p:end);
    c_pool(p+1:end+1) = c_pool(p:end);
    r(p+2:end+1) = r(p+1:end);
    r(p+1) = p_pool(p);

    cm(p_pool(p)) = data(p_pool(p));

          % Update cumulative min for points before the new peak
    range = r(p)+1:r(p+1)-1;
    cm_tmp = cummin(data(range(end:-1:1))); cm_tmp = cm_tmp(end:-1:1);
    cm(range) = max(cm(range), cm_tmp);
    if ~isempty(range)
	    [c_pool(p), p_pool(p)] = max(data(range)-cm(range));
	    p_pool(p) = p_pool(p) + r(p);
  	else
	    c_pool(p) = -Inf;
    end

        % Update cumulative min for points after the new peak
    range = r(p+1)+1:r(p+2)-1;
    cm(range) = max(cm(range), cummin(data(range)));
    if ~isempty(range)
	    [c_pool(p+1), p_pool(p+1)] = max(data(range)-cm(range));
	    p_pool(p+1) = p_pool(p+1) + r(p+1);
    else
	    c_pool(p+1) = -Inf;
    end

    % Avoid repeatedly reporting peaks
    %c_pool(find(ismember(p_pool,peaks))) = -Inf;

    if 0
      figure(123); plot(1:length(data), data, 1:length(data), cm);
      fprintf('r ');      fprintf('%d ', r);      printf('\n');
      fprintf('peaks ');  fprintf('%d ', peaks);  printf('\n');
      fprintf('conf ');   fprintf('%d ', conf);   printf('\n');
      fprintf('p_pool '); fprintf('%d ', p_pool); printf('\n');
      fprintf('c_pool '); fprintf('%d ', c_pool); printf('\n');
      keyboard;
    end
  end

  % calculate separation
  if nargout >= 3
    [pk, idx] = sort(peaks);
    sep(num) = 0;       % pre-allocate
    for i = 2:num
      d = data(pk(i-1)+1:pk(i));
      term1 = sum(d < 0.1 * data(pk(i)));
      term1 = term1 / (pk(i)-pk(i-1));

      term2 = sort(d);
      term2 = 1 - term2(ceil(end/4)) / data(pk(i));
      sep(idx(i)) = term1 + term2;
    end
  end

  % calculate biggest rectangle fitting between peaks
  % (This shows the gap is wide and deep)
  if nargout >= 4
    rect(2*num,3) = 0;        % pre-allocate
    for i = 2:num
      d = data(pk(i-1):pk(i));
      best = find_rect (d);
      rect(idx(i),:) = [best(1:2)+pk(i-1), best(3)];
    end

    for i = 2:num
      d = data(min (peaks(i-1:i)):max (peaks(i-1:i)));
      best = find_rect (d);
      rect(i + num,:) = [best(1:2)+min(peaks(i-1:i)), best(3)];
    end
  end
end
