function rectangles = find_drift (rectangles, cv, overwrite)
% Find gradual delay in the turn-on and turn-off times of device timers
% Assumes this drift is uniform for all rectangles.
  if isempty (rectangles)
    return;
  end

  if nargin < 3
    overwrite = false;
  end

  if size (cv, 2) < 365
    a = [rectangles.on_off];
    if any ([a(1,:), a(3,:)] > size (cv, 2) + 1) ...
        && ~any ([a(1,:), a(3,:)] < 365 - size (cv, 2))
      cv = [zeros(size (cv, 1), 365 - size (cv, 2)), cv];
    end
  end

  % TODO: Copy to generate training/tweak
  oo = [rectangles.on_off];
    % restore sanity
  oo(1:2,:) = max (oo(1:2,:), 1);
  oo([1,3],:) = min (oo([1,3],:), size (cv, 2) + 1);
  oo([2,4],:) = min (oo([2,4],:), size (cv, 1) + 1);
  oo([1,3],:) = sort (round (oo([1,3], :)));
  for i = 1:length (rectangles)
    rectangles(i).on_off = oo(:, i);
  end

  local_time  = cell (2, length (rectangles));
  reliability = cell (2, length (rectangles));
  all_times = cell (1, 2 * length (rectangles));
  i = 1;
  fprintf ('finding drift...\n');
  while i <= length (rectangles)    % WHILE because rectangles grows.
%  if size (cv, 2) < 365
%    a = [rectangles.on_off];
%    if any ([a(1,:), a(3,:)] > size (cv, 2) + 1) ...
%        && ~any ([a(1,:), a(3,:)] < 365 - size (cv, 2))
%      cv = [zeros(size (cv, 1), 365 - size (cv, 2)), cv];
    [times1, times2, rel1, rel2, weighted1, weighted2, new_rect] = fit_edge (rectangles(i), cv);
    while ~isempty (new_rect)
      disp (rectangles(i).on_off');
      rectangles(i) = new_rect(1);
      rectangles(end+1) = new_rect(2);
      [times1, times2, rel1, rel2, weighted1, weighted2, new_rect] = fit_edge (rectangles(i), cv);
      if ~isempty (new_rect)
%        keyboard
      end
    end

    local_time{1,i} = times1;
    reliability{1,i} = rel1;
    all_times{1, 2*i - 1} = weighted1;

    local_time{2,i} = times2;
    reliability{2,i} = rel2;
    all_times{1, 2*i} = weighted2;

    % for turn-on time
    %   for middle third
    %     Find locally best jump time for each day
    %     Estimate reliability of each
    %   for first and last third
    %     find days when the "middle" sample changes, based on the above
    %     for each region, find locally best jump and reliability
    %   Scale by reliability
    %   shift to mean 0, scale by reliability
    %   add to all_times
    % repeat for turn-off time

    i = i + 1;
  end
  oo = [rectangles.on_off];

  weighted = [all_times{:}];
  positive = weighted(2,:) > 0;
  negative = weighted(2,:) < 0;
  slope = (mean (weighted(1, positive), 'omitnan') / mean (weighted(2, positive), 'omitnan') ...
         + mean (weighted(1, negative), 'omitnan') / mean (weighted(2, negative), 'omitnan')) / 2;

  residual = mean (abs (weighted(1,:) - slope*weighted(2,:)), 'omitnan')...
           / mean (abs (weighted(2,:)), 'omitnan');

  reject_H0 = (residual / sqrt (size (weighted, 2)) < 0.13 * abs (slope));
  statistic = residual / sqrt (size (weighted, 2)) / abs (slope);
  display (statistic)

  % Robust least squares fit
  % hypothesis test "slope != 0" against "slope == 0"
  if reject_H0 || overwrite
    oo = [rectangles.on_off];
    for i = 1:length (rectangles)
      old_slope = 0;
      if isfield (rectangles(i), 'drift')
        old_slope = sum (rectangles(i).drift);
      end
      rectangles(i).drift = slope;
      rectangles(i).on_off = oo(:,i) ...
                             - (slope-old_slope) * 0.5 * (oo(3,i) - oo(1,i))...
                                     * [0; 1; 0; 1];
      % Set turn-on and turn-off times to those at the start of the rectangle
    end
  end

  if reject_H0
figure(1); show_rectangles (rectangles, cv);
    for i = 1:length (rectangles)
      % Try to merge with following rectangle
      following = find (oo(1,:) > oo(3,i) - 5 & oo(1,:) < oo(3,i) + 15);
      following = following(following ~= i);
      drift = rectangles(i).drift * (oo(1,following) - oo(1,i));
      next = (max (abs (oo([2,4],following) - oo([2,4],i) - drift)) < 1);
      if ~any (next)
        continue
      end

      next = following (next);
      [~, idx] = min (oo(3, next));
      next = next(idx);
      % skip if another timer settings in between
      if any (oo(1,:) >= oo(3,i) & oo(3,:) <= oo(1,next))
        continue;
      end

      % check that the intervening period seems to be pool pump too
      max_start = ceil (max (oo(2,i)+drift, oo(2,next)));
      min_end  = floor (min (oo(4,i)+drift, oo(4,next)));
      interior = cv(max_start:min_end, oo(3,i):oo(1,next));
      if any (interior(:) < rectangles(i).power)
        continue;
      end

      if any (max (abs (oo([2,4],next) - oo([2,4],i) - drift)) > 0.5)
        % TODO: look for actual step
        continue;
      end

      oo(:, next) = [oo(1,i); oo(2,i); oo(3, next); oo(4,i)];
      rectangles(next).on_off = oo(:,next);
      oo(1,i) = NaN;
    end
    rectangles = rectangles(oo(1,:) > 0 & oo(1,:) ~= oo(3,:));

%    figure(1); show_rectangles(rectangles, cv);
%    figure(2); imagesc (cv);
%    keyboard
  end


%     % for turn-on time
%     %   for middle third
%     %     Find locally best jump time for each day
%     %     Estimate reliability of each
%     %     Shift to mean 0
%     %     Scale by reliability
%     %     Least squares fit of slope
%     %   for first and last third
%     %     find days when the "middle" sample changes, based on the above
%     %     for each region, find locally best jump and reliability
%     %     shift ot mean 0, scale by reliability, least squares fit
%     % repeat for turn-off] time
%     % test hypothesis "slope != 0" against "slope == 0"?

end

function [times1, times2, reliability1, reliability2, weighted1, weighted2, new_rect] = fit_edge (rect, cv)
  on_off = rect.on_off;
  new_rect = [];
  if isfield (rect, 'drift')
    drift = rect.drift;
  else
    drift = 0;
  end
  % Initially, fit line to entire edge
  [times1, reliability1, weighted1, jumps1] = raw_times (on_off, drift, cv, true);
  [times2, reliability2, weighted2, jumps2] = raw_times (on_off, drift, cv, false);

  most_jump = sum ([jumps1; jumps2] ~= 0) >= 3;
  split = any (most_jump);

  if split
    % for now, just try to split in two.
    f = find (most_jump);
    [~, idx] = max (sum (abs ([jumps1(:, f); jumps2(:, f)])));
    split = f(idx) + 1;
    if split < 4 || split > on_off(3) - on_off(1) - 4
      split = false;
    else
      % Fit line to both before and after split.
      % If slopes are similar, and the discontinuity is large
      % then assume separate rectangles
      a = [times1(1:split-1); times2(1:split-1)];
      m1 = median(a, 2);
      [s1, md1] = slope (bsxfun (@minus, a, m1));

      a = [times1(split:end); times2(split:end)];
      m2 = median(a, 2);
      [s2, md2] = slope (bsxfun (@minus, a, m2));

      mean_diff = (m2 - m1) / (0.5 * length (times1));
      % Hypothesis test: is mean_diff == s1 (== s2)
      c1 = corr (jumps1');
      c2 = corr (jumps2');
      c = abs ([c1(1,2), c2(1,2)]);
      c(isnan (c)) = 0;
      c = c(1) + c(2);
      if c < 1.7 || ...
          (c < 1.9 && sum (abs ([jumps1(:,split); jumps2(:,split)])) < 7)
        mac = max (abs ([c1(1,2), c2(1,2)]));
        if abs (s1 - s2) > (md1 + md2 + 0.001) / 2 ...
            && (mac < 0.9 || min (max ([s1, md1; s1, md2])) > 0.001)
          split = false;
        elseif c < 0.7 && mac < 0.6
          split = false;
        else
          ratio = (abs (mean_diff - (s1 + s2)/2)) ./ (0.001 + abs (s1 - s2) + [md1; md2]) * c;
          ratio = (ratio < 2);
          if all (ratio)
            split = false;
          elseif any (ratio)
            oo = rect.on_off([2,4]);
            m1(ratio) = oo(ratio);
            m2(ratio) = oo(ratio);
          end
        end
      end

%{
      % Divide into four regions: each side of split is divided in half.
      % Estimate drift from the two pairs not spanning split.
      % Check if jump is bigger than drift would suggest
      d1 = floor (split / 2);
      d2 = split;
      d3 = split + floor ((on_off(3) - on_off(1) - split) / 2);
      m1 = mean ([times1(1:d1-1);  times2(1:d1-1)],  2);
      m2 = mean ([times1(d1:d2-1); times2(d1:d2-1)], 2);
      m3 = mean ([times1(d2:d3-1); times2(d2:d3-1)], 2);
      m4 = mean ([times1(d3:end);  times2(d3:end)],  2);

      drift = sum (((m2 - m1) / d1 + (m4 - m3) / (d3 - d2))) / 4;
      excess_change = abs (m3 - m2 - drift * (d3 - d1)/2);

      % Find new on/off times
      if sum (excess_change) < 0.5
        split = false;
      end
%}
    end
  end

  if split
    % generate new rectangles.
    figure(1); show_rectangles (evalin ('caller', 'rectangles'), cv)
    new_rect  = [rect; rect];
    new_rect(1).on_off(3) = new_rect(1).on_off(1) - 1 + split;
    new_rect(2).on_off(1) = new_rect(1).on_off(3);

    %new_rect(1).on_off([2,4]) = (m1 + m2) / 2;
    %new_rect(2).on_off([2,4]) = (m3 + m4) / 2;
    new_rect(1).on_off([2,4]) = m1;
    new_rect(2).on_off([2,4]) = m2;
  end

%{
  % See if edge should be split
  t = medfilt1 (times, 3);
  N = (1:length (t));
  x_cumsum = cumsum (t);
  x2_cumsum = cumsum (t.*t);
  cum_mean = x_cumsum ./ N;
  cum_var = x2_cumsum ./ N - cum_mean .* cum_mean;
  cum_sd = sqrt (cum_var ./ (N-1));

  rev = t(end:-1:1);
  rev_cumsum = cumsum (rev);
  rev2_cumsum = cumsum (rev.*rev);
  rev_mean = rev_cumsum ./ N;
  rev_var = rev2_cumsum ./ N - rev_mean .* rev_mean;
  rev_sd = sqrt (rev_var ./ (N-1));

  rev_mean = rev_mean(end:-1:1);
  rev_sd  = rev_sd(end:-1:1);

  gap = abs (cum_mean - rev_mean) - cum_sd - rev_sd;

  [g, idx] = min (gap);
  if g > 0
    % Probably split here
    disp (gap);
    disp (idx);
    disp (g);
    figure(3); plot (N, cum_mean, N, cum_mean + cum_sd, '--', N, cum_mean - cum_sd, '--', N, rev_mean, N, rev_mean + rev_sd, ':', N, rev_mean - rev_sd, ':')
  end
%}
%{
  positive = weighted(2,:) > 0;
  negative = weighted(2,:) < 0;
  slope = (mean (weighted(1, positive)) / mean (weighted(2, positive)) ...
         + mean (weighted(1, negative)) / mean (weighted(2, negative))) /2;
%}
  % TODO
  % Find times that the true edge may lie outside
  % the assumed slot

end

function [s, md] = slope (points)
  % Calculate average slope from left to right of zero-mean rows of points
  x = 1:size (points, 2);
  x = x - mean (x);
  mid = length (x) / 2;
  left = median (points(:, 1:floor(mid)), 2, 'omitnan') ./ mean (x(1:floor(mid)), 2, 'omitnan');
  right = median (points(:, ceil(mid+1):end), 2, 'omitnan') ./ mean (x(ceil(mid+1):end), 2, 'omitnan');
  s = mean ([left; right]);
  md = mean (abs ([left; right] - s));
end

function [times, reliability, weighted, jumps] = raw_times (on_off, drift, cv, is_on)
  if on_off(3) == on_off(1)
    times = [];
    reliability = [];
    weighted = [];
    jumps = [];
    return
  end

  times = zeros (1, on_off(3)-on_off(1));
  all_on = times;     % pixels in whole row when pump is assumed on
  all_off = times;
  all_part = times;
  reliability = times;

  nominal = on_off(2 * (2 - is_on));
  offsets = drift * (1:on_off(3) - on_off(1));
  changes = [0, (find (diff (floor (nominal + offsets)))), length(offsets)];
  % loop over points in different hour slots
  for i = 1:length (changes) - 1
    st = on_off(1) + changes(i);
    en_1 = on_off(1) + changes(i+1) - 1;
    idx = changes(i)+1:changes(i+1);
    nominal_mid = nominal + drift * ((st + en_1)/2);
    [pre, mid, post, has_mid] = edge_neighbours (nominal_mid, size (cv, 1));

    % If the on_off was integer, there was no "mid", so create one.
    if ~has_mid
      mid = [pre, post];
      pre  = mod1 (pre  - 1, size (cv, 1));
      post = mod1 (post + 1, size (cv, 1));
    end
    if ~is_on
      tmp = pre;
      pre = post;
      post = tmp;
    end
    off  = cv(pre,  st:en_1);
    on   = cv(post, st:en_1);
    part = mean (cv(mid, st:en_1), 1);

    normal  = (off <= part) & (part <= on);
    odd     = (off > on);
    clip_on = (part > on)  & ~odd;
    clip_off= (off > part) & ~odd;

    if abs (pre - post) > 3       % wraps midnight
      if pre > post
        pre = pre - size (cv, 1); % no longer used as index, so can be <0
      else
        post = post - size (cv, 1);
      end
    end

    times(idx(clip_on)) = pre + (post > pre);
    times(idx(clip_off)) = post;

    tn = (on(normal) - part(normal)) ./ (on(normal) - off(normal));
    if ~has_mid
      tn = 2 * tn;
    end
    if is_on
      tn = pre + 1 + tn;
    else
      tn = pre - tn;
    end
    times(idx(normal)) = tn;
    times(idx(odd)) = mean (times(idx(~odd)));      % only used for splitting

    all_on  (idx) = on;
    all_off (idx) = off;
    all_part(idx) = part;
    % Base reliability estimate mainly on how low "off" sample is.
    % It could also consider local noise in "on" value
    reliability(idx) = 1 - off / max (off);
    reliability(idx(odd)) = 0;
    reliability(idx(clip_on | clip_off)) = 0.5 * reliability (idx(clip_on | clip_off));
  end

  N = 1:(on_off(3) - on_off(1));
  iter = max (1, round (log(N(end)) / log(2)) - 2);
  gap = max (0, min (1, [all_on-all_part; all_part-all_off] / mean (all_on - all_off)));
  jmp = zeros (2, N(end), iter);
  for j = 1: iter
    filt = 2 ^ j - 1;
    t = medfilt1 (gap, filt, [], 2);
    rise = t(:, filt+1:end) - t(:, 1:end-filt);
    rise(abs (rise) < 0.1) = 0;
    pad = zeros (2, (filt - 1) / 2);
    jmp(:, :, j) = [pad, j * rise, pad, [0;0]];
  end
  j = sum (jmp, 3);
  raw_jumps = j;
  j(abs(j) < abs(mean(j, 2))) = 0;
  idx = j(:, 1:end-1) .* j(:, 2:end) < 0;
  j(idx) = 0;
  j([[false; false], idx]) = 0;
  if iter > 1
    j = medfilt1 (j, 3, [], 2);
  end
  jumps = j;

  reliability(~isfinite (times)) = 0;
  times(~isfinite(times)) = post;       % Get rid of NaNs

  % Give lower weight to less reliable samples
  weighted = [times; 1:(on_off(3) - on_off(1))];
  weighted = bsxfun (@minus, weighted, mean (weighted, 2, 'omitnan'));
  weighted = bsxfun (@mtimes, weighted, reliability);
end

function [pre, mid, post, a] = edge_neighbours (on_off_col, hr_day, a)
  % Find the row before/after an horizontal edge.
  % If the edge is almost an integer, these differ by 1.
  % If the edge is nar the middle of an interval, they differ by 2.
  if nargin < 3
    a = (abs (on_off_col - round (on_off_col)) > 0.2);
  end
  pre(~a) = round (on_off_col(~a)) - 1;
  mid(~a) = pre(~a);
  pre(a) = floor (on_off_col(a)) - 1;
  mid(a) = pre(a) + 1;
  post = mid + 1;

  pre  = mod1 (pre,  hr_day);
  mid  = mod1 (mid,  hr_day);
  post = mod1 (post, hr_day);
end

function m = mod1 (value, modulus)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Modulo operator, returning values in [1, modulus], not [0, modulus-1].
  m = mod (value-1, modulus) + 1;
end
