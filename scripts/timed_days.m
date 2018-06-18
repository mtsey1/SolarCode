function [ st, en ] = timed_days (rectangles, r, cv, extend, shrink)
%TIMED_DAYS  Estimate start/end days of a timed device
  on_off = rectangles(r).on_off;
  burst = rectangles(r).burst;
  st = on_off(1);
  en = on_off(3);
  if nargin < 4
    extend = [1, 1];
  end
  if any (extend)
    oo = [rectangles.on_off];
    [overlap, st_en_offset] = overlap_hrs (on_off, rectangles);
    if extend(1)
      overlap_before = overlap & [oo(3,:)] <= st;
      if any (overlap_before)
        st = max (oo(3, overlap_before) + st_en_offset(overlap_before));
      else
        st = 1;
      end
    end
    if extend(2)
      overlap_after = overlap & oo(1,:) >= en;
      if any (overlap_after)
        en = min (oo(1, overlap_after) + st_en_offset(overlap_after));
      else
        en = size (cv, 2) + 1;
      end
    end
  end
  if nargin < 5
    shrink = [0,0];
  end
  if shrink(1)
    max_st = en;
  else
    max_st = on_off(1);
  end
  if shrink(2)
    min_en = st;
  else
    min_en = on_off(3);
  end
  days = st:en-1;
  if isfield (rectangles, 'alt_days') && ~isempty (rectangles(r).alt_days)
    m = mod1 (days, 7);
    is_alt = bsxfun (@eq, m, rectangles(r).alt_days);
    days = days (~any (is_alt, 1));
    if isempty (days)
      return;
    end
  end
  burst = burst (1:max (1, end-1));
  p = 0.8 * rectangles(r).power;

  if on_off(2) < on_off(4)
    patch = cv(burst, days);
  else
    patch=[cv(ceil (on_off(2)):end, days);
           cv(1:floor (on_off(4))-1, min (days+1,size (cv,2)))];
  end
  num_small = sum (patch < 0.9 * rectangles(r).power, 1);
  num_small = -rolling_min (-rolling_min (num_small'))';
  valid_from = find (num_small(1:on_off(1) - st + 1), 1, 'last');
  if ~isempty (valid_from)
    st = st + valid_from;
  end;
  valid_from = find (num_small(on_off(3) - st + 1:end), 1);
  if ~isempty (valid_from)
    en = on_off(3) + valid_from;
  end
  days = days (days >= st & days < en);
  
  average = mean (cv(burst, days), 1);
  pre  = mean (cv(mod1 (floor (on_off(2)) - (1:2), size (cv, 1)), days), 1);
  post = mean (cv(mod1 (ceil  (on_off(4)) + (0:1), size (cv, 1)), days), 1);
  average = [average; pre; post];
  
  jump = find_step (average);
  [j_max, mx] = max (jump(1,1:1 + max_st - st));
  [j_min, mn] = min (jump(1, min_en - st + 1:end));
  mn = min_en - st + mn;
  mid = floor (length (days) / 2);
  range_early = 1:mid;
  range_late  = mid+1 : length (days);
  if isempty (mn)
    mn = length (pre);
    j_min = 0;
  end
  if isempty (mx)
    mx = 1;
    j_max = 0;
  end
  
  if j_max == 0 && j_min == 0 ...
    return;       % entirely flat, so st and en are maximum values
%  elseif  corr(jump(1,:)', jump(2,:)').^2 +  corr(jump(1,:)', jump(3,:)').^2 < 0.5
%    if sqrt (var (average(1,:))) < 0.5 * mean (average(1,:))
%      return;
%    end
  end
  
  if mx < mn
    left = mean (average(:, 1:mx-1), 2);
    mid = mean (average(:, mx:mn), 2);
    ml = mid - left;
    if corr(jump(1,range_early)', jump(2,range_early)').^2 ...
        +  corr(jump(1,range_early)', jump(3,range_early)').^2 > 0.6 ...
        || sqrt (var (medfilt1 (average(1,range_early), 5))) > 0.25 * mean (average(1,range_early)) ...
        || left(1) < p
      if ((j_max > 1 && ml(1) > p) ...
            || (j_max > 3 && ml(1) > 0.8 * p)) ...
          && (left(1) < p || min (abs (ml(1) - ml(2:3))) > 0.5 * abs (ml(1)))
        st = days(mx);
      else
%        st = on_off(1);
      end
    end
    right = mean (average(:,mn+1:end), 2);
    mr = mid - right;
    if corr(jump(1,range_late)', jump(2,range_late)').^2 ...
        +  corr(jump(1,range_late)', jump(3,range_late)').^2 > 0.6 ...
        || sqrt (var (medfilt1 (average(1,range_late), 5))) > 0.25 * mean (average(1,range_late)) ...
        || right(1) < p
      if ((j_min < -1 && mr(1) > p) ...
            || (j_min < -3 && mr(1) > 0.8 * p)) ...
          && (right(1) < p || min (abs (mr(1) - mr(2:3))) > 0.5 * abs (mr(1)))
        en = days(mn);
      else
%        en = on_off(3);
      end
    end
  end
end
