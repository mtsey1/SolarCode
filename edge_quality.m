function [top, bottom] = edge_quality (rect, cv)
  st = rect.on_off(1);
  en = rect.on_off(3);
  on = rect.on_off(2);
  off = rect.on_off(4);
%   if en - st > 35
%     partitions = 5;
%   else
%     partitions = 2;
%   end
  
  if ~isfield (rect, 'alt_days') || isempty (rect.alt_days)
    top    = find_quality (st, en, on,  cv);
    bottom = find_quality (st, en, off, cv);
  else
    top    = find_quality (st, en, on,  cv, rect.alt_days);
    bottom = find_quality (st, en, off, cv, rect.alt_days);
  end
end

function q = find_quality (st, en, edge, cv, omit)
  range = st:en-1;
  if nargin >= 6
    for i = 1:length(omit)
      range = range (mod (range, 7) ~= omit(i));
    end
  end
  
  [pre, mid, post, has_mid] = edge_neighbours (edge, size (cv, 1));

  if has_mid
    a = cv(pre, range);
    b = cv(mid, range);
    c = cv(post,range);
    jump = c - a;
    jump = (jump + mean (jump)) / 2;  % If a(i)=b(i)=c(i), make pos(i)=0
    pos = (c - b) ./ jump;
  else
  % If the on_off was integer, there was no "mid", so create one.
    a = mod1 ([pre-1, pre, post, post+1], size (cv, 1));
    sgn = sign (mean (cv(post, range) - cv(pre, range)));
    a = cv(a, range);

    jump = [a(4,:)-a(2,:); a(3,:)-a(1,:)];
    [jump, idx] = max (jump * sgn);
    jump = jump / sgn;
    jump = (jump + mean (jump)) / 2;

    pos = (a(2,:) - a(1,:)) ./ jump;
    idx = (idx ~= 2);
    pos(idx) = 1 + (a(3,idx) - a(2,idx)) ./ jump(idx);
    pos = pos / 2;
  end
  valid = pos >= 0 & pos <= 1;
  if ~any (valid)
    [~, valid] = min (abs (pos - 0.5));
  end
  pos = pos(valid);
  
  % Test if positions are nearly uniform.
  % That is evidence that it is just noise, not a real edge
  s = sort (pos);
  offset = (s - (1:length(pos)) / length (pos));
  ks = max (abs (offset));
  figure; plot (s);
  
  % look for multiple possible edge times
  zero_crossings = find (s(1:end-1) .* s(2:end) <= 0 & s(2:end) ~= 0);
  zero_crossings = [1 zero_crossings+1 length(pos)];

  if s (floor ((zero_crossings(1) + zero_crossings(2))/2)) > 0
    start = 1;
  else
    start = 2;
  end
  
  j = 1;
  a = zeros(4, ceil (length(zero_crossings-start-2)/2));
  for i = start+1:2:length(zero_crossings)-1
    [pre, left] = max (s(zero_crossings(i-1) : zero_crossings(i)));
    [post, right] = min (s(zero_crossings(i) : zero_crossings(i+1)));
    a(:, j) = [pre; -post; left+zero_crossings(i-1)-1; right+zero_crossings(i)-1];
    j = j + 1;
  end
  
  q = ks * sum (valid) / length (valid);
end