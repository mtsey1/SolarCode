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