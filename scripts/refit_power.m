basedir = 'C:/Users/Lachlan Andrew/Google Drive/FIT1041/pools/data/';

for house = 1:4000
  % rect_name       = sprintf ([basedir, 'rect%05d.mat'], house);
  saved_rect_name = sprintf ([basedir, 'saved_rect%05d.mat'], house);
  house_name      = sprintf ([basedir, 'house%05d'], house);
  if ~exist (saved_rect_name, 'file')
    continue;
  end
  load (house_name);
  load (saved_rect_name);
  if isempty (rect)
    continue;
  end
  rectangles = rect;

  % image without days for which no measurements exist
  valid_days = find (~isnan (image(1,:)));
  cv = image (:, ~isnan (image(1,:)));
  inv = nan(1,365);
  inv(valid_days) = 1:length (valid_days);
  
  valid = true (1, length (rectangles));
  on_off = [rectangles.on_off];
  for i = 1:length (rectangles)
    on_off([1,3], i) = round (on_off([1,3], i));
    rectangles(i).on_off([2,4]) = mod1 (on_off([2,4], i), size (cv, 1));
    if on_off(3,i) > length (valid_days)
      on_off(3,i) = length (valid_days);
    end
    if on_off(1,i) < 1
      on_off(1,i) = 1;
    end
    if on_off(3,i) - 1 < on_off(1,i)
      fprintf ('Start after end for rectangle %d, house %d.  Skipping\n', i, house);
      valid(i) = false;
      continue;
    end
    rectangles(i).on_off([1,3]) = valid_days(on_off([1,3], i));
  end
  rectangles = rectangles(valid);

  % Update powers.  ~Min-link clustering on ranksum.
  jumps = cell (length (rectangles), 0);
  weights = {};
  classes = zeros (1, length (rectangles));
  on_off = [rectangles.on_off];
  on_off([1,3],:) = round (on_off([1,3],:));
  on_off([2,4],:) = mod1 (on_off([2,4],:), size (cv, 1));
  [pre_up, ~, post_up] = edge_neighbours (on_off(2,:), size (cv, 1));
  [pre_dn, ~, post_dn] = edge_neighbours (on_off(4,:), size (cv, 1));
  valid = true (1, length (rectangles));
  for attempt = 1:2
    for i = 1:length (rectangles)
      up = diff (cv([pre_up(i), post_up(i)], inv(on_off(1,i)):inv(on_off(3,i)-1)));
      dn = diff (cv([post_dn(i), pre_dn(i)], inv(on_off(1,i)):inv(on_off(3,i)-1)));
      if isempty (up)
        valid(i) = false;
        continue;
      end
      done = false;
      for j = 1:length (jumps)
        p = ranksum (jumps{j}, [up, dn]);
        if p > 0.001
          done = true;
          jumps(j) = {[jumps{j}, up, dn]};
          classes(i) = j;
          % TODO: update weights
          break
        end
      end
      if ~done
        jumps(end+1) = {[up, dn]};
        classes(i) = length (jumps);
        % TODO: update weights
      end
    end
    if any (valid)
      break;
    end
    inv = 1:365;
  end

  rectangles = rectangles(valid);
  classes = classes(valid);

  for i = 1:length (jumps)
    jumps{i} = mean (jumps{i});
  end

  old_power = [rectangles.power];
  old_cor = set_cor_pool (rectangles, 1:365);
  
  for i = 1:length (rectangles)
    st = rectangles(i).on_off(1);
    en = rectangles(i).on_off(3);
    slice = cv(rectangles(i).burst(1:end-1),inv(st):inv(en-1));
    missed = st - 1 + find (min(slice,[],1) < 0.8 * jumps{classes(i)});
    if length (missed) < 0.25 * (en - st)
      rectangles(i).power = jumps{classes(i)};
    end
  end
  
  disp (house);
%%{
  disp (old_power);
  disp ([rectangles.power]);
  new_cor = set_cor_pool (rectangles, 1:365);
  figure(1); show_rectangles (rectangles, image);
  figure(2); plot (image);
  figure(3); imagesc (image + old_cor');
  figure(4); imagesc (image + new_cor');
%}
  %keyboard
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
