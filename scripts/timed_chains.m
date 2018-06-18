function [chains, dists, chain_value] = timed_chains (rectangles, pdfs, cv)
  % Identify sequences of rectangles corresponding to the first, second
  % or subsequent times that a pool pump is on each day.
  len = length (rectangles);
  on_off = [rectangles.on_off];
  slots_per_day = size (cv,1);
  
  duration = on_off(4,:) - on_off(2,:);
  at_midnight = (duration < 0);
  duration(at_midnight) = duration(at_midnight) + slots_per_day;
  
  centres = (on_off(4,:) + on_off(2,:)) / 2;
  am = (centres < slots_per_day / 2);
  centres(at_midnight & am) = centres(at_midnight & am) + slots_per_day/2;
  centres(at_midnight &~am) = centres(at_midnight &~am) - slots_per_day/2;

  % Calculate edge weights
  costs = inf (len + 2);  % len+1 = day 0, len+2 is day end+1
  
  [~, idx] = sort (on_off(1,:));
  inv_idx(idx) = 1:length (idx);
  
  ks_p = [rectangles.ks_p];
  ks_p = ks_p(1,:);
  ks_p(~isfinite (ks_p)) = mean (ks_p(isfinite (ks_p)));
  if ~isfinite (ks_p)   % if none finite, assign any equal value in (0,1).
    ks_p(:) = 0.5;
  end
  
  chain_costs = zeros (size (rectangles));
  
  for i = 1:len
    ii = idx(i);
    p = zeros (4, len - i);
    duration_ratio = duration(ii) ./ duration(idx(i+1:end));
    duration_ratio(duration_ratio<1) = 1 ./ duration_ratio(duration_ratio<1);
    p(1,:) = pdfs.duration_ratio (duration_ratio);
    
    day_gaps = on_off(1,idx(i+1:end)) - on_off(3, ii);
    p(2,:) = pdfs.day_gap (day_gaps);
    
    centre_changes = centres (ii) - centres(idx(i+1:end));
    centre_changes(centre_changes > slots_per_day/2) ...
      = centre_changes(centre_changes > slots_per_day/2) - slots_per_day;
    centre_changes = abs (centre_changes);
    p(3, day_gaps > 0)= pdfs.centre_change_gap   (centre_changes(day_gaps > 0));
    p(3, day_gaps<=0) = pdfs.centre_change_no_gap(centre_changes(day_gaps<=0));
    
    % Kolmogorov-Smirnov trust of next rectangle
    p(4,:) = -log (1 - ks_p(idx(i+1:end)));
    
    costs(ii,idx(i+1:end)) = mean (p, 1, 'omitnan');
    
    % Look for other rectangles that start/end at or near this transition
    other_start = find (abs (on_off(1,:) - on_off(3,ii)) < 3);
    other_end   = find (abs (on_off(3,:) - on_off(3,ii)) < 3 ...
                        & ii ~= 1:size (on_off,2));
    if length (other_start) + length (other_end) > 1
      factor = other_start;
      factor(:) = 0.9;
      
      if length (other_start) >= 2 && length (other_end) >= 1
        a = bsxfun (@minus, on_off(4, [ii, other_end]), on_off(2, other_start)');
        hrs2 = size(cv,1) / 2;
        a(a >  hrs2) = a(a >  hrs2) - size (cv, 1);
        a(a < -hrs2) = a(a < -hrs2) + size (cv, 1);
        
        % TODO: find best pairing
        best = abs (a(1,1) - a(2,2));
        next_best = abs (a(1,2) - a(2,1));
        if best < next_best
          pair_start = other_start(2);
          pair_end = other_end(1);
        else
          [other_start(1), other_start(2)] = deal (other_start(2), other_start(1));
          [best, next_best] = deal (next_best, best);
          pair_start = other_start(1);
          pair_end = other_end(1);
        end
        factor(1:2) = factor(1:2) .* (1 - 0.25 * exp (-[best, next_best]));
        
        dur_pre  = on_off(4, [ii, other_end]) - on_off(2, [ii, other_end]);
        dur_post = on_off(4, other_start)    - on_off(2, other_start);
        dur_pre (dur_pre  < 0) = dur_pre (dur_pre  < 0) + size (cv, 1);
        dur_post(dur_post < 0) = dur_post(dur_post < 0) + size (cv, 1);
        factor(1:2) = factor(1:2) * (1 - 0.25 * exp (-abs (dur_pre (1) - dur_pre (2))));
        factor(1:2) = factor(1:2) * (1 - 0.25 * exp (-abs (dur_post(1) - dur_post(2))));
        
        run_pre = on_off(3, [i, other_end]) - on_off(1, [i, other_end]);
        run_post= on_off(3, other_start)    - on_off(1, other_start);
        factor(1:2) = factor(1:2) * (1 - 0.25 * exp (-abs (run_pre (1) - run_pre (2))));
        factor(1:2) = factor(1:2) * (1 - 0.25 * exp (-abs (run_post(1) - run_post(2))));
      end
      costs(ii, other_start) = factor .* costs(ii, other_start);
    end
    
    % No need to permute these by idx
    %costs(end-1, i) = -log (pdfs.day_gap (on_off(1,i)) ...
    %                        * (1-rectangles(i).ks_p));
    %costs(i, end)   = -log (pdfs.day_gap (size(cv,2) + 1 - on_off(3,i)));
  end
  costs(end-1,1:end-2) = pdfs.day_gap (on_off(1,:)) - log (1 - ks_p(1));
  costs(1:end-2, end)  = pdfs.day_gap (size(cv,2) + 1 - on_off(3,:));
  
  % This is currently needed due to
  % a bug in duration_ratio_pdf in has_pool_pump.m
  if any (costs(:) < 0)
    fprintf ('timed_chains.m:  ERROR!!!  No edge costs should be negative\n');
    costs = max (0, costs);
  end
  
  [s, d] = find (isfinite (costs));
  g = digraph (s, d, costs(isfinite (costs(:))));
  
  unallocated = len;
  i = 1;
  chains = cell(1);                   % Avoid Matlab "preallocate" warning
  dists = zeros(1);                   % AVoid Matlab "preallocate" warning
  while unallocated > 0
    [ch, dists(i)] = shortestpath (g, len+1, len+2, 'method', 'acyclic');
    if length (ch) < 3
      dists = dists (1:i-1);
      break;
    end
    chains{i} = ch(2:end-1);
    
    % See how much each rectangle reduced the cost of the chain by.
    j = 2:length(ch)-1;
    x = size (costs, 2);
    if length (j) > 1
      chain_value(ch(j)) = sqrt (length (ch)) ...
                           * (costs(ch(j-1) + x*(ch(j+1)-1)) ...
                            - costs(ch(j-1) + x*(ch(j)-1)) ...
                            - costs(ch(j)   + x*(ch(j+1)-1)));
    else
      chain_value(ch(j)) = sqrt (length (ch)) * (pdfs.day_gap (size (cv,2)) ...
                                                 - dists(i));
    end
    unallocated = unallocated - (length (ch) - 2);
    i = i + 1;
    s = [];
    d = [];
    for j = 2:length(ch) - 1
      ii_ch = inv_idx(ch(j));
      s = [s, ch(j(ones(1, len+2 - ii_ch))), idx(1:ii_ch), len+1];
      d = [d, idx(ii_ch:end), len+2,         ch(j(ones(1, ii_ch+1)))];
    end
    g = rmedge (g, s, d);
  end
  if isempty (chains{1})
    chains = cell (0);
  end
end
