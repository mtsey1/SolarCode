ground_dir = '../pools/pools_keep/';
slots_per_day = 48;

files = dir ([ground_dir, '*.mat']);
len = length (files);

mean_duration = zeros (1, len);
var_duration  = zeros (1, len);

mean_hours_per_day = zeros (1, len);
var_hours_per_day  = zeros (1, len);

list_day_gaps = cell (1, len); % done

list_duration_ratio = cell (1, len); % done

list_change_centres_no_gap = cell (1, len); % done
list_change_centres_gap = cell (1, len); % done
var_change_centres = zeros (1, len);

list_runlengths = cell (1, len);

for i = 1:len
  a = load ([ground_dir, files(i).name]);
  rectangles = a.r;
  on_off = [rectangles.on_off];
  if isempty (on_off)
    continue
  end
  on_off([1,3], :) = round (on_off([1,3], :));

  % Calculate per-rectangle quantites
  duration = on_off(4,:) - on_off(2,:);
  at_midnight = (duration < 0);
  duration(at_midnight) = duration(at_midnight) + slots_per_day;
  
  centres = (on_off(4,:) + on_off(2,:)) / 2;
  am = (centres < slots_per_day / 2);
  centres(at_midnight & am) = centres(at_midnight & am) + slots_per_day/2;
  centres(at_midnight &~am) = centres(at_midnight &~am) - slots_per_day/2;
  
  % If pump is on twice per day, split into the two chains
  [~, idx] = sort (on_off(1,:));
  if length (idx) > 1
    overlap = on_off(1,idx(2:end)) < on_off(3, idx(1:end-1));
  else
    overlap = false;
  end
  if any (overlap(:))
    % Identify chains 1 and 2
    penultimate = length (idx) - 1;
    % An initial run of rectangles in chain 1 with no match in chain 2
    j = find (overlap, 1);
    if j == 1
      chain1 = idx(1);
      chain2 = idx(2);
    else
      chain1 = idx (1:j-1);

      % The first two overlapping rectangles
      centre_changes = centres (idx(j-1)) - centres(idx(j:j+1));
      centre_changes(centre_changes > slots_per_day/2) ...
        = centre_changes(centre_changes > slots_per_day/2) - slots_per_day;
      centre_changes = abs (centre_changes);
      if centre_changes(1) < centre_changes(2)
        chain1 = [chain1, idx(j)];
        chain2 = idx(j+1);
      else
        chain1 = [chain1, idx(j+1)];
        chain2 = idx(j);
      end
    end
    j = j + 1;
    while j < penultimate
      j = j + 1;
      if abs (on_off(1,idx(j)) - on_off(1,idx(j+1))) < 3
        % append closest in time of day
        centre_changes = bsxfun (@minus, centres([chain1(end), chain2(end)])', ...
                                         centres(idx(j:j+1)));
        centre_changes = mod (centre_changes, slots_per_day);
        centre_changes(centre_changes > slots_per_day/2) ...
          = centre_changes(centre_changes > slots_per_day/2) - slots_per_day;
        centre_changes = abs (centre_changes);
        if centre_changes(1) + centre_changes(4) ...
            < centre_changes(2) + centre_changes(3)
          chain1 = [chain1, idx(j)];
          chain2 = [chain2, idx(j+1)];
        else
          chain1 = [chain1, idx(j+1)];
          chain2 = [chain2, idx(j)];
        end
        j = j + 1;
      elseif on_off(1,idx(j)) <  on_off(3, chain2(end)) - 5 ...
          && on_off(1,idx(j)) >= on_off(3, chain1(end))
        chain1 = [chain1, idx(j)];
      elseif on_off(1,idx(j)) <  on_off(3, chain1(end)) - 5 ...
          && on_off(1,idx(j)) >= on_off(3, chain2(end))
        chain2 = [chain2, idx(j)];
      else
        % Match to the one closer in time-of-day
        centre_changes = centres (idx(j)) - centres([chain1(end), chain2(end)]);
        centre_changes(centre_changes > slots_per_day/2) ...
          = centre_changes(centre_changes > slots_per_day/2) - slots_per_day;
        centre_changes = abs (centre_changes);
        if centre_changes(1) < centre_changes(2)
          chain1 = [chain1, idx(j)];
        else
          chain2 = [chain2, idx(j)];
        end
        figure(101); show_on_off(on_off(:,chain1), randn(48, 365));
        figure(102); show_on_off(on_off(:,chain2), randn(48, 365));
        %keyboard;
      end
    end
    if j == penultimate
      j = j + 1;
      centre_changes = centres (idx(j)) - centres([chain1(end), chain2(end)]);
      centre_changes(centre_changes > slots_per_day/2) ...
        = centre_changes(centre_changes > slots_per_day/2) - slots_per_day;
      centre_changes = abs (centre_changes);
      if centre_changes(1) < centre_changes(2)
        chain1 = [chain1, idx(j)];
      else
        chain2 = [chain2, idx(j)];
      end
    end
    if length (chain2) > length (chain1)
      [chain1, chain2] = deal (chain2, chain1);
    end
    
    % If there are really more than two chains, skip ot avoid confusion.
    if any (on_off(1,chain1(2:end)) < mean (on_off([1,3], chain1(1:end-1))))
      continue;
    end

  else
    [~, chain1] = sort (on_off(1,:));
    chain2 = [];
  end
  figure(101); show_on_off(on_off(:,chain1), zeros(48, 365));
  figure(102); show_on_off(on_off(:,chain2), zeros(48, 365));
  figure(103); show_on_off(on_off, zeros(48, 365));
  %keyboard;
  
  % Calculate per-transition quantities
  duration_ratio = [duration(chain1(2:end))./duration(chain1(1:end-1)), ...
                    duration(chain2(2:end))./duration(chain2(1:end-1))];
  duration_ratio(duration_ratio < 1) = 1 ./ duration_ratio(duration_ratio < 1);
  list_duration_ratio(i) = { duration_ratio };
  % Duration ratio pdf: (1/3) delta(r-1) + (2/3) exp(-3.25 (r-1))
  
  centre_changes = [centres(chain1(2:end))-centres(chain1(1:end-1)), ...
                    centres(chain2(2:end))-centres(chain2(1:end-1))];
  gaps = [(on_off(1,chain1(2:end)) > on_off(3, chain1(1:end-1))), ...
          (on_off(1,chain2(2:end)) > on_off(3, chain2(1:end-1)))];
  centre_changes = mod (centre_changes, slots_per_day);
  backwards = (centre_changes > slots_per_day/2);
  centre_changes(backwards) = centre_changes(backwards) - slots_per_day;
  list_change_centres_no_gap(i) = { centre_changes(~gaps) };
  list_change_centres_gap(i)    = { centre_changes(gaps)  };
  
  list_day_gaps1(i) = { [on_off(1,chain1(2:end))-on_off(3,chain1(1:end-1))] };
  list_day_gaps2(i) = { [on_off(1,chain2(2:end))-on_off(3,chain2(1:end-1))] };
  
  list_runlengths(i) = { [on_off(3,chain1)-on_off(1,chain1), ...
                          on_off(3,chain2)-on_off(1,chain2)] };

end

ldr = sort ([list_duration_ratio{:}]);
plot (log (ldr), log (1 - (1:length (ldr)) / length (ldr)), [0, 0.03], [0, -0.39], [0, 1.5], [-0.3, -4.4]);
iii = (ldr > 1.03);
plot (log(ldr(iii)), log (1 - (1:sum (iii)) / sum(iii)), [0.02, 2.5], [0,-6.5]);
keyboard;

iii = (ldr > 1.03);
plot (log(ldr(~iii)), log (1 - (1:sum (~iii)) / sum(~iii)), [0, log(1.03)], [0, -2]);
keyboard

% hldr = 1 ./ (ldr(9:end) - ldr (1:end-8));
% iii = (ldr(5:end-4) > 1.03);
% plot (ldr([false(1,4), iii]), log (1 - cumsum(hldr(iii)) / sum(hldr(iii))), [1, 3.7], [0, -12]);
% keyboard;

lcc_n = sort ([list_change_centres_no_gap{:}]);
hlcc_n = 1 ./ (lcc_n(9:end) - lcc_n (1:end-8));
plot ((lcc_n(5:end-4) - lcc_n(end-4:-1:5))/2, log (hlcc_n + hlcc_n(end:-1:1)));
keyboard;

alcc_n = sort (abs ([list_change_centres_no_gap{:}]));
halcc_n = 1 ./ (alcc_n(9:end) - alcc_n(1:end-8));
plot (alcc_n(5:end-4), log (halcc_n));
plot (alcc_n(9:end), log (halcc_n));

alcc_n = sort (abs ([list_change_centres_no_gap{:}]));
halcc_n = 1 ./ (alcc_n(19:end) - alcc_n(1:end-18));
plot (log (1 + alcc_n(19:end)), log (halcc_n));

lcc_g = sort ([list_change_centres_gap{:}]);
hlcc_g = 1 ./ (lcc_g(9:end) - lcc_g (1:end-8));
plot ((lcc_g(5:end-4) - lcc_g(end-4:-1:5))/2, log (hlcc_g + hlcc_g(end:-1:1)));
keyboard;

alcc_g = sort (abs ([list_change_centres_gap{:}]));
halcc_g = 1 ./ (alcc_g(9:end) - alcc_g(1:end-8));
plot (alcc_g(5:end-4), log (halcc_g));
plot (alcc_g(9:end), log (halcc_g));

alcc_g = sort (abs ([list_change_centres_gap{:}]));
halcc_g = 1 ./ (alcc_g(19:end) - alcc_g(1:end-18));
plot (log (1 + alcc_g(19:end)), log (halcc_g));

ldg1 = sort ([list_day_gaps1{:}]);
ldg1 = ldg1(ldg1 >= 0);
plot (log (1+ldg1), sqrt (-log (1 - (1:length (ldg1)) / length (ldg1))), [0, 5.5], [0.65, 2.5]);  
keyboard;

ldg2 = sort ([list_day_gaps2{:}]);
ldg2 = ldg2(ldg2 >= 0);
plot (log (1+ldg2), sqrt (-log (1 - (1:length (ldg2)) / length (ldg2))), [0, 5], [1, 1.9]);
keyboard;

% hldg = 1 ./ (ldg(9:end) - ldg(1:end-8));
% iii = (ldg(5:end-4) > 0);
% plot (ldg([false(1,4), iii]), log (1 - cumsum(hldg(iii)) / sum (hldg(iii))));
% keyboard;
