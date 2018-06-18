function [pump, cor_pool, power, profile, pattern, timer, cor_am, rect_am, rectangles] = has_pool_pump (userData, vampires, meta, cust, solar_cor)
  %record = [];
  issolar = any (solar_cor ~= 0);
  pump = 0;
  cor_pool = zeros(size(userData));
  power = 0;
  profile = 0;
  pattern = 0;
  timer = [];
  rectangles = [];

%%%
% TODO:
% - discount "reliability" of runs near wake-up time, especially winter
% - select which rectangles are "HWS" (or at least grouped, using AIC)
% - if reliability high, coerce match_jump to select an existing run
% - Align rectangles.  How?
% - Confidence of rectangles.  How?
% - day-of-week
% - Adaptive: if "spike then flat" found in some runs, look for it elsewhere
%             if day-of-week found in some runs, look for it elsewhere
%             if overnight found in some runs, look for it elsewhere
%             if morning/afternoon operation found, look for it elsewhere
%             if miss one block, look for it in other blocks of the same day
% - fill gaps

% Terms:
%   run -- consecutive days with a given jump
%   reliability -- confidence that a run is from a timed device
%   rectangle -- consecutive days with a given turn-on and matching turn-off
%   quality/trust -- confidence that a rectangle is actually a timed device
%   quality -- confidence that the end of a rectangle is the true end?

% Columns in "runs"
% 1. start day  (first day of run)
% 2. stop day   (first day after run)
% 3. half-hour slot
% 4. half-hour slot, interpolated
% 5. jump (increase in power)
% 6. jump, as estimated from interpolated time
% 7. reliability
% 8. trend (average increase based on previous and subsequent slots)


% Final data structure:
% - List of timer settings, with start day+time, end day+time, on/off/on/off
% - List of resulting "rectangles"

% "rectangle" data structure:
% - on_off  -- quadruple of (on-day, on-hour, off-day, off-hour)
% - power   -- estimated power of device
% - edges   -- triple of estimated jump up, power, jump down
% - missed  -- list of days that rectangle device is not on during rectangle
% - burst   -- cached set of 30 min slots constituting this rectangle
% - trust   -- numeric estimate of how likely this is to be a pump not noise
% After estimating timer settings, rectangles may also get
%               TODO refactor: Does this need adjusting?
% - alt_days -- days-of-week (1/Jan being 0) for which timer settings differ
% - alt_hrs  -- first half-hour and 1+last half-hour of new time.

% Fields in rectangles.trust
% 1. top/bottom Edge quality
% 2. Body quality
% 3. Miss quality
% 4. Do we start when another rectangle stops?
% 5. How many other rectangles are on at the same time of day?
% 6. Q_p -- adherence to common power
% 7. Q_L -- left edge
% 8. Q_R -- right edge

%known_pools = { '06347', '12202', '14366', '23024', '90004', '90011', '90027', '90028', '90046', '90047', '90052', '90055', '90057', '90059', '90063', '90064', '90069', '90070', '90072', '90074', '90084', '90087', '90089', '90091', '90095', '90096', '90097', '90115', '90116', '90119', '90128', '90131', '90144', '90164', '90171', '90173' };
%known_pools = {'814', '1501', '1858', '2083', '2105', '2452', '2620', '2903', '573'};
 known_pools = {};

  c = stripSpikes (squeeze(userData)');

  valid_days = find(~isnan(c(1,:)));
  cc = c (:,valid_days);

  if isempty(cc)
    pump = 0;
  if any (strcmp (cust, known_pools)), keyboard ('No valid days> '); end
    return
  end

  vampires = vampires(valid_days);
  cv = cc - repmat(vampires,  [size(cc,1), 1]);

  data = squeeze (evalin ('caller', 'data(i,:,:)'))';
  % FIXME: Maybe pass raw data with spikes, since AM is often spikey.
  if 0
    [cor_am, rect_am] = morning (c, meta, valid_days, cust, data);
  else
    cor_am = [];
    rect_am = [];
  end
%return

  % Try to guess when the pump is on.
  % Some pumps have two (or more?) powers, so record time and value

  [runs, d, dd, band] = runs_from_cc (cc, meta);

  summer = meta.summer;
  if (length(valid_days) < size(c,2))
    [~, summer] = ismember (summer, valid_days);
    summer = summer(summer~=0);
  end

  if isempty(runs)
    on_all_day = full(vampires > 0.3);
    pump = (sum(on_all_day(summer)) > 0.5 * sum(~isnan(c(1,meta.summer))));
 if any (strcmp (cust, known_pools))
  if pump
    keyboard ('always> ');
  else
    keyboard ('no pump> ');
  end
 end
    return
  end

  [~, t_idx] = sort(runs(:,7));
  max_reliability  = runs(t_idx(end),7);

  % if no reliable runs detected, and there isn't a high base
  % level of load (maybe from an always-on pool pump), say "no pump"
  %       (always-on pumps should be detected in has_untimed_pump)
  if max_reliability < 15 && median(vampires(~isnan(vampires))) < 0.3
  if any (strcmp (cust, known_pools)), keyboard ('no reliable runs> '); end
    return
  end


  % If one looks like a heater turned on each winter morning, and
  % not much else, skip expensive processing
  main = (runs(:,7) > 15); %20);
  if ~any(main)
      main = t_idx(end);      % if no reliable runs, take the best
  end
  tr = @(x)(x * meta.SamPerDay/48);
  if isequal (meta.hemisphere, 'south')
    is_winter = runs(main,1) > 81 &  runs(main,2) < 281;
  else
    is_winter = runs(main, 1) > 281 | runs(main, 2) < 81;
  end
  if  max(runs(main,5)) > 0 && ...
    all(   runs(main,3) < tr(20) ...
        & (runs(main,3) < tr(18) | runs(main,5) < 0) ...
        &  runs(main,3) > tr(8) ...
        &  is_winter)
    pump = -1;
if any (strcmp (cust, known_pools)), keyboard ('heater> '); end
    return
  end

  % For testing, ignore electric hot water
  % That should be stripped out before this is called.
  if all((runs(main,3) == meta.SamPerDay-1 & runs(main,5) > 0) | (runs(main,3) <= 2 & runs(main,5) < 0)) && ismember(meta.SamPerDay-1, runs(main,3))
    pump = -2;
if any (strcmp (cust, known_pools)), keyboard ('HWS> '); end
    return
  end

  % Save for reference when splitting dubious rectangles
  orig_runs = runs;

  % Merge runs oscillating between neighbouring times
  runs = merge_oscillating_runs (runs, cc, dd, band);

  % Merge adjacent runs if that seems suitable
  runs = merge_adjacent_runs(runs, cc, dd, band);

  % Remove runs whose step size is far from the weighted mean step size
  mid = sum(abs(runs(:,6)).*runs(:,7)) / sum(runs(:,7));
  idx = (abs(log(abs(runs(:,6)) / mid)) < 0.4);
  mid = sum(abs(runs(idx,6)).*runs(idx,7)) / sum(runs(idx,7));

  runs = runs ((abs(log(abs(runs(:,6)) / mid)) < 0.7) | (runs(:,7) > 5), :);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % create rectangles based on well-trusted runs
  for pass = 1:2
    [rectangles, new_runs, cancelled] = runs_to_rectangles (runs, cv, dd, band);
    cor = cor_pool;
    cor(valid_days,:) = cancelled;
    runCount = size(runs,1);

    % Check if the pump caused errors in estimating the solar generation
    if issolar
      on_off = [rectangles.on_off];
      min_ov = 80;       % minimum overlap
      long = find(orig_runs(:,2) - orig_runs(:,1) > min_ov);
      long = long(orig_runs(long,7) > 10);
      ons  = orig_runs(long(orig_runs(long, 3) <= size (cv,1)/2 & orig_runs(long, 5) > 0),:);
      offs = orig_runs(long(orig_runs(long, 3) >= size (cv,1)/2 & orig_runs(long, 5) < 0),:);
      corr = [];
      for i = 1:size (ons, 1)
        off = offs(offs(:,2) > ons(i,1) + min_ov & offs(:,1) < ons(i,2) - min_ov,:);
        off = off(abs (off(:,5)) > 0.7 * ons(i,5) & abs (off(:,5)) < 1.4 * ons(i,5),:);
        if size (off, 1) > 1
          % TODO: eliminate overlap
        end
        for j = 1:size (off, 1)
          st = max (ons(i,1), off(j,1));
          en = min (ons(i,2), off(j,2));
          % check this rectangle doesn't already exist
          match_rect = on_off(1,:) < st + min_ov/10 ...
                       & on_off(3,:) > en - min_ov/10 ...
                       & abs (on_off(2,:) - ons(i,4)) <= 1 ...
                       & abs (on_off(4,:) - off(j,4)) <= 1;
          if any (match_rect)
            continue;
          end
          % check no big runs between them
          runs_between = orig_runs(:, 1) < en - min_ov/2 ...
                         & orig_runs(:, 2) > st + min_ov/2 ...
                         & orig_runs(:, 2) - orig_runs(:,1) > min_ov/2 ...
                         & abs (orig_runs(:, 5)) > 0.7 * ons(i,5) ...
                         & orig_runs(:, 7) > 10 ...
                         & orig_runs(:, 3) > ons(i,3) ...
                         & orig_runs(:, 3) < off(j,3);
          if any (runs_between)
            continue;
          end
          corr = solar_correction (cv(:, st:en-1), ons(i, 4), offs(j, 4));
        end
      end
      if isempty (corr) && ~isempty (on_off)
        midday = find (on_off(2,:) < 0.45 * size(cv, 1) ...
                    & on_off(4,:) > 0.55 * size(cv, 1) ...
                    & on_off(3,:) > on_off(1,:) + 0.05 * (40 + size(cv, 2)));
        if ~isempty (midday)
          % search for dip in any rectangle, starting from the longest
          [~, idx] = sort(on_off(1, midday) - on_off(3, midday));
          for j = idx
            i = midday (j);
            corr = solar_correction (cv(:, on_off(1,i):on_off(3, i)-1), ...
                                     on_off(2,i), on_off(4,i));
            if ~isempty (corr)
              break;
            end
          end
        end
      end
      if isempty (corr)
        break;
      end
      cv = bsxfun (@plus, cv, corr);
    else
      break;
    end
  end

  if isempty(rectangles)
 if any (strcmp (cust, known_pools))
 tmp = zeros(size(c));
 for j = 1:runCount
  if runs(j,3) ~= -1
   tmp(runs(j,3), valid_days(runs(j,1)+1:runs(j,2))) = runs(j,7)*sign(runs(j,5));
  end
 end
 figure(2); imagesc (tmp);
 keyboard ('no rectangles> ');
 end
    return
  end

  % Attributes of a "rectangle":
  %  Start/end date
  %  Start/end hour
  %  Power level
  %  Exception:  Days of week
  %              Start/end hour
  %  Missed days
  %  override-to-on days
  %  Device number (e.g., pool pump, heater, ...)
  %
  %  Meta: Reliability of each edge
  %        Reliability of power level
  %        Reliability of device number


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Abort now if unlikely to be a pool
  [pump, pwr, cor] = seems_like_pool (rectangles, cor, c, meta, valid_days, 1);
  if pump == 0
    return
  end
  pump = 0;   % For early abort returns

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Merge adjacent rectangles
  on_off = [rectangles(:).on_off];
  for r = 1:length(rectangles)
    prev = abs(on_off(2,:)-on_off(2,r)) < 0.5 & ...
           abs(on_off(4,:)-on_off(4,r)) < 0.5 & ...
           abs(on_off(3,:)-on_off(1,r)) < 0.5;
    if any(prev)
      prev = find(prev);
      prev = prev(1);
      % Should take weighted means. Will be fixed by get_weighted_edge
      rectangles(r).on_off = [on_off(1,prev)
                         mean(on_off(2,[r,prev]))
                              on_off(3,r)
                         mean(on_off(4,[r,prev]))];
      on_off(:,r) = rectangles(r).on_off;
      rectangles(r).power = mean([rectangles([r,prev]).power]);
      rectangles(prev).power = 0;
    end
  end
  rectangles = rectangles([rectangles.power] ~= 0);

  % Match up rectangles, trying to identify days at which the timer changed
  rectangles = align_rectangles (rectangles, cv, d);

  if isempty(rectangles)
    return
  end

  rectangles = assign_trust (rectangles, cv, pwr, issolar, false);
  if isempty (rectangles)
    return
  end

  [~, idx] = sort(-skipNaN(@sum,[rectangles(:).trust]));
  rectangles = rectangles(idx);
  %on_off = round([rectangles.on_off]);

  rectangles = remove_noise_rectangles (rectangles, size(cv,2));
  if isempty (rectangles)
    return
  end

  rectangles = rectangles(any([rectangles(:).trust],1));
  cor = set_cor_pool (rectangles, cor, valid_days);

  % Find actual on-time, not just the nearest hour
  % Find days pump not turned on at all
  % Find days pump left on all day
  % Find days pump is on at an unusual time
  % Find days pump is on, but no matching "run" was found
  % Find "actual" power(s?), and remove noise rectangles

  pump = seems_like_pool (rectangles, cor, c, meta, valid_days, 2);
  if pump == 0
    return
  end

  pwr = common_power(rectangles);
  % [pwr, other_rectangles] = common_power (rectangles);  % ??
  [~, rectangles] = timer_settings(rectangles, cv, issolar);
  % figure(10); show_rectangles (rectangles, cv);
  rectangles = extend_rectangles  (rectangles, cv, pwr, issolar, valid_days);
  %figure(12); show_rectangles (rectangles, cv, []);
  if isempty (rectangles)
    return;
  end

  on_off = [rectangles.on_off];
  rectangles = rectangles(on_off(3,:) - on_off(1,:) > 1 ...
                          & abs (diff (on_off([2,4],:))) > 1);
  rectangles = tweak_on_off_pwr   (rectangles, cv, pwr, issolar, valid_days);
  rectangles = rectangles ([rectangles.power] > 0.1);
  if isempty (rectangles)
    return;
  end
  rectangles = merge_neighbour_rects (rectangles, cv);
%   rectangles = assign_trust (rectangles, cv, pwr, issolar, false);
  rectangles = split_rectangle_days (rectangles, cv, orig_runs, issolar, valid_days);
  rectangles = check_missing_days (rectangles, cv, pwr, issolar, valid_days);
%  rectangles = extend_rectangles (rectangles, cv, pwr, issolar, valid_days);
  on_off = [rectangles.on_off];
  rectangles = rectangles(on_off(3,:) - on_off(1,:) > 1 ...
                          & abs (diff (on_off([2,4],:))) > 1);
  if isempty (rectangles)
    return
  end

  % Very preliminary work to find what is going on in missed days.
  % TODO: cull spurious rectangles from output.
  [rectangles, others] = set_missed_precise (rectangles, mean ([rectangles.power]), ...
                                   cv, cc, orig_runs, issolar);
  rectangles = [rectangles, others];

  [timer, rectangles] = timer_settings(rectangles, cv, issolar);
badness = [rectangles.badness];
figure(6); show_heatmap (1./sqrt(badness(1,:)), rectangles, cv);
figure(7); show_heatmap (1./sqrt(badness(2,:)), rectangles, cv);

  % Trim noise that is low power and has "bad" edges
  badness = [rectangles.badness];
  [sorted_bad, idx1] = sort (1./sqrt(badness(1,:)) + 1./sqrt(badness(2,:)));
  on_off = [rectangles.on_off];
  weights = (on_off(3,:) - on_off(1,:)) .* mod1 (on_off(4,:) - on_off(2,:), size (cv, 1));
  idx2 = [rectangles.power] < 0.6 * sum (weights .* [rectangles.power]) / sum (weights);
  idx2(idx1(end/2:end)) = false;
  idx2(idx1(sorted_bad > mean (sorted_bad) / 2));
  rectangles = rectangles(~idx2);
  if isempty (rectangles)
    return
  end
  figure(1); show_rectangles (rectangles, cv, []);

  if isfield (rectangles(1), 'alt_days')
    % timer_settings has identified day-of-week variations.
    % Now we can get much better estimates of jump sizes and missed days
    cor = set_cor_pool_precise (rectangles, cor, valid_days);
    rectangles = powers_reliabilities_DoW (rectangles, cv, cor);
    rectangles = set_missed_precise (rectangles, pwr, cv, cc, orig_runs, issolar);
  end

  % Get rid of overlap
  % Extend again  (loop?)



pdfs.duration_ratio = @duration_ratio_pdf;
pdfs.day_gap = @day_gap_pdf;
pdfs.centre_change_gap = @centre_change_gap_pdf;
pdfs.centre_change_no_gap= @centre_change_no_gap_pdf;
% Re-assign trust to rectangles based on Kolmogorov-Smirnov test.
rectangles = SNR_trust (rectangles, cv);
if isempty (rectangles)
  return
end

[chains, dists, chain_values] = timed_chains (rectangles, pdfs, cv);
qualities = scalar_quality (rectangles, chains, dists, chain_values, cv, issolar);

figure(3); show_rectangles (rectangles(qualities > 0.25 * max (qualities)), cv);
figure(4); show_rectangles (rectangles(qualities < 0.25 * max (qualities)), cv);


  cor = set_cor_pool_precise (rectangles, cor, valid_days);
  figure(2); imagesc (cor');

  [pump, power, cor_pool] = seems_like_pool (rectangles, cor, c, meta, valid_days, 3);
  if nargout >= 3
    profile = -sum(cor_pool, 1)';
    if nargout >= 5
      pattern = 1;
    end
  elseif pump == 0
    profile = 0;
    pattern = 0;
  end
  if pump
    rect_dir = '~/rsrch/power/UnitedEnergy/pools';
    if (~exist (rect_dir, 'file'))
      mkdir (rect_dir);
    end
    rect_file = [rect_dir, '/rect_', cust, '.mat'];
    if exist (rect_file, 'file')
      try
        load (rect_file);
        rectangles = r;
      catch
      end
    end
%    [t, r, ok] = pool_ground_truth (cv, timer, rectangles);
    ok = false | false;
    if ok
      save ([meta.metaDataFile, rect_file], 'r', 't');
    end
  end
  if isfield (rectangles, 'alt_hrs')
    %figure(1); imagesc (userData');
    %figure(2); imagesc (-cor');
    %fprintf('pump %g  power %g\n', pump, power(1));
    %keyboard
  end
  figure(9); imagesc (cv);
end

function [runs, d, dd, band ] = runs_from_cc (cc, meta)
  % Extract time/start-day/end-day runs from 2-D array of jump sizes
  d1 = [0; diff(cc(:))];
  d1(1) = d1(1+meta.SamPerDay*(size (cc,2)>1));    % guess first jump
  d  = reshape (d1, size (cc));      % half-hour differences
          % create a band around midnight for identifying steps
  band = 2;
  dd = [d; d(1:2*band, :)];
  variance_by_hour = var (diff (cc'));
  med_lengths = [5, 7, 9, 11, 13, 15, 17, 19, 21];
  length_bins = min (1+floor (variance_by_hour*5), length (med_lengths));
  smoothed_c = zeros (size (cc));
  for i = 1:length (med_lengths)
    idx = (length_bins == i);
    smoothed_c(idx,:) = -rolling_min (-rolling_min (cc(idx,:)', med_lengths(i)), med_lengths(i))';
  end
  [~, times, jp] = find_jumps (smoothed_c, 0.05, 1);

  sp = sparse (mod (round (times),size (cc,1))+1, ...
               floor (round (times)/size (cc,1))+1, ...
               jp, ...
               size (smoothed_c,1), size (smoothed_c,2));
  sp = sparse (medfilt1 (full (sp),9,[],2));
  runs = runs_from_sparse (sp, dd, cc, band);
end

function [hcor, rectangles] = morning (c, meta, valid_days, cust, uncorrected)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  hcor = zeros (size (c'));
  rectangles = [];

  c1 = c(floor(6/48*meta.SamPerDay):floor(18/48*meta.SamPerDay), :);
  c2 = c(floor(8/48*meta.SamPerDay):floor(20/48*meta.SamPerDay), :);
  c_dls(:,  meta.daylight_saving) = c1(:,  meta.daylight_saving);
  c_dls(:, ~meta.daylight_saving) = c2(:, ~meta.daylight_saving);

  [v_wk, vd_idx] = intersect (valid_days, meta.weekdays);
  c_wk = c_dls (:, v_wk);
  if isempty (c_wk)
    return
  end
  [hruns, ~, hdd, hband] = runs_from_cc (c_wk, meta);
  hruns = hruns(hruns(:,3) ~= 1, :);  % skip phantoms due to truncated day
%  global debug
%  if debug.flags & 1
%    save (debug.fname, 'hruns');
%  end

  if isempty (hruns)
    return
  end
  % Merge runs oscillating between neighbouring times
  hruns = merge_oscillating_runs (hruns, c_wk, hdd, hband);

  % Merge adjacent runs if that seems suitable
  hruns = merge_adjacent_runs(hruns, c_wk, hdd, hband);

  %?? Map back to normal coordinates so missing weekends are found
  %   by runs_to_rectangles?  Or do r_to_r only for weekdays, and
  %   check weekends separately?
  %hruns(:,1:2) = v_wk(hruns(:,1:2));
  %hruns(:,3:4) = hruns(:,3:4) + 7;
  rectangles = runs_to_rectangles (hruns, c_wk, hdd, hband);

  % Questions:
  % - Split rectangles spanning daylight saving.

  if isempty (rectangles)
    return
  end
  % Map on_off to index into valid days
  for i = 1:length (rectangles)
    r = rectangles(i);
    on_off = r.on_off;
    rectangles(i).on_off(1) = vd_idx(on_off(1));
    rectangles(i).on_off(3) = vd_idx(on_off(3)-1)+1;
  end

  dls_changes = find (diff (meta.daylight_saving));
  dls(1) = find (valid_days >= dls_changes(1), 1);
  dls(2) = find (valid_days >= dls_changes(2), 1);

  % Split rectangles that span the daylight savings transitions
  for n = 1:2
    oo = [rectangles(:).on_off];
    spn = (oo(1,:) < dls(n)) & (oo(3,:) > dls(n));
    if any (spn)
      e = length (rectangles);
      spn = find (spn);
      for j = 1:length (spn)
        i = spn(j);
        new_rect = rectangles(i);

        % split rectangle, ensuring no "missed" days at the ends
        rectangles(i).on_off(3) = dls(n);
        missed = intersect (rectangles(i).missed, oo(1,i):dls(n)-1);
        if ~isempty (missed) && missed(end) == dls(n)-1
          offset = find (diff (missed(end:-1:1)) - 1, 1) - 1;
          if ~isempty (offset)
            rectangles(i).on_off(3) = rectangles(i).on_off(3) - offset;
          end
          missed = missed(1:end-offset);
        end
        rectangles(i).missed = missed;

        new_rect.on_off(1) = dls(n);
        missed = intersect (new_rect.missed, dls(n):oo(3,i)-1);
        if ~isempty (missed) && missed(1) == dls(n)
          offset = find (diff (missed)-1, 1);
          if ~isempty (offset)
            new_rect.on_off(1) = new_rect.on_off(1) + offset - 1;
          end
          missed = missed(offset:end);
        end
        new_rect.missed = missed;
        rectangles(e+j) = new_rect;
      end
    end
  end

  % Shift to "normal" times
  oo = [rectangles(:).on_off];
  dls_rect =  (oo(3,:) < dls(1) ...
        | oo(1,:) > dls(2));
  oo(:, dls_rect) = bsxfun (@plus, oo(:, dls_rect), [0;5;0;5]);
  oo(:,~dls_rect) = bsxfun (@plus, oo(:,~dls_rect), [0;7;0;7]);
  for i = 1:length (rectangles)
    rectangles(i).on_off = oo(:,i);
  end

  % Delete rectangles that are only on on days others are too
  % Delete artefacts caused by "wrapping" on the time window
  before = bsxfun (@gt, oo(1,:), oo(1,:)');
  after  = bsxfun (@lt, oo(3,:), oo(3,:)');
  [~, covered] = find (before & after);
  split = (oo(2,:) > oo(4,:));
  rectangles(unique ([covered; find(split)'])) = [];
  if isempty (rectangles)
    return
  end

  % If two rectangles overlap, delete the smaller
  % TODO: should use quality, not size
  % TODO: consider trimming instead of deleting
  oo = [rectangles(:).on_off];
  overlap = bsxfun (@gt, oo(3,:), oo(1,:)') ...
          & bsxfun (@lt, oo(1,:), oo(3,:)') ...
          & bsxfun (@gt, oo(4,:), oo(2,:)') ...
          & bsxfun (@lt, oo(2,:), oo(4,:)');
  % clear diagonals -- everything overlaps with itself
  overlap (1:(length (overlap)+1):numel(overlap)) = 0;
  [a, b] = find (overlap);
  if any (a)
    sizes = (oo(3,:) - oo(1,:)) .* (oo(4,:) - oo(2,:));
    [~, idx] = sort (sizes(a), 'descend');
    for i = idx(:)'
      if sizes(a(i)) < sizes(b(i))
        sizes(a(i)) = 0;
      else
        sizes(b(i)) = 0;
      end
    end
    rectangles(sizes == 0) = [];
  end

  % - What happens on weekends?  Same, nothing, other?
  cv = c(:, valid_days);
  vd_we = 1:length (valid_days);
  vd_we (vd_idx) = 0;
  vd_we = vd_we (vd_we > 0);    % weekend days within valid_days
  for i = 1:length (rectangles)
    r = rectangles(i);
    on_off = r.on_off;
    rectangles(i).missed = vd_idx(r.missed);
    v_we = vd_we (vd_we >= on_off(1) & vd_we < on_off (3));
      if isempty (v_we)               % no weekends?  Treat as non-uniform
        rectangles(i).profile = 3;
        continue
      end
    v_week = vd_idx (vd_idx >= on_off(1) & vd_idx < on_off (3));
    burst_plus = round(on_off(2))-1:round (on_off(4));
    % Classify:
    % Uniform profile? (sum of deviations from median)
    % Scaled profile? (normalize peak, then sum of deviations from med)
    % Flat profile?
    %
    % If uniform profile, weekends different if profile differs
    % If non-uniform profile, weekends different if profile much smaller

    shifted = bsxfun (@minus, cv(burst_plus, v_week), mean (cv(burst_plus, v_week), 1));
    med_wk = median (shifted, 2);
    dev_wk = abs (bsxfun (@minus, shifted, med_wk));
    if mean (dev_wk(:)) < 0.7 * mean (abs (med_wk - median (med_wk)))
      rectangles(i).profile = 0;
      % uniform profile
      % Look for non-flat profile (profile 1)?
      %fprintf('uniform profile');
    else
      % look for profile that is uniform up to scaling
      shifted = bsxfun (@rdivide, shifted, max (shifted));
      shifted(isnan (shifted)) = 0;
      med_wk = median (shifted, 2);
      dev_wk = abs (bsxfun (@minus, shifted, med_wk));
      if mean (dev_wk(:)) < 0.6 * mean (abs (med_wk - median (med_wk)))
        rectangles(i).profile = 2;
        % scaled uniform
        %fprintf('scaled profile');
      else
        rectangles(i).profile = 3;
        % non-uniform profile
        %fprintf('non-uniform profile');
      end
    end

    shifted = bsxfun (@minus, cv(burst_plus, v_we), mean (cv(burst_plus, v_we), 1));
    med_we = median (shifted, 2);
    %med_we = med_we - min (med_we([1,end]));

    if rectangles(i).profile < 2
      df = norm (med_wk - med_we) / (max ([med_wk; med_we]) - min ([med_wk; med_we]));
    else
      df = (mean (med_wk) - mean (med_we)) / mean ([med_wk; med_we]);
    end
    if df > 0.5
      rectangles(i).alt_days = meta.weekends(1:2)';
      % Dummy!!  Should look for actual times on weekends
      rectangles(i).alt_hrs = [on_off(2), on_off(2)+1e-4];
      rectangles(i).profile = rectangles(i).profile + 4;
    end
  end

  % Is this actually a pool pump?  (Duration?)
  % - Number of days
  % - Quality of rectangles.
  oo = [rectangles.on_off];
  sizes = oo(3,:) - oo(1,:) .* (oo(4,:) - oo(2,:));
  power = (sizes * [rectangles.power]') / sum (sizes);
  if power == 0
    return
  end
  issolar = false;  % FIXME -- should be passed in
  rectangles = end_days (rectangles, cv);
  rectangles = assign_trust (rectangles, cv, power, issolar, true);
  if isempty (rectangles)
    return
  end

  oo = [rectangles.on_off];
  lengths = oo(3,:) - oo(1,:);
  width_ratios = lengths ./ (oo(4,:) - oo(2,:));

  if sum (lengths) < 0.3 * length (valid_days)
    t = [rectangles.trust];
    r = (width_ratios > 10 & t(1,:) > 50);
    if ~any (r)
      %fprintf ('My guess: spurious\n');
      return
    end
  end

  % Estimate correction, not making residual the min for the  day
  hcor = set_cor_pool_precise (rectangles, hcor, valid_days);

  % Fudge for over-estimated power consumption.
  % Should also check/fix power consumption calculation
  corrected = bsxfun (@max, c + hcor', max (min (c), 0));
  hcor = (corrected - c)';

if 1
  figure(4); imagesc (hcor');
  figure(5); imagesc (c + hcor');
  fprintf('\ncust = %s', cust);
  figure(6); imagesc (uncorrected);
  figure(7); plot (uncorrected(:, meta.July));
  figure(8); plot (uncorrected(:, meta.July) + hcor(meta.July,:)');
  figure(9); show_rectangles (rectangles, cv);

  rectangles = extend_rectangles (rectangles, cv, [], valid_days);
%   [~, ~, ok] = pool_ground_truth (cv, [], rectangles);
%   if (~ok)
%    keyboard
%   end
%  for i = 1:length (rectangles)
%    fprintf ('proflle %g\n', rectangles(i).profile);
%    figure(10);
%    %plot ()
%  end
end
if 1+1 == 0
  show_runs (hruns, c_wk, valid_days); % avoid "unused" error
end
%    if ~isempty(rectangles)
%      weight = zeros (1,19);
%      for i = 9:19
%        a = sum (hruns(hruns(:,3) == i & hruns(:,5) > 0, 7));
%        if ~isempty (a)
%          weights(i) = a;
%        end
%      end
%      on_time = 1:19 * weight / sum (weight);
%
%      main = (hruns(:,7) > 20);
%      if isempty (main)
%        main = (hruns(:,7) > 10);
%        if isempty (main)
%          main = 1:size(hruns,1);
%        end
%      end
%      st = min (hruns(main, 1));
%      en = max (hruns(main, 2));
%      figure(4); plot (c_wk(:, st+1:en));
%      figure(5); plot (median (c_wk(:, st+1:en), 2));
%      figure(6); imagesc (hcor');
%      show_runs (hruns, c_wk, valid_days);
%    end
end

function show_runs (runs, c, valid_days)
 tmp = zeros(size(c));
 for j = 1:size(runs,1)
  if runs(j,3) ~= -1
   tmp(runs(j,3), valid_days(runs(j,1):runs(j,2)-1)) = runs(j,7)*sign(runs(j,5));
  end
 end
 figure(2); imagesc (tmp);
 figure(3); plot (c);
end

function rectangles = end_days (rectangles, c)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Tweak start/stop days of rectangles, and assign a confidence to each.
  g = 4;
  on_off = [rectangles.on_off];
        % approximate variance of all elements of a matrix, or NaN if empty
  var1d = @(x)(var ([x(:); mean(x(:))]));
  ep = 0.001;
  for r = 1:length(rectangles)
    power = rectangles(r).power;
    time_olap = on_off(4,:) > on_off(2,r) ...
              & on_off(2,:) < on_off(4,r);
    % Find times with no near-overlap of left edge
    o_lap = (on_off(3,:) > on_off(1,r) - g ...
           & on_off(1,:) < on_off(1,r) ...
           & time_olap);
    burst = rectangles(r).burst;
    if ~any (o_lap)
      mask = [zeros(length(burst),g), power*ones(length(burst),g)];
      %mask(1,g+1:2*g) = power * (on_off(2,r) - floor(on_off(2,r)));
      %mask(length(burst),g+1:2*g) = power * (1-(on_off(2,r) - floor(on_off(2,r))));
    else
      % only consider times that this differs from its neighbour
      burst = setdiff (burst, [rectangles(o_lap).burst]);
      mask = [zeros(length(burst),g), power*ones(length(burst),g)];
    end
    % Find quality of nominal edge, those feasible to the left, those to the
    %      right supported by extended edges.
    quality_L = zeros(1,size(c,2));
    if ~isempty(burst)
      for day = on_off(1,r)-1:-1:max(2, on_off(1,r)-1-g)
        c_range = max(1,1+day-g):min(day+g,on_off(3,r)-1);
        m_range = g-day + c_range;
        quality_L(day-1) = 1 / (ep + var1d(c(burst,c_range) - mask(:,m_range)));
        if any (c(burst,day-1) < 0.9 * power)
          break
        end
      end

      for day = max(2,on_off(1,r)):min(on_off(1,r)-1+g,on_off(3,r)-1)
        c_range = max(1,1+day-g):min(day+g,on_off(3,r)-1);
        m_range = g-day + c_range;
        quality_L(day-1) = 1 / (ep + var1d(c(burst,c_range) - mask(:,m_range)));
        if any (c(burst,day-1) < 0.9 * power)
          break
        end
      end
    end
    rectangles(r).quality_L = quality_L;

    % Repeat for right edge with image flipped.
    o_lap = (on_off(3,:) > on_off(3,r) ...
           & on_off(1,:) < on_off(3,r) + g ...
           & time_olap);
    burst = rectangles(r).burst;
    if ~any (o_lap)
      mask = [power*ones(length(burst),g), zeros(length(burst),g)];
      %mask(1,g+1:2*g) = power * (on_off(2,r) - floor(on_off(2,r)));
      %mask(length(burst),g+1:2*g) = power * (1-(on_off(2,r) - floor(on_off(2,r))));
    else
      % only consider times that this differs from its neighbour
      % If there are none, quality is equal for all day offsets
      burst = setdiff (burst, [rectangles(o_lap).burst]);
      mask = [power*ones(length(burst),g), zeros(length(burst),g)];
    end
    % Find quality of nominal edge, those feasible to the left, those to the
    %      right supported by extended edges.
    quality_R = zeros(1,size(c,2));
    if ~isempty(burst)
      for day = on_off(3,r):min(on_off(3,r)-1+g,size(c,2))
        c_range = max(on_off(1,r),1+day-g):min(day+g,size(c,2));
        m_range = g-day + c_range;
        quality_R(day-1) = 1 / (ep + var1d(c(burst,c_range) - mask(:,m_range)));
        %jump_R (day) = mean(c(burst,day) - c(burst,day+1));
        %quality_R (day) = 1 / (0.1 + var1d(c(burst,max(1,day-g):day)) ...
        %                           + var1d(c(burst,day+1:min(day+g, on_off(3,r)))));
        if day >= size(c,2) || any (c(burst,day+1) < 0.9 * power)
          break
        end
      end

      for day = on_off(3,r)-2:-1:max(on_off(1,r)+1, on_off(3,r)-1-g)
        c_range = 1+max(on_off(1,r)-1,1+day-g):min(day+g,size(c,2));
        m_range = g-day + c_range;
        quality_R(day-1) = 1 / (ep + var1d(c(burst,c_range) - mask(:,m_range)));
        %jump_R (day) = mean(c(burst,day) - c(burst,day+1));
        %quality_R (day) = 1 / (0.1 + var1d(c(burst,max(1,day-g):day)) ...
        %                           + var1d(c(burst,day+1:min(day+g, on_off(3,r)))));
        if any (c(burst,day+1) < 0.9 * power)
          break
        end
      end
    end
    rectangles(r).quality_R = quality_R;
  end
end

function [power, trust] = common_power(rectangles)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Mean power of rectangles, weighted by their trust values
 % (What if two speeds of pump?)
  power = [rectangles(:).power] .* skipNaN(@sum,[rectangles(:).trust], 1);
  power = sum(power) / sum(skipNaN(@sum,[rectangles(:).trust]));

  if nargout >= 2
    trust = zeros(1,length(rectangles));
    for r = 1:length(rectangles)
      trust(r) = -abs(log(rectangles(r).power/power));
    end
  end
end

function [t, wgts, p3, jumps] = get_weighted_edge (time, data, power, g, offset, turn_on)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % A transition part way through an interval causes a mid-level power
  % Is floor(time) the mid point or end of the transition?
  if abs(time - round(time)) < 0.2
    mid = round(time) + [-1, 0];
  else        % time middle of the interval: floor(time) is the mid-point
    mid = floor(time);
  end
  if any(mid <= 0)
    mid = mid + size(data,1)/2;
  end

  % f
  wgts = zeros (length (mid), size (data, 2));
  for m = 1:length(mid)
    %transition = data(mid(m));
    after  = (mid(m)+1:   mid(m)+g);    % periods before and after.  Element 1
    before = (mid(m)-1:-1:mid(m)-g);    % is closest to the transition in both
    if before(end) <= 0
        before = before + size(data,1)/2;
    end

    if turn_on
      is_off = data(before,:);
      is_on  = data(after,:);
    else
      is_off = data (after,:);
      is_on  = data (before,:);
    end

    wgts(m,:) = 1 ./(offset + (is_off(1,:)) / power) ...
              + 0.5./(1 + (max(is_off, [], 1) - min(is_off, [], 1)) / power) ...
              + 0.5./(1 + (max(is_on,  [], 1) - min(is_on,  [], 1)) / power);
  end

  w = sum(wgts, 2);
  [~, mx] = max (w);  % pick time that looks most like a transition
  %w    = w   (mx);
  mid  = mid (mx);
  wgts = wgts(mx,:);

  p3    = wgts * [data(mid + 1,:)                     % after
                  data(mid - 1 + size(data,1)/2, :)   % before
                  data(mid,:)]';                      % middle
  frac = (p3(3) - p3(1)) / (p3(2) - p3(1));
  t = mid + max(0, min(1, frac));
  if t < 0.5
    t = t + size(data,1)/2;
  elseif t > size(data,1)/2 + 0.5
    t = t - size(data,1)/2;
  end
  if mid > 1
    jumps = data(mid+1,:) - data(mid-1,:);
  else
    jumps = data(mid+1,:) - data(mid-1+size(data,1)/2,[2:end, end]);
  end
end

function [pump, power, cor_pool] = seems_like_pool (rectangles, cor_pool, c, meta, valid_days, stage)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (isempty(rectangles))
    pump = 0;
    power = 0;
    cor_pool = 0;
    return
  end
  final = 3;
  switch stage   % mc = mean correction (estimate of power)
    case 1
      on_frac = 0.4; mc_min = 0; mc_max = 100; big_min = 0.26; %allday_frac = 0.2;
    case 2
      on_frac = 0.6; mc_min = 0.1; mc_max = 2; big_min = 0.3; %allday_frac = 0.5;
    case final
      on_frac = 0.8; mc_min = 0.1; mc_max = 2; big_min = 0.39; %allday_frac = 0.5;
  end


  % Guess from rectangles and cor_pool whether this is a pool pump.
  pump = 0;
  %-% A. Most days have some timed device, within wider power range
  days_covered = sum(full(sum(cor_pool, 2) ~= 0));
  days_covered_summer = sum (full (sum (cor_pool(meta.summer,:), 2)) ~= 0);
  mean_cor = full(-mean(cor_pool(cor_pool(:) ~= 0)));
  summer = meta.summer;
  if (length(valid_days) < size(c,2))
    [~, summer] = ismember (summer, valid_days);
    summer = summer (summer ~= 0);
  end
  hours_per_day = sum (cor_pool(:) ~= 0) / days_covered / 2;
  if (days_covered > on_frac * sum(~isnan(c(1,:))) ...
       || days_covered_summer     > 0.8 * max (length (summer), 20)...
     ) && mean_cor > mc_min ...
       &&(mean_cor < mc_max/2 ...
          || (mean_cor < mc_max && hours_per_day > 4 && hours_per_day < 15))

    power = full(-mean(cor_pool(cor_pool(:) ~= 0)));
    pump = 1;
  else
    %-% B. "big" rectangles, not typical of noise, within limited power range
    on_off = [rectangles.on_off];
    st = on_off(1,:); turn_on = on_off(2,:);
    en = on_off(3,:); turn_off = on_off(4,:);
    durations = turn_off - turn_on;
    overnight = (durations < 0);
    durations(overnight) = durations(overnight) + meta.SamPerDay;
    big_rectangles = ((durations > 6 | (durations >= 3 & turn_on > 20)) ...
                      & (en - st > 40 ...
                         | (en - st > 20 & sum ([rectangles(:).trust],1) > 40)));
    if isequal (meta.hemisphere, 'south')
      summer_rect = big_rectangles & ((valid_days(st) > 250 | valid_days(en-1) < 150) ...
                                     |(valid_days(st) > 225 & en - st > 60) ...
                                     |(valid_days(en-1) < 175 & en - st > 60));
    else
      summer_rect = big_rectangles & (valid_days(st) < 330 & valid_days(en-1) > 70);
    end
    pool_rect = summer_rect & ([rectangles.power] > big_min ...
                             & [rectangles.power] < 0.8);
    for r = find (pool_rect)
%should fail:
%21.889071    5.133345    0.592593   -2.000000    0.607425   -0.031258 19.444903   31.844259
%should pass:
%54.77089    4.15663    0.79661   -2.00000    1.38793   -0.62901    5.86031 2.21109

                            % TODO: set these based on stage
      thresh = bsxfun(@gt, rectangles(r).trust, ...
                      [30;1.5;2.9;-1.1;10;-0.5; 20; 20]);
      if skipNaN(@sum,[thresh(1); thresh]) >= 4
        pump = 2;
        pool_rects = rectangles(pool_rect);
        power = mean ([pool_rects.power]);
      end
    end
  end

      %-% C. pump not on a timer, but turned on for days at a time
      %%%% TODO: Move this to has_untimed_pump
  if pump == 0
    % Take 1: just minimum daily power
                    % TODO: Why recalculate vampires?
    power = full(-mean(cor_pool(cor_pool(:) ~= 0)));
    on_all_day = (min(c) > max(0.9 * power, 0.2));
    may_match = [];

    if any (on_all_day)
      % Take 2: jumps in minimum and mean daily power
      [~, days, jp] = find_jumps(medfilt1 (min (c)), 0.05, 1);
      days = round (days);
      candidates = (abs (jp) > 0.7 * power);
      if any (candidates)
        day_mins = [min(c)'; 0];
        mean_diffs = day_mins (round (days(candidates)+1)) ...
                   - day_mins (round (days(candidates)));
        candidates(candidates) = candidates(candidates) ...
                                & abs  (mean_diffs) > 0.5 * power ...
                                & sign (mean_diffs) == sign (jp(candidates));
        % replaced intervals between jump candidates by a single mean
        candidates = find (candidates);
        if ~isempty (candidates)
          alternating = zeros (length (candidates), 1);
          alternating(1) = mean (day_mins(1:days(candidates(1))));
          for i = 1:length (candidates)-1
              alternating(i+1) = mean (day_mins(days(candidates(i))+1 ...
                                              :days(candidates(i+1))));
          end
          alternating(length (candidates)+1) ...
                 = mean (day_mins (days(candidates(end))+1:end));
          alternating(~isfinite (alternating)) = 0;
          diffs = diff(alternating);
          days = [1; days; size(cor_pool, 2)];

          may_match = zeros (length (candidates), 5);
          matches = 0;
          for i = 1:length(candidates)
            if sign (jp(candidates(i))) == sign (diffs(i))
              [time, height, my_jp, mate_jp] ...
                                  = match_jump(alternating, i, power);
              if my_jp > 0 && time > i
                st = i;
                en = round (time);
                if (days(en) > days(st) + 3)
                                  % ignore short runs due to aircon on hot nights
                  matches = matches + 1;
                  may_match(matches,:) = [st, en, height, my_jp, mate_jp];
                end
              end
            end
          end     % for i
          may_match = may_match (1:matches, :);
        end       % if ~isempty (candidates)
      end         % if any(candidates)

      if ~isempty (may_match)
        may_match = may_match(may_match(:,end)<0, :);
        day_power = min (abs (may_match(:,3:5)), [], 2);
        idx = (day_power > 0.3);
        may_match = may_match(idx,:);
        day_power = day_power(idx,:);

        if ~isempty (may_match)
          may_match(:,2) = round (may_match(:,2));
          for i = 1:size (may_match, 1)
            cor_pool(days(may_match(i,1)):days(may_match(i,2)),:) = -day_power(i);
          end
          pump = 3;
        end
      end       % if ~isempty (may_match)
    end         % if any on all day

%    vampires = find_vampires(c' + cor_pool);
%    cor_pool(on_all_day,:) = ...
%           repmat(-min(vampires(on_all_day), power), [1, size(cor_pool,2)]);
%
%    if sum(on_all_day(summer)) > 0.5 * max (length (summer), 20)...
%       && mean_cor > mc_min ...
%       &&(mean_cor < mc_max/2
%          || (mean_cor < mc_max && hours_per_day > 4 && hours_per_day < 15))
%        pump = 1;
%    end
  end           % if pump == 0

  if 0 && stage > 1
    fprintf('on_frac == %g mean_cor == %g\n', days_covered / sum(~isnan(c(1,:))), mean_cor);
    figure(10);
    imagesc(max(0,c + cor_pool'));
    p = 'p';
    fprintf('Power: %g\n', power);
    %mid = [mid, power]
    if stage == final
        ny = 'NY';
    else
        ny = 'NU';          % "unknown" if not yet passed all check
    end
    fprintf('My guess is %c\n', ny(pump+1));
    while p == 'p' || p == 'P'
      p = input ('Pool pump? (Y/N/U/P) ', 's');
      if isempty(p), p = ny(pump+1); end
      if length(p) > 1,  p=p(1); end
      if p == 'p'
        figure(1);
        c1 = c; c1(isnan(c1)) = 0;
        plot(c1(:,meta.summer));
        keyboard;
      end
    end
    switch upper(p)
      case {'Y','YES'}, pump = 1;
      case {'N','NO'},  pump = 0;
      case {'U'},       pump = -1;
      otherwise
          fprintf('Unknown option.  Assuming  Unknown\n');
          pump = -1;
    end
  end
  if pump == 0
    power = 0;
  end
end

function cor_pool = set_cor_pool (rectangles, cor_pool, valid_days)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cor_pool = zeros(size(cor_pool));
  for r = 1:length(rectangles)
    on_off = rectangles(r).on_off;
    st = on_off(1); turn_on = on_off(2); en = on_off(3); turn_off = on_off(4);
    ht = rectangles(r).power;
    disp ([r, on_off(:)', ht]);
    missed = rectangles(r).missed;
    if turn_on < turn_off
      burst = ceil(turn_on):floor(turn_off);        % turn_off must be last
    else                                        % as (1:end-1) is used
      burst = [ceil(turn_on):size(cor_pool,2), 1:floor(turn_off)];    % below
    end
    b = burst(1:end-1);
    cor_pool(valid_days(st:en-1),b) = cor_pool(valid_days(st:en-1),b)-ht;
    cor_pool(valid_days(missed), b) = cor_pool(valid_days(missed), b)+ht;
  end
end

function cor_pool = set_cor_pool_precise (rectangles, cor_pool, valid_days)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Calculate the power consumed by the pump at each time
 % Consider turn on/off within a slot (fractional on/off times) and
 % alternative day-of-week settings.
  cor_pool = zeros(size(cor_pool));
  for r = 1:length(rectangles)
    on_off = rectangles(r).on_off;
    st = on_off(1); turn_on = on_off(2); en = on_off(3); turn_off = on_off(4);
    ht = rectangles(r).power;
    missed = rectangles(r).missed;

    days = st:en-1;
    alt_days = [];
    if isfield (rectangles(r), 'alt_days') && sum (rectangles(r).alt_days)
      idx = ismember (mod1 (days, 7), rectangles(r).alt_days);
      alt_days = days(idx);
      days = days(~idx);
    end
    if turn_on < turn_off
      burst = ceil(turn_on):floor(turn_off);        % turn_off must be last
    else                                        % as (1:end-1) is used
      burst = [ceil(turn_on):size(cor_pool,2), 1:floor(turn_off)];    % below
    end
    b = burst(1:end-1);
    d = valid_days (setdiff (days, missed));
    if isempty (alt_days)
      cor_pool(d,b) = min (cor_pool(d,b), -ht);
      frac = ht * (turn_on - floor (turn_on) - 1);
      if frac ~= 0
        if turn_on < 1
          turn_on = turn_on + size(cor_pool, 2);
        end
        cor_pool(d, floor (turn_on)) = min (cor_pool(d, floor (turn_on)), frac);
      end
      frac = ht * (floor (turn_off) - turn_off);
      if frac ~= 0
        if turn_off < 1
          turn_off = turn_off + size(cor_pool, 2);
        end
        cor_pool(d, floor (turn_off)) = min (cor_pool(d, floor (turn_off)), frac);
      end
    else
      % TODO: fix for the case that alt_days is on fewer hours.
      cor_pool(d,b) = min (cor_pool(d,b), -ht);

      turn_on  = rectangles(r).alt_hrs(1);
      turn_off = rectangles(r).alt_hrs(2);
      if turn_on ~= turn_off
        if turn_on < turn_off
          burst = ceil(turn_on):floor(turn_off);
        elseif turn_off > turn_on
          burst = [ceil(turn_on):size(cor_pool,2), 1:floor(turn_off)];
        end
        b = burst(1:end-1);
        b(b==0) = size(cor_pool,2);
        d = valid_days (setdiff (alt_days, missed));
        cor_pool(d, b) = min (cor_pool(d,b), -ht);
      end
    end
  end
end

function [base, spikes] = stripSpikes(c, peakRatio)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if nargin < 2
      peakRatio = 5;
  end
  idx = find(c(2:end-1) - c(1:end-2) > peakRatio * abs(c(3:end) - c(1:end-2)));
  base = c;
  base(idx+1) = (c(idx) + c(idx+2))/2;
  if nargout > 1
    spikes = c -base;
  end
end

function [tm_real, jp, jp_real, reliability, trend] ...
                              = find_reliability(run, d, ~, samPerDay,band)
% (~ is 'c')
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tm = run(3); st = run(1); en = run(2)-1;
  tt = mod1 (tm-band, samPerDay);
  slice = d(tt:tt+2*band, st:en);
  slice = slice(:, ~isnan(slice(1,:)));       % Ignore entirely-NaN days

  %jp = mean(slice(band+1,:));
  %sd = sqrt(var(slice(band+1,:) ));

  % Process a band around this run, to see if the jump is "significant"
  jp_all = mean(slice,2);
  sd_all = sqrt(var(slice,[],2));
  sd = sd_all(band+1);
  jp = jp_all(band+1);
  trend  = sum(jp_all - jp) / (2*band);

  % If jump is in the middle of the measurement slot, it affects two rows
  jp_all (sd_all > abs(jp_all)) = 0;  % ignore this row if jump < noise
  neighbours = jp_all([band, band+2]);        % one above / one below
  [~, better] = min(abs(jp - neighbours));
  jp_better = jp_all(better);
  if   sign(jp_better) == sign(jp) && ...
      (sign(jp_better) ~= sign(trend) || abs(jp_better) > 2*abs(trend)) && ...
      (jp_better+jp) ~= 0
    tm_real = tm + 2*(better-1.5) * (jp_better / (jp_better+jp));
    jp_real = jp + jp_better;
  else
    tm_real = tm;
    jp_real = jp;
  end

  reliability = abs(jp - trend) * sqrt(en-st) / (sd + 1e-12) ;
  if en-st < 5 && reliability > 14    % less than "significance" threshold,
    reliability = 14;               % but still may match a significant jump
  end
end

function [runs, runCount] = runs_from_sparse(sp, dd, c, band)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  runCount = 0;
  runs = zeros(2*size(c,1),3);
  % Find "runs": clusters of consecutive days with a jump at the same time
  for i = 1:size (c,1)
    d = diff(sp(i,:) ~= 0);
    swap = (sign(sp(i,1:end-1)) .* sign(sp(i,2:end)) == -1);
    if ~any(d)
      if (sp(i,1) ~= 0 && ~isnan(sp(i,1)))
        starts = 1;
        ends = size(c,2)+1;
      else
        continue
      end
    else
      starts1 = find(d>0 | swap) + 1;     % first day of run
      ends1   = find(d<0 | swap) + 1;     % first day after run
      if isempty(starts1)
        starts1 = 1;
      elseif isempty(ends1)
        ends1 = size(c,2) + 1;
      end
      if ends1(1) <= starts1(1)
        starts = [1, starts1];  % starts1 and ends1 to avoid Matlab warning
      else
        starts = starts1;
      end
      if length(starts) > length(ends1)
        ends = [ends1, size(c,2) + 1];
      else
        ends = ends1;
      end
      % Delete intervals that are all NaN
      if any(isnan(sp(i,:)))
        allNaN = zeros(size(starts));
        for j = 1:length(starts)
          if all(isnan(sp(i,starts(j)+1:ends(j))))
            allNaN(j) = 1;
          end
        end
        starts = starts(~allNaN);
        ends   = ends  (~allNaN);
      end
    end

    L = length(starts);
    if L > 0
      runs(runCount+(1:L),:) = [starts', ends', repmat(i,[L, 1])];
      runCount = runCount + L;
    end
  end
  runs = runs(1:runCount,:);  % truncate excess pre-allocated space

  % If any run is long, calculate "reliability" of all runs
  if max(runs(:,2) - runs(:,1)) > 14
    for i = 1:runCount
      [tm_real, jp, jp_real, reliability, trend] ...
              = find_reliability(runs(i,:), dd, c, size (c,1), band);
      runs(i,4:8) = [tm_real, jp, jp_real, reliability, trend];
    end
  else
    runs(:,4:8) = 0;        % if all short, skip that step
  end

  if runCount > 0
    runs = runs(~isnan(runs(:,7)),:);
    runCount = size(runs,1);
  else
    runs = zeros(0,8);    % Octave complains if 0x8 is compared with say 0x4
  end
end

function runs = merge_oscillating_runs (runs, cc, dd, band)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Merge runs oscillating between neighbouring times
  runCount = size(runs,1);
  toDelete = zeros(1,runCount);
  for i = 1:runCount
    if runs(i,3) < 0    % if already flagged for deletion
      continue
    end

    % Find runs that occur immediately after run i.
    if runs(i,2) < size(cc,2) + 1
      after1 = runs(:,1) == runs(i,2);
      after2 = after1;
      if runs(i,3) > 2
        after1 = after1 & (abs(runs(:,3)-runs(i,3))==1);
      end
      if runs(i,3) < size(cc,1) + 1
        after2 = after2 & (abs(runs(:,3)-runs(i,3))==1);
      end
      after = find (after1 | after2);
    else
      after = [];
    end

    rm = runs(i,7);
    rt = 0;
    rr = [];
    for k = 1:length(after)
      a = after(k);
      last = runs(a,2);
      aa = find(runs(:,3) == runs(i,3) & runs(:,1) == runs(a,2));
      if ~isempty(aa)
        if length(aa) > 1
          % Multiple overlapping "runs"
          % Merge with most reliable previous one
          [~, r] = max(runs(aa,7));
          aa = aa(r);
        end
        last = runs(aa,2);
      end

      run = [runs(i,1), last, runs(i,3)];
      [x1,x2,x3, r, x4] = find_reliability(run, dd, cc, size (cc,1),band);
      if r > rm && r > runs(a,7) && (isempty(aa) || r > runs(aa,7))
        rr = [a,aa];
        rt = runs(i,3);
        rm = r;
        x = [x1, x2, x3, x4];
      end

      run = [runs(i,1), runs(a,2), runs(a,3)];
      [x1,x2,x3, r, x4] = find_reliability(run, dd, cc, size (cc,1), band);
      if r > rm && r > runs(a,7)
        rr = a;
        rt = runs(a,3);
        rm = r;
        x = [x1, x2, x3, x4];
      end
    end

    if rm > runs(i,7)
      if rr(1) < i
        toDelete(i) = 1;
        runs(i,3) = -1;
        toExtend = rr(1);
      else
        toDelete(rr(1)) = 1;
        runs(rr(1),3) = -1;
        toExtend = i;
      end
      runs(toExtend,1) = runs(i,1);
      runs(toExtend,2) = runs(rr(end),2);
      runs(toExtend,3) = rt;
      runs(toExtend,4:8) = [x(1), x(2), x(3), rm, x(4)];
      if length(rr) > 1
        toDelete(rr(2)) = 1;
        runs(rr(2),3) = -1;
      end
    end
  end
  runs = runs(~toDelete,:);
end

function runs = merge_adjacent_runs(runs, cc, dd, band)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [~, idx] = sort(runs(:,3) + runs(:,1)/size(cc,2));
  runs = runs(idx,:);
  if (size(runs,1) > 1)
    candidates = find(diff(runs(:,3)) == 0 & diff(sign(runs(:,5))) == 0 ...
                & (runs(1:end-1,7) > 10 | runs(2:end,7) > 10));
    for i = candidates(:)'
      st = runs(i,  2);
      en = runs(i+1,1);
      row= runs(i,  3);
      prev = mod1 (row-1, size(cc,1));
      if runs(i,5) < 0
        top = prev;
      else
        top = row;
      end
      if (all(cc(top, st:end) > max(runs(i:i+1,5))) && ...
           en - st < 0.2 * mean(runs(i:i+1,2) - runs(i:i+1,1)) && ...
           sign(mean(cc(row,st:end) - cc(prev,st:end))) == sign(runs(i,5)))
        run = [runs(i,1), runs(i+1,2), row];
        [x1,x2,x3,x4,x5] = find_reliability(run, dd, cc, size (cc,1), band);
        if x4 > runs(i,7)
              % merge into i+1, which is further checked next iter'n
          runs(i+1,1) = runs(i,1);
          runs(i+1,4:8) = [x1, x2, x3, x4, x5];
          runs(i,3) = -1;             % flag for deletion
        end
      end
    end
  end
  runs = runs(runs(:,3) ~= -1,:);   % delete runs flagged above for deletion
end

function [rectangles, runs, cancelled] = runs_to_rectangles (runs, cv, dd, band)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [~, m] = sort(-runs(:,7));
  runs = runs(m,:);
  cancelled = zeros(size(cv'));
  cc1 = max (cv, 0);
  runCount = size(runs,1);
  i = 0;
  partners = cell (runCount, 1);    % rectangles taken from this run
%runCount
%show_runs (runs, cv, 1:size(cv,2))

  rectCount = 0;
  rectangles = zeros (0,0);   % prevent Matlab warning
  if runCount == 0
    return;
  end
  rectangles(runCount).on_off = [];     % Force it to a structure
  rectangles(runCount).power  = [];     % Force it to a structure
  rectangles(runCount).edges  = [];     % Force it to a structure
  rectangles(runCount).missed = [];     % Force it to a structure
  rectangles(runCount).burst  = [];     % Force it to a structure
  rectangles(runCount).trust  = [];     % Force it to a structure

  while i < runCount
    first_partner = 0;    % Haven't split any runs yet
    second_partner = 0;
    i = i + 1;
%fprintf ('i = %d\n', i);
%show_runs (runs, cv, 1:size(cv,2))
    if runs(i,2) <= runs(i,1)  || runs(i,3) == -1
      continue
    end
    if runs(i,2) - runs(i,1) > 14     % long run -> average over day-of-week
      period = 7;
      same_each_day = 0;
    else                        % short run -> average over whole time
      period = 1;
      same_each_day = 1;
    end
    time    = zeros(1,period);
    height  = zeros(1,period);
    interval = cell(1,period);
    if period ~= 1
      for j = 1:period
        interval(j) = {runs(i,1)+j-1:period:runs(i,2)-1};
        x = mean(quantile(cc1(:,interval{j}), [0.05 0.1 0.2 0.4], 2), 2)';
                    % - mean(vampires(interval{j}));
        % Find matching jump...
        [time(j), height(j)] = match_jump (x, runs(i,3), runs(i,6));
      end
      % TODO: should account for wrapping around midnight
      if any(abs(time - mean(time)) > 1)      % times may vary by day of week
        mx = max(time);
        mn = min(time);
        if max(min(abs(time - mx), abs(time-mn))) < 1
          time = NaN;
          same_each_day = 1;
        else
          same_each_day = 1;
        end
      else
        same_each_day = 1;
      end
    end
    % If doesn't vary by day of week, recalculate using whole block.
    if same_each_day
      %period = 1;
      interval = runs(i,1):runs(i,2)-1;
      x = mean(quantile(cc1(:,interval), [0.05 0.1 0.2 0.4], 2), 2)';
      [time, height, my_jp, mate_jp] = match_jump (x, runs(i,3), runs(i,6));
      up   = max(my_jp, mate_jp);
      down = min(my_jp, mate_jp);
    end

    if height(1) < 0  % what does this correspond to?
      continue
    end
    mate = find_mate (runs, i, time);

    % If this is reliable, but no match was found, try harder
    if isempty (mate) && runs(i, 7) > 15
%fprintf('Trying %d again\n', i);
      % Merge with nearby runs at the same time
      nearby = (runs(:,3) == runs(i,3) & runs(:,1) < runs(i,2) + 7 ...
                & runs(:,2) > runs(i,1) - 7 & runs (:,7) > 12 ...
                & abs (runs(:,5) - runs(i,5)) < 0.2 ...
                & runs (:,3) ~= -1 & (1:size (runs,1))' > i);
      if any (nearby(i+1:end))    % merge all nearby runs
        runs(i,[1 2]) = [min(runs(nearby, 1)), max(runs(nearby, 2))];
        interval = runs(i,1):runs(i,2)-1;
        nearby(i) = false;
        runs(nearby, 3) = -1;  % flag for deletion
      end
      x = mean(quantile(cc1(:,interval), [0.05 0.1 0.2 0.4], 2), 2)';
      [time, height, my_jp, mate_jp] = match_jump (x, runs(i,3), runs(i,6));
      up   = max(my_jp, mate_jp);
      down = min(my_jp, mate_jp);
      mate = find_mate (runs, i, time);

      if (isempty (mate))
        % See if "time" specifies run that was not detected
        % Try with several quantiles
        % Try ignoring all seven sets of two days
        D = [1 2 3 4 5; 2 3 4 5 6; 3 4 5 6 7; 4 5 6 7 1; 5 6 7 1 2; 6 7 1 2 3];
        v = zeros (1, size (D,1));        % variancees
        e = zeros (1, size (D,1));    % reliability of edges
        t = zeros (1, size (D,1));    % times of edges
        tt = zeros (1, size (D,1));   % "real" times
        jp = zeros (1, size (D,1));   % jumps
        jp_r = zeros (1, size(D,1));   % jumps corresponding to "real" times
        tr = zeros (1, size (D,1));   % trends
        m = cell (1, size (D,1));     % mates
        for j = 1:size (D,1);
          days = D(j,:);
          ival = interval(ismember (mod1 (interval, 7), days));
          if isempty (ival)
            continue
          end
          x = mean(quantile(cc1(:,ival), [0.05 0.1 0.2 0.4], 2), 2)';
          [time, height, my_jp, mate_jp] = match_jump (x, runs(i,3), runs(i,6));
          up   = max(my_jp, mate_jp);
          down = min(my_jp, mate_jp);
          region = index1 (cc1, round(up), round(down), ival);
          v(j) = var (region(:));
          t(j) = time;
          m{j} = find_mate (runs, i, time);


          if ~isempty (m{j})
            [e(j), idx] = max (runs([m{j}], 7));
            tmp = m{j};
            m{j} = tmp(idx);
          elseif length (ival) < 2
            continue
          else
            c_ival = cv(:,ival);
            d_ival = dd(:,ival);
            run = [1, length(ival), round(t(j))];
            [tt(j), jp(j), jp_r(j), e(j), tr(j)] ...
                     = find_reliability(run, d_ival, c_ival, size (cv,1), band);
          end
          %fprintf ('%f ', mate_jp, time, e(j), tt(j), jp(j), jp_r(j), tr(j));
          if isempty (m{j})
            %fprintf('no match\n');
            m{j} = zeros(1,0);          % needed for [m{:}] below
          else
            %mj = m{j};
            %fprintf ('matched %d\n', mj(1));
            %runs(mj(1),:)
          end
        end
        mate = [m{:}];
        if length (mate) > 1
          % Choose "best" match.  One with highest e()?  For now, take first
          mate = mate(1);
        elseif isempty (mate)
          % More hacks
        end
      end
%     elseif isempty (mate)
%       quality = runs(i, 7)
    end

    if length(mate) == 1
      if runs(i,6) > 0
        turn_on = runs(i,4);
        turn_off = time(1);
      else
        turn_on = time(1);
        turn_off = runs(i,4);
      end
      if turn_on < 0.5
        turn_on = turn_on + size(cv,1);
      elseif turn_on >= size(cv,1)+0.5
        turn_on = turn_on - size(cv,1);
      end
      if turn_off < 0.5
        turn_off = turn_off + size(cv,1);
      elseif turn_off >= size(cv,1)+0.5
        turn_off = turn_off - size(cv,1);
      end
      ON  = ceil(turn_on);
      OFF = floor(turn_off);
      if ON < OFF
        burst = ON:OFF;
      else
        burst = [ON:size(cc1,1), 1:OFF];
      end

      % Align ends of runs
      tweak_second = [i, mate];
      st = runs(i,1);
      new_en = []; new_en_val = 0;        % for delaying truncation, below
      while runs(i,1) ~= runs(mate,1)
        % check interval in one run but before the other
        ival = min(runs(i,1),runs(mate,1)):max(runs(i,1),runs(mate,1))-1;
        if abs(runs(i,1) - runs(mate,1))<6   % close; tweak the ends
          x = mean(quantile(cc1(:,ival), [0.05 0.1 0.2 0.4], 2), 2)';
                      %- mean(vampires(ival));
          tm = match_jump(x, runs(i,3), runs(i,6));
          if (abs(tm - time) < 1)
            st = ival(1);
          else
            st = ival(end);
          end
          % Align starts, bring along any abutting runs
          runs(runs(:,2) == runs(i,1) & runs(:,3) == runs(i,3), 2) = st;
          runs(runs(:,2) == runs(mate,1) & runs(:,3) == runs(mate,3), 2) = st;
          runs([i, mate],1) = st;
          break;
        else                        % far; maybe a switch between
          st = ival(end)+1;
          if runs(i,1) > runs(mate,1)
            grow=i; shrnk=mate;
          else
            shrnk=i; grow=mate;
          end
              % Record which rectangle caused split.
              % At the end, we'll try to widen the rectangle
              % if  shrnk  is still "unused".
              % Positive means "extend to the left"
          first_partner = shrnk;

          before = find(runs(:,2) == runs(grow,1) & runs(:,3) == runs(grow,3));       % contiguous?
          if ~isempty(before) && sign(runs(before, 5)) == sign(runs(grow,5))
            runs(grow,1) = runs(before, 1);
            runs(grow,7) = max(runs(grow,7), runs(before,7));
            runs(before,3) = -1;
            st = runs(grow,1);      % in case while loop ends
            continue;
          end
                      % truncate this run, after the next phase
          new_en = shrnk; new_en_val = st;
          tweak_second = grow;       % don't re-tweak  fix
          break;
        end
      end

      en = runs(i,2);
      while runs(i,2) ~= runs(mate,2)
        ival = min(runs(i,2),runs(mate,2)):max(runs(i,2),runs(mate,2));
        if abs(runs(i,2) - runs(mate,2))<6   % close; tweak the ends
          x = mean(quantile(cc1(:,ival-1), [0.05 0.1 0.2 0.4], 2), 2)';
                      %- mean(vampires(ival));
          tm = match_jump(x, runs(i,3), runs(i,6));
          if (abs(tm - time) < 1)
            en = ival(end);
          else
            en = ival(1);
          end
          % Align starts, bring along any abutting runs
          for j = tweak_second(:)'
            runs(runs(:,1) == runs(j,2) & runs(:,3) == runs(j,3), 1) = st;
          end
          runs(tweak_second,2) = en;
          break;
        else                        % far; maybe a switch between
          en = ival(1);
          if runs(i,2) < runs(mate,2)
            grow=i; shrnk=mate;
          else
            shrnk=i; grow=mate;
          end
              % Record which rectangle caused split.
              % At the end, we'll try to widen the rectangle
              % if  shrnk  is still "unused".
              % Negative means "extend to the right"
          second_partner = shrnk;

          after = find(runs(:,1) == runs(grow,2) ...
                     & runs(:,3) == runs(grow,3));   % contiguous?
          if length (after) > 1
            after = after(end);
          end
          if ~isempty(after) && sign(runs(after, 5)) == sign(runs(grow,5))
            runs(grow,2) = runs(after, 2);
            runs(grow,7) = max(runs(grow,7), runs(after,7));
            runs(after,3) = -1;
            continue;
          end
                      % don't make big adjustments to both ends.
                      % If that is called for, just leave old run
          if any(shrnk == tweak_second)
            runs(shrnk,1) = en;
          else
            % Split run, rather than leaving it...
            new_run = runs(shrnk,:);
            new_run(2) = new_en_val;
            runs(shrnk,1) = en;

            runCount = runCount + 1;
            if mate > i
              mate = mate + 1;
            end
            runs = [runs(1:i,:); new_run; runs(i+1:end,:)];
            partners = partners([1:i, i:end]);    % keep partners(i)
            if first_partner > i                  % ... and new...
              first_partner = first_partner + 1;
            end
            if second_partner > i                 % ... partners
              second_partner = second_partner + 1;
            end

            % ...and don't shrink the original below
            new_en = []; new_en_val = 0;
          end
          break;
        end
      end
      runs(new_en,2) = new_en_val;        % do delayed truncation now

 if runs(i,2) < runs(i,1) || runs(mate,2) < runs(mate,2)
 fprintf('Error!!!!!!!!!!!\n');
 keyboard;
 end

      overlap = (runs(:,1) < en & runs(:,2) >= st & ismember(runs(:,3),burst) & (1:runCount)' ~= i & (1:runCount)' ~= mate);
      if any(overlap)
        runs(overlap & runs(:,1) <  st, 2) = st;
        runs(overlap & runs(:,2) >= en, 1) = en;
        runs(overlap & runs(:,2) <= en & runs(:,1) >= st-1, 3) = -1;
        if any((runs(overlap,1) > runs(overlap,2)) & runs(overlap,3) ~= -1)
          fprintf('***************Error\n');
          keyboard;
        end
      end
      %ht = min(height(1), max(runs(i,6), runs(mate,6)));
      ht = min(height(1), max(up, -down));
      slice = cc1(burst(1:end-1),st:en-1);
      missed = (min(slice,[],1) < 0.8 * ht);
      if ~all(missed)
        ht = min(ht, min(min(slice(:, ~missed))));
      end
      old_st_en = [st, en];
      if ~isempty (missed) && missed(1)
        offset = find (~missed, 1);
        if ~isempty (offset)
          st = st + offset - 1;
        end
        missed = missed(offset:end);
      end
      if ~isempty (missed) && missed(end)
        offset = find (~missed(end:-1:1), 1) - 1;
        en = en - offset;
        missed = missed(1:end-offset);
      end
          % if mostly "missed", probably not really a block
          % More reliable runs can miss a higher fraction to be trusted
      if  sum(missed) > (en - st) * (1-exp(-min(runs([i,mate],7))/20))
        continue
      end

      clear new_rectangle;
      new_rectangle.on_off = [st; turn_on; en; turn_off];
      new_rectangle.power  = ht;
      new_rectangle.edges  = [runs(mate,6); height(1); runs(i,6)];
      new_rectangle.missed = st + find(missed) - 1;
      new_rectangle.burst   = burst;
      new_rectangle.trust  = (runs(i,7) * runs(mate,7)) * (en-st) * (mod(turn_off - turn_on, size(cv,1))+1);

      rectCount = rectCount + 1;
      rectangles(rectCount) = new_rectangle;
      if first_partner > 0
        partners{first_partner} = [partners{first_partner}, rectCount];
      end
      if second_partner > 0
        partners{second_partner} = [partners{second_partner}, -rectCount];
      end

      b = burst(1:end-1);
      missed = new_rectangle.missed;
      cancelled(st:en-1,b) = cancelled(st:en-1,b)-ht;
      cancelled(missed, b) = cancelled(missed, b)+ht;
      cc1(b,st:en-1) = cc1(b,st:en-1) - ht;
      cc1(b,missed)  = cc1(b,missed)  + ht;

      % delete this run if it is entirely matched, or keep matching it
      ht_ratio = abs(runs(i,6)) / ht;
      if all (runs(i,1:2) == old_st_en)
        if ht_ratio > 0.5 && ht_ratio < 2
          runs(i,3) = -1;
        end
      else
        i = i - 1;
      end

      % delete this run if it is entirely matched
      ht_ratio = abs(runs(mate,6)) / ht;
      if runs(mate,1) == st && runs(mate,2) == en && ht_ratio > 0.5 ...
         && ht_ratio < 2
        runs(mate,3) = -1;
      end
    else
      % What do we do if length(mate)~=1?
    end         % if length(mate) == 1
  end           % while i < runCount

  rectangles = rectangles (1:rectCount);

  % Delete runs that have been set to empty
  partners = partners(runs(:,3) >= 0);
  runs = runs(runs(:,3) >= 0,:);

  % Try to extend rectangles where one run was much longer than the other
  trusted = find (runs(:,7) > 10);
  changed = 0;
  for i = trusted'
    %fprintf ('%d before: ', i)
    for p = partners{i}
      %rectangles(abs(p)).on_off'
      if p > 0  % extend to the left
        on_off = rectangles(p).on_off;
        leftmost = max(1, runs(i,1)-2);
        range = leftmost:on_off(1)-1;
        if length (range) > 1
          m = index1 (cv, round (on_off(2)), round (on_off(4)), range);
          m = mean (m, 1);
          mean_LR = cumsum (m) ./ (1:length (m));
          mean_RL = cumsum (m(end:-1:1)) ./ (1:length(m));
          mean_RL = mean_RL(end:-1:1);
          difference = mean_RL(2:end) - mean_LR(1:end-1);
          [~, switchover] = min (abs (difference - rectangles(p).power));
          if rectangles(p).on_off(1) ~= switchover + leftmost
            changed = 1;
          end
          rectangles(p).on_off(1) = switchover + leftmost;
        end
      else      % extend to the right
        on_off = rectangles(-p).on_off;
        rightmost = min (runs(i,2)-1, size (cv,2));
        range = on_off(3)-1:rightmost;
        if length(range) > 1
          m = index1 (cv, round (on_off(2)), round (on_off(4)), range);
          m = mean (m, 1);
          mean_LR = cumsum (m) ./ (1:length (m));
          mean_RL = cumsum (m(end:-1:1)) ./ (1:length(m));
          mean_RL = mean_RL(end:-1:1);
          difference = mean_LR(2:end) - mean_RL(1:end-1);
          [~, switchover] = min (abs (difference - rectangles(-p).power));
          if rectangles(-p).on_off(3) ~= switchover-1 + on_off(3)
            changed = 1;
          end
          rectangles(-p).on_off(3) = switchover-1 + on_off(3);
        end
      end
      %rectangles(abs(p)).on_off'
    end
    %fprintf ('after: ')
  end

  % Recalculate cancellation if required
  if changed
    %[~, r, ok] = pool_ground_truth (cv, [], rectangles);
    ok = false | false;
    if ok
      rectangles = r;
      cancelled = zeros (size (cv'));
      for i = 1:length (rectangles)
        b = rectangles(i).burst(1:end-1);
        missed = rectangles(i).missed;
        st = rectangles(i).on_off(1);
        en = rectangles(i).on_off(3);
        ht = rectangles(i).power;
        cancelled(st:en-1,b) = cancelled(st:en-1,b)-ht;
        cancelled(missed, b) = cancelled(missed, b)+ht;
      end
    end
  end
end

function mate = find_mate (runs, me, time)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Given an estimate  time  of the off/on paired with run  me, find a run
 % that seems to be that pair.
  mate = find(runs(:,1) < runs(me,2) & runs(:,2) > runs(me,1) ...
              & abs(runs(:,3) - time(1)) <= 1 & (1:size (runs,1))' ~= me ...
              & sign(runs(:,5)) == -sign(runs(me,5)) & runs(:,1) ~= runs(:,2));
  if length(mate) > 1
    % find possible matching transition with greatest overlap with
    % the match.
    [o_lap_start, os] = max(runs(mate & runs(mate,2) < runs(me,2), 2) - runs(me,1));
    [o_lap_end, oe] = max(runs(me,2) - runs(mate & runs(mate,1) > runs(me,1), 1));
    if o_lap_start > o_lap_end
      mate = mate(os);
    else
      mate = mate(oe);
    end
  end
end

function corr = solar_correction (data, on, off)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Return either an empty matrix, or a quadratic compensation for
 % under-estimated solar generation.
  corr = [];
  on = round (on);
  off= round (off);

  % Max/min may not be quite at edges of rectangle, so search for them.
  y = median (data(on:off, :), 2);
  [~, pos] = max (y(1:round(end * 0.25)));
  on = on + pos - 1;
  [~, pos] = max (y (round (end * 0.75):end));
  off = off - pos + 1;
  y = median (data(on:off, :), 2)';    % TODO: avoid recalculation of median

  % duration
  d = off - on;
  d2 = d * d;

  % ordiniates
  x = 0:d;
  x2 = x .* x;
  x3 = x .* x2;
  x4 = x2 .* x2;

  % Fit a quadratic to data(on:off);
  sx3 = sum (x3);
  A = inv ([2*sum(x4), 2*sx3; 2*sx3, 2*sum(x2)]);
  Xd = A * [d2; d];
  XE = A * [2*sum(x2 .* (y(1) - y));
            2*sum(x  .* (y(1) - y))];

  % Lagrange multiplier
  lambda = (y(end) - y(1) + [d2 d] * XE) / ([d2 d] * Xd);

  ab = Xd * lambda - XE;

  % if power dips during the rectangle
  if ab(1) > 0 && ab(2) < 0 && ab(2) >= -2 * d * ab(1)
    % check if dip is substantial
    fit = ab(1) * x2 + ab(2) * x;
    fit_err = mean (abs (y - y(1) - fit));
    if min (fit) + fit_err < 0
      % Dips is substantial, so correct between 10am and 6pm.
      first = min (on,  round (size(data, 1) / 2.4));
      last  = max (off, round (size(data, 1) * 0.75));
      x = (first:last) - on;
      corr = (ab(1) * x + ab(2)) .* x + y(1);
      corr = max (0, min (corr(1), corr(end)) - corr);
      corr = [zeros(first-1, 1); corr'; zeros(size(data,1)-last, 1)];

      % If there are jumps explaining the dip,
      % disregard the possibility of solar.
      [mx, pos] = max (corr(on:off-1));
      d = diff (y);
      if -2 * min (d(1:pos)) > mx && 2 * max (d(pos:end)) > mx
        corr = [];
      end
    end
  end
end

function rectangles = align_rectangles (rectangles, cv, d)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Match up rectangles, trying to identify days at which the timer changed
 % TODO: Estimate *time* of day when timer changed.  Some rectangles will
 %       be before the change and some after.  Postprocess, or now?
  %show_deletions = false;
  days = size(d,2);
  gap = 5;
  rectangles = end_days (rectangles([rectangles(:).power] > 0), cv);
  if isempty (rectangles)
      return
  end

  on_off = round([rectangles.on_off]);
  evDays = union(on_off(1,:), on_off(3,:));
  %toTweak = find(diff(evDays) < gap);
  [g, toTweak] = sort(diff(evDays));     % process smallest gaps first
  toTweak = toTweak (g < gap);
  for i = toTweak(:)'
    e  = evDays(i);
    g = gap;
    while (true)
      ending  = ((on_off(3,:) >=e ) & (on_off(3,:) <= e+g));
      starting= ((on_off(1,:) >=e ) & (on_off(1,:) <= e+g));
      if ~any(ending & starting)
        break;
      end
      g = min (on_off(3,starting)) - e - 1;
    end

    %both = (ending & starting);         % ignore settings of a few days
    ending   = find(ending);
    starting = find(starting);
    % Ignore rectangles for which the end position is clear
    e1 = max([on_off(3,ending), on_off(1,starting)]);
    if isempty(e1) || e1 == e      % this gap closed in a previous step?
      continue
    end

    locks = zeros (1, length (ending) + length (starting));
    lock_count = 0;
    for r = length(ending):-1:1
      q = rectangles(ending(r)).quality_R;
      [mx, pos] = max(q);
      q(pos) = -1e300;
      mx1 = max(q);
      if mx > 2.5*mx1
            % flag for "position at most 1 away from  pos"
        rectangles(ending(r)).on_off(3) = -pos-1;
        lock_count = lock_count + 1;
        locks(lock_count) = pos;
      end
    end
    for r = length(starting):-1:1
      q = rectangles(starting(r)).quality_L;
      [mx, pos] = max(q);
      q(pos) = -1e300;
      mx1 = max(q);
      if mx > 2.5*mx1
            % flag for "position at most 1 away from  pos"
        rectangles(starting(r)).on_off(1) = -pos-1;
        lock_count = lock_count + 1;
        locks(lock_count) = pos;
      end
    end
    locks = locks (1:lock_count);

    mask_before = zeros(size(cv,1),1);
    mask_after  = zeros(size(cv,1),1);
    tr = [rectangles([ending,starting]).trust];
    pwr = sum([rectangles([ending, starting]).power] .* tr) / sum ([rectangles([ending,starting]).trust]);
    for j = (ending(:)')
      if on_off(4,j) > on_off(2,j)
        mask_before(on_off(2,j)+1:on_off(4,j)-1) = rectangles(j).power;
      else
        mask_before(1:on_off(4,j))   = rectangles(j).power;
        mask_before(on_off(2,j):end) = rectangles(j).power;
      end
    end
    for j = (starting(:)')
      if on_off(4,j) > on_off(2,j)
        mask_after(on_off(2,j)+1:on_off(4,j)-1) = rectangles(j).power;
      else
        mask_after(1:on_off(4,j))   = rectangles(j).power;
        mask_after(on_off(2,j):end) = rectangles(j).power;
      end
    end
    mask_diff = mask_after - mask_before;
    idx_diff  = (mask_diff ~= 0);
        % idx_diff can be empty if =height rectangles separated by gap
    wgt_diff  = max(ceil(sum(idx_diff)/2), 1);
    int = max(1,min(days, (e:e1)-1));       % e:e1, clipped %%%%%%%%%%%
                    % TODO:  was cc.  Should it be??
    mis_before_penalty = [0,cumsum(any(bsxfun(@lt, cv(:,int), 0.9*mask_before)))];
    mis_after_penalty  = [0,cumsum(any(bsxfun(@lt, cv(:,int(end:-1:1)), 0.9*mask_after)))];
    mis_after_penalty = mis_after_penalty(end:-1:1);

    costs = zeros (1+e1-e,1);
    for j = 1:1+e1-e
      pre = max(e+j-2, 1);                                %%%%%%%%%%%
      post = min(e+j-1, days);                            %%%%%%%%%%%
                  % horizontal edges
      cost = sum(abs(d((on_off(2,ending)),pre)-pwr) ...
                +abs(d((on_off(4,ending)),pre)+pwr)) ...
           + sum(abs(d((on_off(2,starting)),post)-pwr) ...
                +abs(d((on_off(4,starting)),post)+pwr));
                  % vertical edges
                  % TODO:  was cc.  Should it be??
      cost = cost + (sum(abs(mask_diff(idx_diff) - (cv(idx_diff,post)-cv(idx_diff,pre))))/wgt_diff);
      cost = cost + mis_before_penalty(j)+mis_after_penalty(j+1);
      costs(j) = cost;
    end

    [best_cost, best_pos] = sort(costs);
    range = max(locks) - min(locks);
    clip = @(x) (max(e, min(e1, x)));
    if ~isempty(range)
      if any(range == [0, 1])
        if length(locks) == length(starting)+length(ending)
          % all start/end on the same day.
          % Find the time of the switch
          for j = ending(:)'
            rectangles(j).on_off(3) = clip(-rectangles(j).on_off(3));
          end
          for j = starting(:)'
            rectangles(j).on_off(1) = clip(-rectangles(j).on_off(1));
          end
        else  % not all rectangles locked
        %%%  TODO:  should do the following:
          % Those that don't have locked ends, check quality of
          % locked end dates. If at most twice the min, use that.
          % Otherwise, use max
          if best_cost(1) > 0.5 * best_cost(2)
            pos = e-1+best_pos(1);
          else
            pos = clip(floor(median(locks)) + 1);
          end
          for j = ending(:)'
            if rectangles(j).on_off(3) < 0
              rectangles(j).on_off(3) = clip(-rectangles(j).on_off(3));
            else
              rectangles(j).on_off(3) = pos;
            end
          end
          for j = starting(:)'
            if rectangles(j).on_off(1) < 0
              rectangles(j).on_off(1) = clip(-rectangles(j).on_off(1));
            else
              rectangles(j).on_off(1) = pos;
            end
          end
        end
      else      % range not 0 or 1
        %%%  TODO:  should do the following:
        % Fix the fixed ones.  For others, see which fixed one they
        % are closest to.
        % Find the times of the switches
        if best_cost(1) > 0.5 * best_cost(2)
          pos = e-1+best_pos(1);
        else
          pos = clip(floor(median(locks)) + 1);
        end
        for j = ending(:)'
          if rectangles(j).on_off(3) < 0
            rectangles(j).on_off(3) = clip(-rectangles(j).on_off(3));
          else
            rectangles(j).on_off(3) = pos;
          end
        end
        for j = starting(:)'
          if rectangles(j).on_off(1) < 0
            rectangles(j).on_off(1) = clip(-rectangles(j).on_off(1));
          else
            rectangles(j).on_off(1) = pos;
          end
        end
      end       % end if range == [0, 1]
    elseif ~isempty (best_cost)
      best_pos = best_pos(1);
      evDays(i)   = e-1+best_pos;
      evDays(i+1) = e-1+best_pos;
      on_off(3,ending)   = evDays(i);
      on_off(1,starting) = evDays(i);
      % Should be a way to write field of all elements of struct array
      for j = ([ending(:); starting(:)]')
          rectangles(j).on_off([1,3]) = on_off([1,3],j);
      end
    end
  end
  on_off = [rectangles(:).on_off];
  rectangles = rectangles(on_off(3,:) > on_off(1,:));
end

function rectangles = assign_trust (rectangles, cv, power, issolar, spike_ok)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Estimate the trustworthiness of the rectangles.
 % * Assign weights to days, high if the "off" side is low power
 %   ** weight =   1   / [1 + (off - vampires) / power)]
 %               + 0.5 / [1 + (max(n'hood off)-min(n'hood off) / power)]
 %               + 0.5 / [1 + (max(n'hood on) -min(n'hood on)  / power)]
 % * Look for threshold in weights and only consider good ones?
 %   ** Find weighted median to recalculate power
 %   ** Q_e = (1/sqrt(length)) \sum weight / (1 + |(|jump|)/power - 1 |)
 % * Identify fractional ON and OFF times
 % * Assess reliability of end times, especially if changing number of runs
 % * Higher reliability if noise over whole rectangle is small

 % Side effect: reassigns on time and off time, power, missed.

 % First three measures of trust based on the rectangle in isolation
 % (quality of jumps and a measure of the local noise)
 %  The rest are based on the compatibility of a rectangle with other
 % "trusted" rectangles.  Thus, each rating affects subsequent ones.
  %show_deletions = false;
  if isempty (rectangles)
    return;
  end
  if length (rectangles(1).trust) > 1
    for r = 1:length (rectangles)
      rectangles(r).trust = rectangles(r).trust(1);
    end
  end

  [~, idx] = sort(-([rectangles(:).trust]));
  rectangles = rectangles(idx);

  g = 4;
  %guard = zeros(2*size(cv,1), g);
  cvw = [cv; cv(:,2:end), cv(:,end)];        % wrap

  for r = 1:length(rectangles)
    on_off = rectangles(r).on_off;
    st = round(on_off(1)); turn_on = on_off(2); en = on_off(3); turn_off = on_off(4);

    run = st:en-1;
    if isempty(run)
      run = st;
    end

    % Weight the reliability of turn-off jump for each day,
    % and estimate time within the half hour that it turned off
    gg = round(turn_off) - round(turn_on);
    if gg < 1
      gg = gg + size(cv,1);
    end
    gg = min ([g, gg]);

    if (~spike_ok || turn_off > turn_on + 1)
      [turn_on, wgts_on, p3_on, jumps_on]  = get_weighted_edge (turn_on, cvw(:,run), power, gg, 1, true);
      [turn_off,wgts_off,p3_off,jumps_off] = get_weighted_edge (turn_off,cvw(:,run), power, gg, 1, false);
      if ((turn_off <= turn_on + 1 && turn_off >= turn_on))
        rectangles(r).trust = [0;0;0];
        continue
      end
    else
      mid    = floor(turn_on);
      after  = mid + 1;
      before = mid - 1;
      if before < 0
        before = before + size(cvw,1)/2;
      end

      d_before = cvw(before,run);
      is_on    = cvw(mid,   run);
      d_after  = cvw(after, run);

      jumps_on = is_on - d_before;
      jumps_off = d_after - is_on;

      power = mean (is_on - (d_before + d_after)/2);
      rectangles(r).power = power;

      wgts_on  = 1 ./ (1 + d_before / power);
      wgts_off = 1 ./ (1 + d_after  / power);

      p3_on  = wgts_on  * [cvw(mid, run)    % after
                           cvw(before, run) % before
                           cvw(mid, run)]'; % middle
      p3_off = wgts_off * [cvw(after, run)  % after
                           cvw(mid, run)    % before
                           cvw(mid, run)]'; % middle
    end
    if turn_on < turn_off
      rectangles(r).burst = ceil(turn_on):floor(turn_off);
    else
      rectangles(r).burst = [ceil(turn_on):size(cv,1), 1:floor(turn_off)];
    end

    % power: weighted sum of on transitions and off transitions
    on_jump  = p3_on (1) - p3_on (2);
    off_jump = p3_off(2) - p3_off(1);
    power = (on_jump + off_jump) / sum([wgts_on, wgts_off]);
    if power < 0
        rectangles(r).trust = [0;0;0];
        continue
    end

    % Find "missed" days; assign zero weight to those jumps.
    if turn_off < turn_on
      wrapped_off = turn_off + size(cv,1);
    else
      wrapped_off = turn_off;
    end
    body = cvw(floor(turn_on+1):floor(wrapped_off-1),run);
    if ~isempty (body)
      if ~issolar
        mins = min(body,[],1);
        %missed = (mins < power);
        missed= (mins < power * 0.8);   % may be some error in vampire estimate
                                        % so allow a small margin
      else
        mins = quantile(body,0.2,1);
        missed = (mins < power * 0.6);
        missed = missed & (jumps_on < 0.8*power) & (jumps_off < 0.8*power);
      end
    else
      missed = [];
    end
    wgts_on (missed) = -wgts_on (missed);     % if jumps also there on missed
    wgts_off(missed) = -wgts_off(missed);     % days, more likely noise
    Q_m = 0.3/(0.1+sum(missed,2)/(en-st));
    if isempty(Q_m)
      Q_m = 0;
    end

    % edge quality = weighted sum of jump qualities
    %jump = abs([p3_on(1)-p3_on(2), p3_off(1)-p3_off(2)]);
    %q = 1+abs(jump/power-1);
    %Q_e = sum([wgts_on/q(1), wgts_off/q(2)]);
    Q_e = sum([wgts_on, wgts_off]);

    %-% Calculate quality of body, by its uniformity
    % ? 2*fraction of points within 10% of min + 1*fraction within 20
    slice = cvw(floor(turn_on+1):floor(wrapped_off-1), run(~missed));
    %extra = bsxfun(@minus, slice, min(cvw(floor(turn_on)-1+size(cv,1),run(~missed)), cvw(floor(turn_off)+1,run(~missed))));
    %Q_b = (2 * sum(slice(:) < 1.1*power) + sum(slice(:) < 1.2*power)) / sqrt((en-st)*(wrapped_off-turn_on));
    %Q_b = sum(0.1 ./ (0.1 + abs(extra(:) / power - 1))) / length(extra(:));
    %Q_b = geo_mean(0.1+abs(extra(:)/power - 1)) * sqrt(en-st);
    var1d = @(x)(var ([x(:); mean(x(:))]));
    Q_b = 0.1 / (0.1 + var1d(slice(:))) * sqrt(en-st);


    rectangles(r).on_off = [st; turn_on; en; turn_off];
    rectangles(r).power = power;
    rectangles(r).missed = st-1 + find(missed);
    rectangles(r).trust = [Q_e; Q_b; Q_m];
  end

  rectangles = rectangles(any([rectangles(:).trust],1));

  if isempty (rectangles)
      return
  end
  %%%%
  % More trustworthy if we start when somebody else ends and vice versa
  on_off = round([rectangles.on_off]);
  on_off(2,on_off(2,:) < 1) = size(cv,1);
  on_off(4,on_off(4,:) < 1) = size(cv,1);
  on_off(2,on_off(2,:) > size(cv,1)) = 1;
  on_off(4,on_off(4,:) > size(cv,1)) = 1;
  neighbour_matches = zeros(1,size(cv,2)+1);
  for r = 1:size(rectangles)
      neighbour_matches([on_off(1,r), on_off(3,r)]) = ...
                      neighbour_matches([on_off(1,r), on_off(3,r)]) + 1;
  end
  use_times = zeros(size(cv,1), 1);
  for r = 1:length(rectangles)
    rectangles(r).trust(4) = -(neighbour_matches(on_off(1,r)) == 1) ...
                             -(neighbour_matches(on_off(3,r)) == 1);
    % recaculate burst, as on/off times changed
    ON  =  ceil (on_off(2,r));
    OFF = floor (on_off(4,r));
    if ON < OFF
      burst = ON:OFF;
    else
      burst = [ON:size(cv,1), 1:OFF];
    end
    if ~isempty(burst)
      use_times(burst) = use_times(burst) + skipNaN(@sum,rectangles(r).trust) * (en-st);
      rectangles(r).burst = burst;
    else
      rectangles(r).trust = 0*rectangles(r).trust;
    end
  end
  rectangles = rectangles(any([rectangles(:).trust],1));
  if isempty (rectangles)
    return
  end

  %%%%
  % More trustworthy if we occur at a "typical" time.
  use_times = filter ([0.25 0.25 0.25 0.25], 1, [use_times; use_times; use_times]);
  use_times = filter ([0.25 0.25 0.25 0.25], 1, use_times(end:-1:1));
  use_times = use_times(end:-1:1);
  for r = 1:length(rectangles)
    rectangles(r).trust(5) = mean(use_times(rectangles(r).burst(ceil(end/2)))) / size(cv,2);
  end

  [~, Q_p] = common_power(rectangles);

  on_off = round([rectangles.on_off]);
  for r = 1:length(rectangles)
    rectangles(r).trust(6:8) = [Q_p(r)
                rectangles(r).quality_L(max(1,min(size(cv,2), on_off(1,r)-1)))
                rectangles(r).quality_R(max(1,min(size(cv,2), on_off(3,r)-1)))];
  end
end

function rectangles = ds_trust (rectangles, cv, power, issolar, spike_ok)
  % Use Dempster-Shafer to calculate trust
  for r = 1:length (rectangles)
    % m_edge(Y) = exp(-(0.01 + sqrt (badness)))
    %             m_edge(Y+N) = 1 - m_edge(Y)
    %             m_edge(N) = 0
    % m_len(Y) = sqrt (len / size(cv,2))    if len > 6
    %            0 if len <= 6              if len <= 6
    %            m_len(Y+N) = 1 - m_len(Y)  if len > 6
    %                         len / 6 * (1-sqrt (6 / size(cv, 2))) if len <= 6
    %            m_len(N) = 0               if len <= 6
    %                       1 - m_len(Y+N)  if len <= 6
    % m_miss (Y) = 0
    %              m_miss(N) = misses / len
    %              m_miss(Y+N) = 1 - m_miss(N)
    %
    % m_power(Y) = 0;
    %              m_power(Y+N) = min(1, power / 2nd smallest power) *
    %                             min(1, (power / 0.25)^2)
    %              m_power(N) = 1 - m_power(Y+N)
  end
  % Assign initial trust
  %   Quality of horizontal edge
  %   Length
  %   Missed days
  %   Power less than half(?) the second-lowest power
  % Iteratively
  %   Calculate trust based on latest trust of other rectangles
  %     Size of overlap with trusted rectangles
  %     Start and/or stops on same day (or near) a trusted rectangle
  %     Location near the nearest trusted rectangle on either side
  %     Duration similar to the nearest trusted rectangle on either side
  %     Duration similar to all trusted rectangles
  %     Location near all trusted rectangles
  %     Power similar to trusted rectangles
  %     If two "chains" of rectangles:
  %       Calculate belief of being in each chain (1, 2, 1+2, ~1, ~2, ~(1+2))
  %       Change from previous reflects change in other chain
  %         (up by similar amount, up by amount the other is down
  %       Similar size to other chain
  %       Same on/off times as other chain
  %   Combine with initial trust
end

function rectangles = remove_noise_rectangles (rectangles, days)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Try to remove noise rectangles.
 % Usually on a fixed number of times per day (one or two)
 % Remove "low trust" rectangles that are in excess of the "typical" number
 % (TODO: Problem if, e.g., on twice in summer and once in winter).
  on_off = [rectangles.on_off];
  evDays = union(on_off(1,:), on_off(3,:));

  count = zeros (1, length (evDays) - 1);
  for i = 1:length(evDays)-1
    overlap = (on_off(1,:) < evDays(i+1) & on_off(3,:) > evDays(i));
    count(i) = sum(overlap);
  end
  usual_times = max(1, round (sum (count .* diff(evDays)) / days));


  [~, first] = ismember(on_off(1,:), evDays);
  [~, last]  = ismember(on_off(3,:), evDays);

  for r = length(rectangles):-1:3
    overlap = max(count(first(r):last(r)-1));
              % if 'too many' bursts today, and they have low 'trust', delete
    %thresh = bsxfun(@gt, rectangles(r).trust, [20;1.5;2.9;-1.1;10;-0.5; 20; 20]);
    thresh = [rectangles(r).trust([1,1]); rectangles(r).trust] > [220; 70; 20;1.5;2.9;-1.1;10;-0.5; 20; 20];
    if overlap > usual_times && 4 >= skipNaN(@sum,thresh)
      if on_off (3,r) - on_off (1,r) > 40
        % find sum of days exceeding usual_times
        de = diff (evDays);
        exceeds = (count(first(r):last(r)-1) > usual_times) * de(first(r):last(r)-1)';
        if exceeds / (on_off(3,r) - on_off(1,r)) < 0.5
          continue;
        end
      end
      rectangles(r).power = -1;
      count(first(r):last(r)-1) = count(first(r):last(r)-1) - 1;
    end
  end
  pwr = [rectangles(:).power];
  rectangles = rectangles(pwr > 0);
end

function qualities = scalar_quality (rectangles, chains, dists, chain_value, cv, issolar)
  % TODO: Should use Demster-Shafer
  ks_p = [rectangles.ks_p];
  idx = isnan (ks_p (1, :));
  ks_p(:, idx) = repmat (mean (ks_p, 2, 'omitnan'), 1, sum (idx));
  ks_p(isnan (ks_p)) = 0.9;   % if all NaN, like for a single long run

  on_off = [rectangles.on_off];
  %nr = zeros (2, length (rectangles));
  nr_on = zeros (3, length (rectangles));
  nr_off = nr_on;
  for r = 1:length (rectangles)
    % TODO: fix for DoW.
    [~, nr_on(:, r)] = not_ramp (on_off(2, r),    1, on_off(1,r):on_off(3,r)-1, cv);
    [~, nr_off(:,r)] = not_ramp (on_off(4, r)-1, -1, on_off(1,r):on_off(3,r)-1, cv);
    %nr(1, r) = not_ramp (on_off(2, r),    1, on_off(1,r):on_off(3,r)-1, cv);
    %nr(2, r) = not_ramp (on_off(4, r)-1, -1, on_off(1,r):on_off(3,r)-1, cv);
  end
  %nr = nr(1, :) .* nr(2, :);
  nr = dempster_shafer_2_combine (nr_on, nr_off);
  nr(1:2,:) = nr(1:2,:) / 2;
  nr(3,:) = 1 - (nr(1,:) + nr(2,:));

  chain_dists = zeros (size (rectangles));
  chain_dists([chains{:}]) = repelem (dists .^ 1.25, cellfun (@length, chains));
  chain_ds(3,:) = 0.5 * exp (-chain_value ./ chain_dists);
  chain_ds(1,:) = 1 - 2 * chain_ds(3,:);
  chain_ds(2,:) = chain_ds(3,:);
  
  Kolmogrov_Smirnoff_ds(1,:) = (1 - exp (log (ks_p(1,:)) / 20)); % * 0.75;
  Kolmogrov_Smirnoff_ds(3,:) = 1 - Kolmogrov_Smirnoff_ds(1,:);
  Kolmogrov_Smirnoff_ds(2,:) = 0;
  
  base = -rolling_min (-rolling_min (cv'))';
  [rectangles, certainty] = end_certainties (rectangles, cv, base);
  clear_edges = sum (certainty, 2)';
  clear_edges_ds(1,:) = clear_edges/2;
  clear_edges_ds(3,:) = (1 - clear_edges_ds(1,:)) / 2;
  clear_edges_ds(2,:) = clear_edges_ds(3,:);

%{
  % Times known to give false positives
  on_off = [rectangles.on_off];
  winter_evening(3,:) = (on_off(1,:) > 100 & on_off(3,:) < 300 ...
                        & on_off(2,:) > 35 & on_off(4,:) > 35 ...
                        & on_off(4,:) < 46) / 2;
  winter_evening(2,:) = winter_evening(3,:);
%}

%fprintf ('on_off chain, ks_p, nr\n');
%[on_off; -chain_value./chain_dists; log(ks_p(1,:))/20; log(1-nr); -clear_edges]'
%fprintf('\n');
  %qualities = (1 - exp ((-chain_value ./ chain_dists + log (ks_p(1,:))/20 + log (1-nr)) - clear_edges));
  qualities_ds = dempster_shafer_2_combine ( ...
                  dempster_shafer_2_combine (nr, chain_ds), ...
                  dempster_shafer_2_combine (Kolmogrov_Smirnoff_ds, ...
                                             clear_edges_ds) ...
                );
  qualities = qualities_ds(1,:) + 0.5 * qualities_ds(3,:);

  % weighted mean of power
  [pwrs, ia, ic] = unique ([rectangles.power]);
  probs = exp (-(pwrs - 0.7).^2 / 0.5);
  wgts = qualities .* probs(ic);
  wgts(isnan (wgts)) = 0;
  if sum (wgts) == 0
    wgts(:) = 1;
  end
  mean_pwr = [rectangles.power] * wgts' / sum (wgts);
  probs = -([rectangles.power] - mean_pwr).^2 / 0.3;
  probs = (probs - max (probs));
  power_ds(2,:) = sqrt (-probs);
  power_ds(3,:) = 1 - power_ds(2,:);
  power_ds(1,:) = 0;
  
  qualities_ds = dempster_shafer_2_combine (qualities_ds, power_ds);
  qualities = qualities_ds(1,:) + 0.5 * qualities_ds(3,:);
  %qualities = (1 - min (1, exp ((-chain_value ./ chain_dists + log (ks_p(1,:))/20 + log (1-nr) - probs ) - clear_edges)));

  figure(5); show_heatmap (qualities, rectangles, cv);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stub
% Should find weighted mean based on quality of each jump
%function power = common_power(rectangles, cv)
%    power = median([rectangles.power]);
%end

function [rectangles, details] = extend_altday (rectangles, cv, ad, fcn, power, issolar, all_on_off)
  % Fill alt-days with neighbouring days, then process with fcn.
  cv_ish = [cv(:,3), cv, cv(:,end-2)];
  r = ad(1):7:size(cv,2);
  cv_ish(:, r+1) = cv_ish (:, r);
  r = ad(2):7:size(cv,2);
  cv_ish(:, r+1) = cv_ish (:, r+2);
  [rectangles, details] = fcn (rectangles, cv_ish(:,2:end-1), power, issolar, all_on_off);
end


function rectangles = extend_rectangles(rectangles, cv, power, issolar, valid_days)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %figure (8); show_rectangles (rectangles, cv, []);
  rectangles = merge_neighbour_rects (rectangles, cv);

  % group rectangles into chains, representing periods of const timer settings
  pdfs.duration_ratio = @duration_ratio_pdf;
  pdfs.day_gap = @day_gap_pdf;
  pdfs.centre_change_gap = @centre_change_gap_pdf;
  pdfs.centre_change_no_gap= @centre_change_no_gap_pdf;
  % Re-assign trust to rectangles based on Kolmogorov-Smirnov test.
  rectangles = SNR_trust (rectangles, cv);
  if isempty (rectangles)
    return
  end
  idx = ([rectangles.ks_p] > 0.5);
  idx = ~ idx(1,:);

  [chains, dists, chain_values] = timed_chains (rectangles, pdfs, cv);
  [rectangles, transition, ~, change_on_off] = check_pairings (rectangles, chains, cv);
  if change_on_off
    [chains, dists, chain_values] = timed_chains (rectangles, pdfs, cv);
  end
  qualities = scalar_quality (rectangles, chains, dists, chain_values, cv, issolar);
  qq = mat2cell (qualities, 1, ones(1, length (qualities)));
  [rectangles.quality] = qq{:};
%   for j = 1:length (chains)
%     figure(100 + j);
%     show_rectangles (rectangles(chains{j}), cv);
%     fprintf ('%g\n', (dists(j)));
%   end
  % Delete rectangles that don't seem to fit into a chain.
  % TODO: check for matches between start/end dates
  % TODO: check individual rectangles within chains -- input to Dempster/Shafer?
  mask = true (size (rectangles));
  for j = 2:length (chains)
    if (dists(j) > 1000 && dists(j)/dists(j-1) > 3) ...
        || dists(j)/dists(j-1) > 50 || dists(j) > 1500

      %figure(100); show_rectangles (rectangles([chains{j:end}]), cv);
      mask([chains{j:end}]) = false;
      break
    end
  end
  thresh = zeros (size (rectangles));
  for j = 1:length (chains)
    %thresh(chains{j}) = 1 - 0.9^(dists(j) / 160);
    thresh(chains{j}) = 1 - exp (-0.8e-3 * (j-0.5)^2 * dists(j));
  end
%  mask = mask & (qualities > thresh);
  figure(8); show_rectangles (rectangles(~mask), cv);
  on_off = [rectangles.on_off];
  mask = mask & (on_off(3,:) > on_off(1,:));
  rectangles = rectangles(mask);

  days = size (cv, 2);
  on_off = [rectangles.on_off];

  ifp = isfield (rectangles, 'profile');
  ifa = isfield (rectangles, 'alt_days');
      % proper rectangles
  if ifp
    a = ([rectangles.profile] == 0);
%  elseif ifa
%    a = cellfun (@isempty, {rectangles.alt_days});
  else
    a = true (1, length (rectangles));
  end
    %a = cellfun (@eq, {rectangles.profile}, {0});
    %a = ([rectangles.profile] == 0);
  [rectangles(a), detail(a,:,:), new_r, new_d] = extend_flat (rectangles(a), cv, power, issolar, [rectangles.on_off]);
  rectangles = [rectangles, new_r];
  detail = [detail; new_d];
  a = [a, true(1, length (new_r))];
  b = a;

      % variable power
  if ifp
    %a = cellfun (@eq, {rectangles.profile}, {3});
    a = ([rectangles.profile] == 3);
    [rectangles(a), detail(a,:,:)] = extend_irregular (rectangles(a), cv, power, issolar);
    b = b | a;
  end

%   if ifa
%       % different on weekends
%     if ifp
%       %a = cellfun (@eq, {rectangles.profile}, {1});
%       a = ([rectangles.profile] == 4);
%     else
%       a = ~cellfun (@isempty, {rectangles.alt_days});
%     end
%     ad = sort ([rectangles.alt_days]);      % sort within pairs.
%     ar = find (a);
%       % if alt days are the same for all rectangles
%     while any (ar)
%       pair = sort (rectangles(ar(find (ar(1,:), 1))).alt_days);
%       match_pos = (ad(1,:) == pair(1) & ad(2,:) == pair(2));
%       match = ar(match_pos);
%
%       if (length (power) > 1)
%         p = power (match);
%       else
%         p = power;
%       end
%       if (all (pair == [1; 7]))
%         pair = [7; 1];
%       end
%
%       [rectangles(match), detail(match, :, :)] ...
%           = extend_altday (rectangles(match), cv, pair, @extend_flat, p, issolar, [rectangles.on_off]);
% %      % Replace weekends by the nearest weekday (Sat by Fri, Sun by Mon)
% %      cv_ish = [cv(:,3), cv, cv(:,end-2)];
% %      r = ad(pos,1):7:size(cv,2);
% %      cv_ish(:, r+1) = cv_ish (:, r);
% %      r = ad(pos,2):7:size(cv,2);
% %      cv_ish(:, r+1) = cv_ish (:, r+2);
% %      rectangles(aa) = extend_flat (rectangles(aa), cv_ish(:,2:end-1), power);
%       ar(match_pos) = 0;
%     end


%    if ~(any (diff (ad (1,:))) || any (diff (ad (2,:))))
%      % Replace weekends by the nearest weekday (Sat by Fri, Sun by Mon)
%      cv_ish = [cv(:,3), cv, cv(:,end-2)];
%      r = ad(1):7:size(cv,2);
%      cv_ish(:, r+1) = cv_ish (:, r);
%      r = ad(2):7:size(cv,2);
%      cv_ish(:, r+1) = cv_ish (:, r+2);
%      rectangles(a) = extend_flat (rectangles(a), cv_ish(:,2:end-1), power);
%    else
%      % What do we do now?
%      warning ('Misaligned alt days');
%      ar = a;
%      ad = sort (ad);
%      wrap = (ad(1,:) == 1 && ad(2,:) == 7);
%      ad(wrap,:) = repmat ([7 1], 1, length (wrap));
%      while any (ar)
%      %   Process all.
%        ar = ar & ~match;
%      end
%    end
%     b = b | a;
%
%     if ifp
%       a = ([rectangles.profile] == 7);
%       ad = [rectangles.alt_days];
%       if ~(any (diff (ad (1,:))) || any (diff (ad (2,:))))
%         % Replace weekends by the nearest weekday (Sat by Fri, Sun by Mon)
%         cv_ish = [cv(:,3), cv, cv(:,end-2)];
%         r = ad(1):7:size(cv,2);
%         cv_ish(:, r+1) = cv_ish (:, r);
%         r = ad(2):7:size(cv,2);
%         cv_ish(:, r+1) = cv_ish (:, r+2);
%         [rectangles(a), detail(a,:,:)] ...
%                = extend_irregular (rectangles(a), cv_ish(:,2:end-1), power);
%       end
%       b = b | a;
%     end
%   end

  if any (~b)
    warning ('Missed %d rectangles', sum (~b));
  end

  % TODO:  Delete rectangles with long stretches below power.
  %        Save their edges, and overwrite min(detail) if max(detail) is high.
  on_off = [rectangles.on_off];
  detail_sup = detail;
  for i = 1:length (rectangles)
    slice = cv (rectangles(i).burst(1:end-1), on_off(1,i):on_off(3,i)-1);
    % TODO: Handle alt_days.
    rectangles(i).missed = on_off(1,i) - 1 + find (min(slice,[],1) < 0.7 * rectangles(i).power);
    if length (rectangles(i).missed) > 0.6 * (on_off(3,i) - on_off(1,i))
      new_r = find_rect  (cv, on_off(1,i), on_off(3,i), ...
                                   rectangles(i).power, i, rectangles);
      if ~isempty (new_r) && ~any (abs (on_off(2,:) - new_r.on_off(2)) < 1 ...
                                 & abs (on_off(4,:) - new_r.on_off(4)) < 1 ...
                                 & on_off(1,:) < new_r.on_off(3) ...
                                 & on_off(3,:) > new_r.on_off(1) ...
                                 & 1:size(on_off, 2) ~= i)
          rectangles(i) = new_r;
          details(i,:,:) = top_and_bottom (new_r, cv, rectangles(i).power);
      elseif rectangles(i).power < 0.8   % > 0.8 too demanding sometimes
        detail_sup(abs (on_off(2,:) - on_off(2,i)) < 0.5, on_off(1,i):on_off(3,i)-1, 2) = 1;
        detail_sup(abs (on_off(4,:) - on_off(4,i)) < 0.5, on_off(1,i):on_off(3,i)-1, 1) = 1;
        on_off (1,i) = -1;
      end
    end
  end
  % Mark zero-length rectangles for deletion
  on_off (1, on_off(1,:) == on_off(3,:)) = -1;

  rectangles = rectangles (on_off(1,:) >= 0);
  if isempty (rectangles)
    return
  end
  
  detail_sup = detail_sup (on_off(1,:) >= 0, :, :);
  detail = detail (on_off(1,:) >= 0, :, :);

  % Remove exactly overlapping rectangles.
  on_off = [rectangles.on_off];
  overlap_all = bsxfun (@eq, on_off(1,:), on_off(1,:)') ...
              & bsxfun (@eq, on_off(3,:), on_off(3,:)') ...
              & abs (bsxfun (@minus, on_off(2,:), on_off(2,:)')) < 0.5 ...
              & abs (bsxfun (@minus, on_off(4,:), on_off(4,:)')) < 0.5;
      % remove one of mutually overlapping pair
  mutual = triu (overlap_all, 1);
  overlap_r = any (mutual);
  rectangles = rectangles (~overlap_r);
  detail_sup = detail_sup (~overlap_r, :, :);
  detail = detail (~overlap_r, :, :);

  % Re-assign trust to rectangles based on Kolmogorov-Smirnov test.
  rectangles = SNR_trust (rectangles, cv);
  if isempty (rectangles)
    return
  end
  idx = ([rectangles.ks_p] > 0.5);
  idx = ~ idx(1,:);

  on_off = [rectangles.on_off];
  %figure(11); show_on_off (on_off(:, ~idx), cv);

  rectangles = rectangles(idx);
  detail_sup = detail_sup(idx, :, :);
  detail = detail(idx, :, :);

  pdfs.duration_ratio = @duration_ratio_pdf;
  pdfs.day_gap = @day_gap_pdf;
  pdfs.centre_change_gap = @centre_change_gap_pdf;
  pdfs.centre_change_no_gap= @centre_change_no_gap_pdf;
  % Re-assign trust to rectangles based on Kolmogorov-Smirnov test.
  rectangles = SNR_trust (rectangles, cv);
  if isempty (rectangles)
    return
  end
  idx = ([rectangles.ks_p] > 0.5);
  idx = ~ idx(1,:);

  [chains, dists, chain_values] = timed_chains (rectangles, pdfs, cv);
  qualities = scalar_quality (rectangles, chains, dists, chain_values, cv, issolar);
  qq = mat2cell (qualities, 1, ones(1, length (qualities)));
  [rectangles.quality] = qq{:};
%   for j = 1:length (chains)
%     figure(100 + j);
%     show_rectangles (rectangles(chains{j}), cv);
%     fprintf ('%g\n', (dists(j)));
%   end
  % Delete rectangles that don't seem to fit into a chain.
  % TODO: check for matches between start/end dates
  % TODO: check individual rectangles within chains -- input to Dempster/Shafer?
  mask = true (size (rectangles));
%   for j = 2:length (chains)
%     if (dists(j) > 1000 && dists(j)/dists(j-1) > 3) ...
%         || dists(j)/dists(j-1) > 50 || dists(j) > 1500
%
%       %figure(100); show_rectangles (rectangles([chains{j:end}]), cv);
%       mask([chains{j:end}]) = false;
%       break
%     end
%   end
  thresh = zeros (size (rectangles));
  for j = 1:length (chains)
    %thresh(chains{j}) = 1 - 0.9^(dists(j) / 160);
    thresh(chains{j}) = 1 - exp (-0.8e-3 * (j-0.5)^2 * dists(j));
  end
  mask = mask & (qualities > thresh);
  figure(8); show_rectangles (rectangles(~mask), cv);
  rectangles = rectangles(mask);
  detail = detail(mask, :, :);
  detail_sup = detail_sup(mask, :, :);



%   % group rectangles into chains, representing periods of const timer settings
%   pdfs.duration_ratio = @duration_ratio_pdf;
%   pdfs.day_gap = @day_gap_pdf;
%   pdfs.centre_change_gap = @centre_change_gap_pdf;
%   pdfs.centre_change_no_gap= @centre_change_no_gap_pdf;
%   [chains, dists] = timed_chains (rectangles, pdfs, cv);
% %   for j = 1:length (chains)
% %     figure(100 + j);
% %     show_rectangles (rectangles(chains{j}), cv);
% %     fprintf ('%g\n', (dists(j)));
% %   end
%   % Delete rectangles that don't seem to fit into a chain.
%   % TODO: check for matches between start/end dates
%   % TODO: check individual rectangles within chains -- input to Dempster/Shafer?
%   for j = 2:length (chains)
%     if (dists(j) > 900 && dists(j)/dists(j-1) > 3) ...
%         || dists(j)/dists(j-1) > 50 || dists(j) + 100*j > 1600
% %       figure(100); show_rectangles (rectangles([chains{j:end}]), cv);
% %       rectangles = rectangles ([chains{1:j-1}]);
% %       detail = detail([chains{1:j-1}],:,:);
% %       detail_sup = detail_sup([chains{1:j-1}],:,:);
%       break
%     end
%   end


  % Look for "new" rectangles, with the same on/off times as existing ones.
  % - both on and off are positive, for more than 20 days.
  % - evidence: on and off cross zero at nearly the same time
  % - evidence: interval between this rectangle and the original is "noisy"
  % - evidence: on and off are much larger than 0
  % - evidence: no overlapping rectangle, especially of higher "quality"
  % - evidence: new rectangle is close to the existing one that triggered it.
  on_off = [rectangles.on_off];
  %dtl = min (detail, [], 3);   % len(rectangles)xsize(cv,2), min of on and off jumps
  dtl = max (max (min (detail(:,:,1), detail_sup(:,:,2)), ...
                  min (detail(:,:,2), detail_sup(:,:,1))), [], 3);
  new_rectangles = rectangles; % arbitrary initialization with right fields
  new_count = 0;
  alt_hrs = isfield (rectangles, 'alt_hrs');
  detail_idx = [];
  cvw = [cv; cv(:,2:end), cv(:,end)];        % wrap

  for i = 1 : size (dtl, 1)
    if rectangles(i).quality < 0.35
      continue;
    end
    if rectangles(i).ks_p(1) > 0.45
%      continue;
foo = new_count;
    end
    a = ranges (find (dtl (i,:) > 0));
    if isempty (a)
      continue;
    end

%{
      % Truncate any overlapping with the rectangle that triggered this.
    overlap = (a(1,:) < on_off(3,i) & a(2,:) > on_off(1,i));
    a(1, overlap & a(1,:) >= on_off(1,i)) = on_off(3,i);
    a(2, overlap & a(2,:) <= on_off(3,i)) = on_off(1,i);
    % TODO Check whether original or overlap is more plausible,
    %      for start and end.
%}

    remove = (a(2,:) - a(1,:) < 6) & a(2,:) < size (cv, 2) & a(1,:) > 1;
    a = a (:, ~remove);
    if isempty (a)
      continue;
    end

    hr_olap = @(x) ~isempty (intersect (x, rectangles(i).burst));
    for j = 1:size (a,2);
      % First, discard dubious candidates.
      dist = min (abs ([a(1,j) - on_off(3,i), a(2,j) - on_off(1,i)]));
      olap = a(1,j) < on_off(3,:) & a(2,j) > on_off(1,:) ...
             & cellfun (hr_olap, {rectangles(:).burst});
      if (a(1,j) == 1) || a(2,j) == size (cv, 2)
        len = 20;
      else
        len = a(2,j) - a(1,j);
      end

      % If two of: this candidate overlaps, is too far away, or too small,
      % then treat it as spurious
      if (any (olap)) + (dist > 70) + (len < 20) >= 2
          continue;
      end

      r = [a(1,j); on_off(2,i); a(2,j); on_off(4,i)];
      r = split_at_overlap (r, on_off, cv, cvw, r(1), i, rectangles, true, issolar);
      if isempty (r) || isequal (r, on_off(:,i))
        continue;
      end
      r(2,:) = r(2,1);
      r(4,:) = r(4,1);

      for n = 1:size(r, 2)
        % Next, add reputable candidates.
        new_count = new_count + 1;
        new_rectangles(new_count).on_off = r(:,n);
        new_rectangles(new_count).missed = [];
        new_rectangles(new_count).burst = rectangles(i).burst;
        new_rectangles(new_count).trust = 1;
        new_rectangles(new_count).quality = 0;
        new_rectangles(new_count).power = rectangles(i).power;
        if alt_hrs
          new_rectangles(new_count).alt_hrs = rectangles(i).alt_hrs;
          new_rectangles(new_count).alt_days = rectangles(i).alt_days;
        end
        detail_idx (new_count) = i;
      end
    end
if rectangles(i).ks_p(1) > 0.45 && new_count > foo
  oo = new_rectangles(foo+1:new_count);
  %figure(11); show_on_off ([oo.on_off], cv);
  new_count = foo;
end
  end
  new_rectangles = end_days (new_rectangles(1:new_count), cv);
  spike_ok = false;
  new_rectangles = assign_trust (new_rectangles, cv, power, issolar, spike_ok);
  rectangles = [rectangles, new_rectangles];
  detail_idx = [1:size(dtl,1), detail_idx];

prev = [];
count = 0;
cv_rev = cv(:, end:-1:1);
detail_rev = detail(:, end:-1:1, :);
while (~isequal (prev, [rectangles.on_off]) && count < 4)
  count = count + 1;
  % Identify overlapping rectangles
      % First, exactly overlapping rectangles.
  on_off = [rectangles.on_off];
  prev = on_off;
  overlap_all = bsxfun (@eq, on_off(1,:), on_off(1,:)') ...
              & bsxfun (@eq, on_off(3,:), on_off(3,:)') ...
              & abs (bsxfun (@minus, on_off(2,:), on_off(2,:)')) < 0.5 ...
              & abs (bsxfun (@minus, on_off(4,:), on_off(4,:)')) < 0.5;
      % remove one of mutually overlapping pair
  mutual = triu (overlap_all, 1);
  overlap_r = any (mutual);
  rectangles = rectangles (~overlap_r);
  detail_idx = detail_idx (~overlap_r);

  % Next, partially overlapped.
  on_off = round ([rectangles.on_off]);

  [on_off, r] = align_pair_edge (rectangles, cv, on_off, detail, detail_idx, transition);
  r_oo = round ([rectangles.on_off]);
  for i = find (r_oo(2,:) ~= on_off(2,:) | r_oo(4,:) ~= on_off(4,:))
    rectangles(i).on_off = r(i).on_off;
  end

  % Perform deletion
  rectangles = rectangles(on_off(2,:) >= 0 & on_off(3,:) > on_off(1,:));
  detail_idx = detail_idx(on_off(2,:) >= 0 & on_off(3,:) > on_off(1,:));
  on_off     = on_off(:,  on_off(2,:) >= 0 & on_off(3,:) > on_off(1,:));

  % Repeat for left edges, by reflecting cv and on_off
  on_off([1,3],:)   = size (cv, 2) + 2 - on_off([3,1],:);

  [on_off, r] = align_pair_edge (rectangles, cv_rev, on_off, detail_rev, detail_idx, transition(end:-1:1));
  r_oo = round ([rectangles.on_off]);
  for i = find (r_oo(2,:) ~= on_off(2,:) | r_oo(4,:) ~= on_off(4,:))
    rectangles(i).on_off([2,4]) = r(i).on_off([2,4]);
  end

  rectangles = rectangles(on_off(2,:) >= 0 & on_off(3,:) > on_off(1,:));
  detail_idx = detail_idx(on_off(2,:) >= 0 & on_off(3,:) > on_off(1,:));
  on_off     = on_off(:,  on_off(2,:) >= 0 & on_off(3,:) > on_off(1,:));

  on_off([1,3],:)   = size (cv, 2) + 2 - on_off([3,1],:);
  oo = [rectangles.on_off];
  for i = find (any (oo([1,3],:) ~= on_off([1,3],:)))
    rectangles(i).on_off([1,3]) = on_off([1,3], i);
  end

  % Again, exactly overlapping rectangles.
  on_off = [rectangles.on_off];
  overlap_all = bsxfun (@eq, on_off(1,:), on_off(1,:)') ...
              & bsxfun (@eq, on_off(3,:), on_off(3,:)') ...
              & abs (bsxfun (@minus, on_off(2,:), on_off(2,:)')) < 0.5 ...
              & abs (bsxfun (@minus, on_off(4,:), on_off(4,:)')) < 0.5;
      % remove one of mutually overlapping pair
  mutual = triu (overlap_all, 1);
  overlap_r = any (mutual);
  rectangles = rectangles (~overlap_r);
  detail_idx = detail_idx (~overlap_r);
end
  % Re-assign trust to rectangles based on Kolmogorov-Smirnov test.
  if ~isempty (rectangles)
    rectangles = SNR_trust (rectangles, cv);
    if ~isempty (rectangles)
      idx =  ([rectangles.ks_p] > 0.005);
      idx = ~ idx (1,:);
      %figure(7); show_rectangles (rectangles (~idx), cv, []);
      rectangles = rectangles(idx);
    end
  end

%figure(12); show_rectangles (rectangles, cv, []);
return; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Find pairs of rectangles such that the day-range of one is
  % entirely within that of the other
  % and the hour ranges overlap.
  on_off = [rectangles.on_off];
  overlap_part = bsxfun (@ge, on_off(1,:)', on_off(1,:) - 2) ...
               & bsxfun (@le, on_off(3,:)', on_off(3,:) + 2) ...
               & bsxfun (@lt, on_off(2,:)', on_off(4,:)) ...
               & bsxfun (@gt, on_off(4,:)', on_off(2,:));
  % Don't double count multual overlaps, or count overlap with itself.
  overlap_part = overlap_part & ~(triu (overlap_part))';

  clear ol;
  [ol(1,:), ol(2,:)] = find (overlap_part);
  used = (0 ~= set_cor_pool (rectangles, cv', 1:size (cv, 2)))';
  on_off = [rectangles(ol(:)).on_off];
  for i = 1:size(ol,2)
    show_rectangles (rectangles ([ol(1,i), ol(2,i)]), cv);
    0;
  end
  for i = 1:size(ol, 2)
    % look for expected jumps (based on cor_pools)
      %top
    pos = 2*i-1;
    left = vert_diff (used, cv, on_off(2, pos), on_off(pos,1), on_off(pos,3));
    right= vert_diff (used, cv, on_off(4, pos), on_off(pos,1), on_off(pos,3));
    top = hor_diff (used, cv, on_off(1, pos), rectangles(ol(1, i)).burst);
    bot = hor_diff (used, cv, on_off(3, pos), rectangles(ol(1, i)).burst);

    pos = pos + 1;
    l2 = vert_diff (used, cv, on_off(2, pos), on_off(pos,1), on_off(pos,3));
    r2= vert_diff (used, cv, on_off(4, pos), on_off(pos,1), on_off(pos,3));
    t2 = hor_diff (used, cv, on_off(1, pos), rectangles(ol(2, i)).burst);
    b2 = hor_diff (used, cv, on_off(3, pos), rectangles(ol(2, i)).burst);
    % Choose rectangle with most "missing" jumps,
    % or least relible jumps, or something...
    % Give extra weight to an edge that follows on from another rectangle
  end


  % Recalculate transitions of timer on/off points.
  % Types of evidence:
  % - multiple timer transitions at the same day
  % - multiple "on" periods shift at the same day
  % - "on" and "off" times shift by the same amount
  % - low variance in surrounding region / within rectangle, and crips drop

  % Steps:
  % For all overlaps in decreasing order
  %   if either rectangle overlaps another on the same edge
  %     continue
  %   if there is one edge to move
  %     if outer rectangle abuts another

  % Remove totally overlapped rectangles.
  % TODO: Check reliability first!
  %figure(11); show_rectangles (rectangles(on_off(1,1:end-1) >= 0),cv, []);
  overlap_all = bsxfun (@ge, on_off(1,:), on_off(1,:)') ...
              & bsxfun (@le, on_off(3,:), on_off(3,:)') ...
              & bsxfun (@ge, on_off(2,:), on_off(2,:)') ...
              & bsxfun (@le, on_off(4,:), on_off(4,:)');
      % don't remove both of mutually overlapping pair
  mutual = triu (overlap_all & overlap_all');
  overlap_r = any (overlap_all & ~mutual);
  rectangles = rectangles (~overlap_r);
  detail_idx = detail_idx (~overlap_r);
  %figure(12); show_rectangles (rectangles, cv, []);
  i;

  for r = 1:length(rectangles)
      % If no next rectangle, or non-overlapping next rectangle, or
      % "unreliable" next rectangle, check if edges extend
  end
end

function rectangles = merge_neighbour_rects (rectangles, cv)
  % Reconcile overlaps
  on_off = round ([rectangles.on_off]);
  pwr = mean ([rectangles.power]);
  d1 = reshape (diff ([cv(end); cv(:)]), size (cv));
  [ol, overlap_hours] = find_overlap (on_off, rectangles);
  for i = 1:size(ol, 2)
    a = ol(1,i);
    b = ol(2,i);
    c = [on_off(1,a), on_off(1,b), on_off(3,a), on_off(3,b)];
    if c(1) > c(2) || (c(1) == c(2) && c(3) > c(4))
      [a, b] = deal (b, a);
      c = [on_off(1,a), on_off(1,b), on_off(3,a), on_off(3,b)];
    end
    d = diff (c);
    if all (d >= 0) && overlap_hours (a, b) > 1
      p = [rectangles([a,b]).power];
      on_off = find_switch_day (on_off, a, b, cv, d1, p, pwr, rectangles);
    end
  end
  for i = find (any (on_off ~= [rectangles.on_off], 1))
    rectangles(i).on_off([1,3]) = on_off([1,3], i);
  end

  % Delete empty rectangles
  idx = (on_off(1,:) ~= on_off(3,:));
  rectangles = rectangles(idx);

  % merge contiguous neighbours
  on_off = [rectangles.on_off];
  a13 = bsxfun (@eq, on_off(1,:), on_off(3,:)');
  a2 = abs (bsxfun (@minus, on_off(2,:), on_off(2,:)')) < 0.5;
  a4 = abs (bsxfun (@minus, on_off(4,:), on_off(4,:)')) < 0.5;
  [x, y] = find (a13 & a2 & a4);
  [~, idx] = sort (-on_off(3,y));
  for i = idx(:)'
    if isfield (rectangles, 'alt_days')
      if isempty (rectangles(x(i)).alt_days)
        rectangles(x(i)).alt_days = rectangles(y(i)).alt_days;
        rectangles(x(i)).alt_hrs = rectangles(y(i)).alt_hrs;
      elseif ~isequal (rectangles(x(i)).alt_days, rectangles(y(i)).alt_days)
        a13(x,y) = false;   % don't delete y
        continue;           % can't merge incompatible alt-days
      else
        % TODO: find new alt_hrs.  For now, just take those of x.
      end
    end
    rectangles(x(i)).on_off(3) = rectangles(y(i)).on_off(3);
  end
  rectangles = rectangles (~any (a13 & a2 & a4, 1));
end

function p = duration_ratio_pdf (ratios)
  % -log (Probability distribution function) for the ratio of durations
  % of consecutive rectangles in a chain,
  % estimated from rectangles manually selected using pool_ground_truth
  idx = (ratios > 1.03);
  p = ratios;
  % TODO: Fix this model!
  p(idx) = 3.733 * log (0.1558 * ratios(idx));
  p(~idx) = 14 * log (ratios(~idx));
end

function p = centre_change_gap_pdf (changes)
  % -log (Probability distribution function) for the change in the midpoint
  % of the on period of timed devices of consecutive rectangles in a chain
  % where there is a period of a few days between them,
  % estimated from rectangles manually selected using pool_ground_truth
  b = 3.6;
  p = abs (changes) * b - log (2 * b);     % TODO: Model spike at 0
end

function p = centre_change_no_gap_pdf (changes)
  % -log (Probability distribution function) for the change in the midpoint
  % of the on period of timed devices of consecutive rectangles in a chain
  % where the second follows immediately after the first,
  % estimated from rectangles manually selected using pool_ground_truth
  b = 2.8;
  p = abs (changes) * b - log (2 * b);     % TODO: Model spikes
end

function p = day_gap_pdf (days)
%   exp_a = -0.6554;
%   b = -0.4373;
%   c_log_days = -0.1131 * log (1 + days(days >= 0));
%   p(days >= 0) = exp_a * (b + 2 * c_log_days) ...
%                  .* (1 + days(days >= 0)) .^ (b -  1 + c_log_days);
%   p(days < 0) = -exp_a * exp (-days(days < 0) .^ 2);
  p(days >= 0) = 5 * days(days >= 0);
  p(days <  0) = 5 * (-days(days < 0)) .^ 1.5;
end

function [rectangles, transition, chain_count, change_on_off] = check_pairings (rectangles, chains, cv)
  % Identify rectangles that start/stop at similar times in
  % different chains.
  % Guess how many chains there "really" are.
  % Reduce "quality" for all at same time as, but not paired with, higher
  % quality rectangles.
  % Split rectangles if the pair is split
  % See also code in timed_chains.
  
  % Evidence:
  % - Similar duration
  % - Start and/or end matches
  % - Jump same as / opposite to each other
  % - Chains match
  
  % Determine certainty about start/end day of each rectangle
  % - Look for hour with noise << rectangle power, and no neighbour
  
  % Find start pairs
  % Find end pairs
  days = size (cv, 2);
  hrs  = size (cv, 1);
  change_on_off = false;
  chain_count = length (chains);
  base = -rolling_min (-rolling_min (cv'))';
  [rectangles, certainty] = end_certainties (rectangles, cv, base);
  
  on_off = [rectangles.on_off];
  transition = zeros (days, 1);
  transition(on_off(1,:)) = 1;
  transition(on_off(3,:)) = 1;
  
  rev_map (transition ~= 0) = 1:sum (transition);
  pairings = (abs (bsxfun (@minus, find (transition), ...
                           [on_off(1,:), on_off(3,:)])) < 3);
  transition(transition > 0) = sum (pairings, 2);
  pairs = [transition(on_off(1,:)); transition(on_off(3,:))];
  pairs = mat2cell (pairs, ones(1, size (pairs, 1)));
  [rectangles.pairings] = pairs{:};

  jump_before = nan (size (rectangles));
  jump_after  = nan (size (rectangles));
  rev_chains = zeros (size (rectangles));  % chain each rect belongs to
  for c = 1:length (chains)
    ch = chains{c};
    rev_chains(chains{c}) = c;
    jumps = on_off([2,4], ch(2:end))-on_off([2,4],ch(1:end-1));
    jumps = mod (jumps + hrs/2, hrs) - hrs/2;
    [~, pos] = max (abs (jumps));
    jumps = jumps (sub2ind (size(jumps), pos, 1:size(jumps, 2)));
    jump_before(ch(2:end))   = jumps;
    jump_after (ch(1:end-1)) = jumps;
  end
  
  % Form initial pairings based on matched start/end dates
  start_quality = cell (size (rectangles));
  end_quality = start_quality;
  start_match = start_quality;
  end_match   = start_quality;
  for r = 1:length (rectangles)
    sm = (pairings (rev_map (on_off(1,r)), 1:end/2));
    em = (pairings (rev_map (on_off(3,r)), end/2+1:end));
    sm(r) = false;
    em(r) = false;
    start_match{r} = find (sm);
    end_match{r}   = find (em);
    
    if ~isempty (start_match)
      start_quality{r} = similarity (r, start_match{r}, rectangles, jump_before, chains, hrs, true);
    end
    if ~isempty (end_match)
      end_quality{r} = similarity (r, end_match{r}, rectangles, jump_after, chains, hrs, false);
    end
  end
  
  % Refine pairings to omit three-way pairings,
  % include overlaps with mismatched ends,
  % and runs of adjacent rectangles.
  % Step 1. Find usual time interval between pairs
  single_start = cellfun (@length, start_match) == 1;
  single_end   = cellfun (@length, end_match)   == 1;
  gaps = mod (on_off(2, [find(single_start), find(single_end)]) ...
            - on_off(2, [start_match{single_start}, ...
                         end_match{single_end}]) + hrs/2, hrs) - hrs/2;
  sd_gaps = sqrt (var (abs (gaps)));
  gaps = sum ([start_quality{single_start}, end_quality{single_end}] .* abs (gaps)) ...
       / sum ([start_quality{single_start}, end_quality{single_end}]);
  
  % Step 2. omit three-way pairings
  multi_start = cellfun (@length, start_match) > 1;
  multi_end   = cellfun (@length, end_match)   > 1;
  
  % Step 3. find overlaps with mismatched edges.
  no_start = find (cellfun (@isempty, start_match) & certainty(:,1)' > 0.5);
  cvw = [cv; cv(:,2:end), cv(:,end)];   % wrap
  for r = no_start(:)'
    % Look for rectangles overlapping the start, about gap away.
    overlap_hr = (on_off(1,:) < on_off(1,r) & on_off(3,:) > on_off(1,r));
    offset = on_off(2, overlap_hr) - on_off(2, r);
    offset = abs (mod (hrs/2 + offset, hrs) - hrs/2) - gaps;
    ok = (abs (offset) < 2 * sd_gaps);
    if any (ok)
      pair = find (overlap_hr);
      pair = pair (ok);
      % TODO: consider all possible pairs.
      % TODO: consider rectangles crossing midnight
      pair = pair(1);
      pre_days = on_off(1, pair):on_off(1,r)-1;
      post_days = on_off(1,r):min (on_off(1,r)+5, on_off(3,r)-1);
      burst = rectangles(pair).burst;
      pre  = max (base(burst, pre_days),  [], 2);
      post = min (base(burst, post_days), [], 2);
      if sum (post(1:end-1) > pre(1:end-1)) > length (burst) / 2
        possible_starts(2) = on_off(1,r);
        if on_off(2,r) < on_off(2, pair)
          st = max (on_off(1,r) - 2, 1);
        else
          st = on_off(1,r)-1;
        end
        jumps = diff (mean (cv(burst, st:st+2)));
        [jump, pos] = max (jumps);
        if jump > rectangles(r).power / 2
          transition(on_off(1,pair)) = transition(on_off(1,pair)) - 1;
          transition(st + pos) = transition(st + pos) + 1;
          change_on_off = true;
          on_off(1,pair) = st + pos;
          rectangles(pair).on_off(1) = st + pos;
        end
      end
    end
    % See if it has a jump near my start, and isn't already paired
    % If so, increase its start time.
  end
  
  % Step 3a. TODO: Repeat for ends
  
  % Step 4. See if we should subdivide long runs
  singles = find (single_start & single_end);
  long =  ([start_match{singles}] ~= [end_match{singles}]);
  long = singles(long);
  new_rect = [];
  for r = long(:)'
    % For each rectangle in the other chain
    ch1 = rev_chains(start_match{r});
    ch2 = rev_chains(end_match{r});
    if ch1 == ch2
      ch = chains{ch1};
      pos1 = find (ch == start_match{r});
      pos2 = find (ch == end_match{r});
      for i = pos1+1 : pos2-1
        % Find deviation of on and off times
        offset1 = on_off([2,4], ch(i)) - on_off([2,4], start_match{r});
        offset2 = on_off([2,4], ch(i)) - on_off([2,4], end_match{r});
        offset = unique_ish (on_off(2, r) ...
                             + [0; offset1(:); offset2(:); ...
                               -offset1(:); -offset2(:); 0]', ...
                             size (cv, 1));
        % eval_edges for on/off +/- each difference and no difference
        dys = get_days (rectangles(ch(i)), on_off(:,ch(i)), false);
        [t, badness] = eval_edges  (offset, ...
                                    dys, on_off(4,r), ...
                                    [on_off(1,ch(i)), on_off(3,ch(i))], ...
                                    rectangles(ch(i)).power, ...
                                    rectangles(ch(i)).power, cvw, ...
                                    [], true, false);
        [~, pos] = min (badness);
        if pos > 0 && badness(pos) < badness(1) / 2
          turn_on = t(pos);
        else
          turn_on = on_off(2,r);
        end
        
        offset = unique_ish (on_off(4, r) ...
                             + [0; offset1(:); offset2(:); ...
                               -offset1(:); -offset2(:); 0]', ...
                             size (cv, 1));
        [t, badness] = eval_edges  (offset, ...
                                    dys, turn_on, ...
                                    [on_off(1,ch(i)), on_off(3,ch(i))], ...
                                    rectangles(ch(i)).power, ...
                                    rectangles(ch(i)).power, cvw, ...
                                    [], false, false);
        [~, pos] = min (badness);
        if pos > 0 && badness(pos) < badness(1) / 2
          turn_off = t(pos);
        else
          turn_off = on_off(4,r);
        end
        
        % split rectangles(r) into three parts:
        %   before the start of rectangles(ch(i)) -> new_2
        %   the same days as rectangles(ch(i))    -> new_1
        %   after the start of rectangles(ch(i))  -> rectangles(r)
        if turn_on ~= on_off(2,r) || turn_off ~= on_off(4,r)
          change_on_off = true;
          new_1 = rectangles(ch(i));
          new_1.on_off([2,4]) = [turn_on, turn_off];
          if new_1.on_off(1) > rectangles(r).on_off(1)
            new_2 = rectangles(r);
            new_2.on_off(3) = new_1.on_off(1);
          else
            new_2 = [];
          end
          on_off(1,r) = new_1.on_off(3);
          rectangles(r).on_off(1) = on_off(1,r);
          new_rect = [new_rect, new_1, new_2];
          % split rectangle.
        end
      end
    end
    rectangles = [rectangles, new_rect];
    
    % Consider merging the matching short runs
  end
  

  % Recalculate pairings if rectangles have changed.
  if change_on_off
    on_off = [rectangles.on_off];
    transition = zeros (days, 1);
    transition(on_off(1,:)) = 1;
    transition(on_off(3,:)) = 1;

    rev_map (transition ~= 0) = 1:sum (transition);
    pairings = (abs (bsxfun (@minus, find (transition), ...
                             [on_off(1,:), on_off(3,:)])) < 3);
    transition(transition > 0) = sum (pairings, 2);
    pairs = [transition(on_off(1,:)); transition(on_off(3,:))];
    pairs = mat2cell (pairs, ones(1, size (pairs, 1)));
    [rectangles.pairings] = pairs{:};
  end

  % Evaluate quality of pairings
  % If start pairing differs from end pairing
  %   Look for break in the longer rectangle; repeat
  %   Consider merging shorter
  %   Consider tweaking start/end time
  
  % Consequences:
  % - Don't extend beyond neighbour in chain
  % - Increased "quality".
end

function s = similarity (a, b, rectangles, jumps, chains, hrs, start)
  on_off = [rectangles(a).on_off, rectangles(b).on_off];
  duration = mod1 (on_off(2,:) - on_off(4,:), hrs);
  diff_len = abs (duration(1) - duration(2:end)) ./ (duration(1) + duration(2:end));
  diff_jump = min (abs ([jumps(a)+jumps(b), jumps(a)-jumps(b)])) ...
              ./ (abs (jumps(a)) + abs (jumps(b)));
  % TODO: Dempster-Shafer?
  s = 1 - mean ([diff_len; diff_jump], 'omitnan');
end

function [rectangles, certainty] = end_certainties (rectangles, cv, base)
  certainty = zeros (length (rectangles), 4);
  on_off = [rectangles.on_off];
  gap = 15;
  mini_gap = 4;
  left  = max (on_off(1, :) - gap, 1);
  right = min (on_off(3, :) + gap, size (cv, 2));
  
  for r = 1:length (rectangles)
    burst = rectangles(r).burst;
    low_burst = [];
    if any (burst > on_off(4,r))
      low_burst = burst(burst < on_off(2, r));
      burst = burst(burst > on_off(4, r) + 1);
    end

    % adjust on_off([1,3],r) if there is a real cliff.
    st = max (left(r), on_off(1,r)-mini_gap);
    en = min (on_off(1,r)+mini_gap, right(r)-1);
    lhs = cummax ([base(burst, st:en);    base(low_burst, st+1:en+1)],    2);
    mid = cummin ([base(burst, en:-1:st); base(low_burst, en+1:-1:st+1)], 2);
    [step, pos] = max (mid(:, end-1:-1:1) - lhs(:, 1:end-1), [], 2);
    h = hist (pos, 1 : (en - st + 1));
    h = find (h > sum(h)/2);
    if ~isempty (h) && mean (step (pos == h)) > rectangles(r).power / 2
      en = st + h;
      if en < rectangles(r).on_off(3)
        on_off(1, r) = en;
        rectangles(r).on_off(1) = en;
      end
    end

    st = max (left(r), on_off(3,r)-mini_gap);
    en = min (on_off(3,r)+mini_gap, right(r)-1);
    rhs = cummax ([base(burst, en:-1:st); base(low_burst, en+1:-1:st+1)], 2);
    mid = cummin ([base(burst, st:en);    base(low_burst, st+1:en+1)],    2);
    [step, pos] = max (mid(:, 1:end-1) - rhs(:, end-1:-1:1), [], 2);
    h = hist (pos, 1 : (en - st + 1));
    h = find (h > sum(h)/2);
    if ~isempty (h) && mean (step (pos == h)) > rectangles(r).power / 2
      st = st + h;
      if st > rectangles(r).on_off(1)
        on_off(3, r) = st;
        rectangles(r).on_off(3) = st;
      end
    end

    % Coarse test of certainty of dates.
    %     Better to test corners?
    %     Check there isn't an underlying shift in power on those dates
    lhs = base(burst, on_off(1,r)-1:-1:left(r));
    rhs = base(burst, on_off(3,r):right(r));
    mid = base(burst, on_off(1,r):on_off(3,r)-1);
    if ~isempty (low_burst)
      lhs = [lhs; base(low_burst, on_off(1,r):-1:left(r)+1)];
      if right(r) >= on_off(3,r)
        rhs = [rhs; base(low_burst, [on_off(3,r)+1:right(r), right(r)])];
      else
        rhs = zeros (length (burst) + length (low_burst), 0);
      end
      mid = [mid; base(low_burst, [on_off(1,r)+1:on_off(3,r)-1, on_off(3,r)-1])];
    end
    cm_lhs = cummax (lhs, 2);
    cm_rhs = cummax (rhs, 2);
    m_mid  = min (mid, [], 2);
    if any (max (cm_lhs, [], 2) < m_mid - rectangles(r).power / 4)
      ml = mean (lhs(:));
      mm = mean (mid(:));
      if ml < mm / 2
        certainty(r, 1) = 1;
      else
        certainty(r, 1) = max (0, 0.7 - ml / mm);
      end
    end
    if any (max (cm_rhs, [], 2) < m_mid - rectangles(r).power / 4)
      mr = mean (rhs(:));
      mm = mean (mid(:));
      if mr < mm / 2
        certainty(r, 3) = 1;
      else
        certainty(r,3) = max (0, 0.7 - mr / mm);
      end
    end
    % TODO: check certainty of times.
    %       Low variance.  Not near night edge.
  end
end

function  rectangles = check_missing_days (rectangles, cv, pwr, issolar, valid_days)
  % find runs of days with no rectangles
  % for each run, call find_rectangles ();
  % If nothing found and run is long,
  %    recursively halve width and try again
  on_off = [rectangles.on_off];
  on_days = sort (on_off(1,:));
  off_days = sort (on_off(3,:));

  % Days with no rectangles
  none = (on_days(2:end) > off_days(1:end-1));
  st = off_days (none);
  en = on_days ([false, none]);

  new_rectangles = [];
  for i = 1:length (st)
    if en(i) - st(i) > 2
      new_r = find_rect (cv, st(i), en(i), pwr, 1, rectangles);
      % TODO:
      % If new_r is empty, try sub-intervals of st(i):en(i)
    else
      % For short gaps, try to extend a neighbouring rectangle.
      new_r = [];
      ends = find (on_off(3,:) == st(i));
      starts = find (on_off(1,:) == en(i));
      done = false;
      if length (starts) > length (ends)
        for j = 1:length (ends)
          if all (all (cv(rectangles(ends(j)).burst, st(i):en(i)-1) ...
                        > rectangles(ends(j)).power))
            done = true;
            rectangles(ends(j)).on_off(3) = en(i);
            break
          end
        end
      end
      if ~done
        for j = 1:length (starts)
          if all (all (cv(rectangles(starts(j)).burst, st(i):en(i)-1) ...
                        > rectangles(starts(j)).power))
            done = true;
            rectangles(starts(j)).on_off(1) = st(i);
            break
          end
        end
      end
    end
    new_rectangles = [new_rectangles, new_r];
  end
  rectangles = [rectangles, new_rectangles];

  % Remove rectangles
  % If (m)any rectangles are in a region that has low power (early morning)
  %   check all days that have high power in that interval.
end

function [rectangles, kept] = tweak_on_off_pwr (rectangles, cv, pwr, issolar, valid_days);
  % Find "best fit" on time and off time for each rectangle
  % Usually, this will be near on_off([2,4]),
  % but sometimes it may be quite a distance away.
  % See assign_trust / get_weighted_edge
  on_off = [rectangles.on_off];
  cvw = [cv; cv(:,2:end), cv(:,end)];        % wrap
  wd = zeros (size(cvw, 1), ceil (size(cvw, 2) * 5/7), 7);
  for i = 1:7
    idx = true (1, size(cv,2));
    idx ([i:7:size(cv, 2), mod1(i+1,7):7:size(cv,2)]) = false;
    wd(:,1:sum(idx),i) = cvw(:, idx);
  end

  durations = on_off(4,:) - on_off(2,:);
  durations(durations < 0) = durations(durations < 0) + size(cv, 1);
  med_duration = median (durations);
  lambda = mean (abs (log (durations / med_duration)));
  if length (durations) == 1
    lambda = 1;
  end

  for i = 1:length (rectangles)
    % Check locations
    range = ceil (on_off(1,i) * 5/7):floor ((on_off(3,i) - 1) * 5/7);
    quantiles = [0.05, 0.1, 0.15, 0.2 0.25];
    y = zeros (size(cvw, 1), length (quantiles));
    for n = 1:7
      y = max (y, quantile (wd(:, range, n), quantiles, 2));
    end
    x = mean (y, 2);
    z = mean (cv(:, on_off(1,i):on_off(3,i)-1)');
    zz = min (cv(:, on_off(1,i):on_off(3,i)-1)');
    % Hack: look for jumps with 1.5 (now 1.1) times larger power, because
    % initial estimates are typicallyl too low.
    p = rectangles(i).power * 1.1;
%    neighbours = (on_off(1,:) == on_off(3,i) | on_off(3,:) == on_off(1,i));
    neighbours = (on_off(1,:) > on_off(3,i)) | (on_off(3,:) < on_off(1,i));

    overlaps = (on_off(3,:) > on_off(1,i) & on_off(1,:) < on_off(3,i));
    overlaps(i) = false;
    full_overlaps = overlaps & (on_off(3,:) >= on_off(3,:) ...
                                & on_off(1,:) <= on_off(1,:));
    days = get_days (rectangles(i), on_off(:,i), false);

    % Find new on time
    turn_off(1) = on_off(4, i);
    t_off= mod1 (round (turn_off(1)),size(cv, 1));
    turn_on(2) = match_jump (x, -t_off, -p);
    if sum (y(:,1) > 0) > size (cv,1) / 2
      turn_on(3) = match_jump (y(:,1), -t_off, -p);
    else
      turn_on(3)  = 0;
    end
    turn_on(4) = match_jump (x, t_off, -p);
    turn_on(5) = match_jump (z, t_off, -p);
    turn_on(6)= match_jump (zz, t_off, -p);
    turn_on  = [turn_on,  on_off(2, neighbours)];
 
    turn_on = unique_ish (turn_on, size(cv, 1));

    [t, badness] = eval_edges (turn_on, days, on_off(4,i), ...
                               [on_off(1,i), on_off(3,i)], ...
                               pwr, rectangles(i).power, cvw, ...
                               on_off(2, full_overlaps)', true, issolar);
    sm = size_match (turn_on, on_off(4,i), size(cv, 1), med_duration, lambda);
    if max (sm) == 0    % Avoid Inf, since it hides all info in badness
      sm(:) = 1;
    end
    [rectangles(i).badness(1), pos] = min (badness ./ sm);
    new_on = t(pos);

    % Find new off time
    t_on = mod1 (round (new_on), size(cv, 1));
    turn_off(2) = match_jump (x, -t_on,  p);
    if sum (y(:,1) > 0) > size (cv,1) / 2
      turn_off(3) = match_jump (y(:,1), -t_on,  p);
    else
      turn_off(3) = 0;
    end
    turn_off(4) = match_jump (x, t_on,  p);
    turn_off(5) = match_jump (z, t_on,  p);
    turn_off(6) = match_jump (zz, t_on,  p);
    turn_off = [turn_off, on_off(4, neighbours)];

    turn_off = unique_ish (turn_off, size(cv, 1));

    [t, badness] = eval_edges (turn_off, days, new_on, ...
                               [on_off(1,i), on_off(3,i)], ...
                               pwr, rectangles(i).power, cvw, ...
                               on_off(4, full_overlaps)', false, issolar);
    sm = size_match (new_on, turn_off, size(cv, 1), med_duration, lambda);
    if max (sm) == 0    % Avoid Inf, since it hides all info in badness
      sm(:) = 1;
    end
    [rectangles(i).badness(2,1), pos] = min (badness ./ sm);

    new_off = t(pos);
    
    % Redo on time if off time changed a lot.
    % TODO: Avoid copy-and-paste from above.
    if abs (new_off - on_off(4, i)) > 5
      turn_off(1) = new_off;
      t_off= mod1 (round (turn_off(1)),size(cv, 1));
      turn_on(2) = match_jump (x, -t_off, -p);
      if sum (y(:,1) > 0) > size (cv,1) / 2
        turn_on(3) = match_jump (y(:,1), -t_off, -p);
      else
        turn_on(3)  = 0;
      end
      turn_on(4) = match_jump (x, t_off, -p);
      turn_on(5) = match_jump (z, t_off, -p);
      turn_on(6)= match_jump (zz, t_off, -p);
      turn_on  = [turn_on,  on_off(2, neighbours)];

      turn_on = unique_ish (turn_on, size(cv, 1));
      [t, badness] = eval_edges (turn_on, days, on_off(4,i), ...
                                 [on_off(1,i), on_off(3,i)], ...
                                 pwr, rectangles(i).power, cvw, ...
                                 on_off(2, full_overlaps)', true, issolar);
      sm = size_match (turn_on, on_off(4,i), size(cv, 1), med_duration, lambda);
      if max (sm) == 0    % Avoid Inf, since it hides all info in badness
        sm(:) = 1;
      end
      [rectangles(i).badness(1), pos] = min (badness ./ sm);
      new_on = t(pos);
    end

    % Reject or accept new settings.
    if (any (mod (abs ([new_on-on_off(2, overlaps), ...
                       new_off-on_off(4, overlaps)]), size(cv, 1)) < 1) ...
	        && (abs (new_on  - on_off(2, i)) > 5  ...
           || abs (new_off - on_off(4, i)) > 5 )) ...
       || (new_on > new_off && new_on < new_off - 0.5)
      rectangles(i).on_off(2) = -1;
      continue;
    end
    rectangles(i).on_off(2) = new_on;
    rectangles(i).on_off(4) = new_off;
    on_off([2,4],i) = [new_on, new_off]';
    if new_off > new_on
      rectangles(i).burst = ceil (new_on):floor (new_off);
    else
      rectangles(i).burst = [ceil(new_on):size(cv,1), 1:floor(new_off)];
    end
  end

  % Remove duplicate rectangles
  on_off = [rectangles.on_off];
  kept = find (on_off(2,:) > 0);
  rectangles = rectangles (kept);
  if ~any (kept)
    return;
  end
  on_off = [rectangles.on_off];
  overlap_all = bsxfun (@eq, on_off(1,:), on_off(1,:)') ...
              & bsxfun (@eq, on_off(3,:), on_off(3,:)') ...
              & abs (bsxfun (@minus, on_off(2,:), on_off(2,:)')) < 0.5 ...
              & abs (bsxfun (@minus, on_off(4,:), on_off(4,:)')) < 0.5;
      % remove one of mutually overlapping pair
  mutual = triu (overlap_all, 1);
  overlap_r = any (mutual);
  rectangles = rectangles(~overlap_r);
  kept = kept(~overlap_r);

  % Update powers.  ~Min-link clustering on ranksum.
  jumps = {};
  weights = {};
  classes = zeros (1, length (rectangles));
  on_off = [rectangles.on_off];
  [pre_up, ~, post_up] = edge_neighbours (on_off(2,:), size (cv, 1));
  [pre_dn, ~, post_dn] = edge_neighbours (on_off(4,:), size (cv, 1));
  for i = 1:length (rectangles)
    up = diff (cv([pre_up(i), post_up(i)], on_off(1,i):on_off(3,i)-1));
    dn = diff (cv([post_dn(i), pre_dn(i)], on_off(1,i):on_off(3,i)-1));
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

  for i = 1:length (jumps)
    jumps{i} = mean (jumps{i});
  end

  for i = 1:length (rectangles)
    st = rectangles(i).on_off(1);
    en = rectangles(i).on_off(3);
    slice = cv(rectangles(i).burst(1:end-1),st:en-1);
    missed = st - 1 + find (min(slice,[],1) < 0.8 * jumps{classes(i)});
    if length (missed) < 0.25 * (en - st);
      rectangles(i).power = jumps{classes(i)};
    end
  end
end

function idx = unique_ish (vec, hrs)
  % Remove near-duplicate values from  vec.
  % Replace retained values by the mean of all the near-duplicates they
  % represent.
  vec = mod1 (vec, hrs);
  idx = (vec >= hrs + 0.5);
  vec(idx) = vec(idx) - hrs;

  idx = zeros (1, hrs);
  idx (round (vec)) = true;
  idx = find (idx);
  for j = 1:length (idx)
    idx(j) = mean (vec (abs (vec - idx(j)) <= 0.5));
  end
end

function sm = size_match (turn_on, turn_off, slots, med_duration, lambda)
  durations = turn_off - turn_on;
  durations(durations < 0) = durations(durations < 0) + slots;
  sm = exp (-abs (log (durations / med_duration)) / lambda);
end

function [score, ds] = not_ramp (time, on_off, days, cv)
 % Score from 1 if cv(time, days) looks like a two-part step,
 % to 0 if it looks like a linear change.
 % time is the nominal time of the step
 % on_off is 1 for turning on, -1 for turning off
 % days is the range of days that the step/ramp occurs
 % cv is either cv or cvw.

 % Turn time into a string of four half-hours,
 % starting from the "off" state and ending one step after the end of the
 % transition.
 time = floor (time * on_off) - 1;
 time = mod1 ((time:time+3) * on_off, size (cv, 1));
 low_to_high = cv(time, days);

 means = mean (low_to_high, 2);
 vars  = var  (low_to_high, 0, 2);

 mean_last_two = (vars(end-1) * means(end-1) + vars(end) * means(end)) ...
                 / (vars(end-1) + vars(end));
 badness1 = sum ((means(end-1:end) - mean_last_two) .^ 2 ./ vars(end-1:end)) / 2;

 % Least square linear regression
 sum_a   = sum (1 ./ vars);
 sum_ia  = sum ((1:length (vars)) ./ vars');
 sum_i2a = sum ((1:length (vars)).^2 ./ vars');
 sum_y   = sum (means ./vars);
 sum_iy  = sum ((1:length (vars))' .* means ./ vars);
  coeffs = [sum_a, sum_ia; sum_ia, sum_i2a] \ [sum_y; sum_iy];
 linear_fit = coeffs(1) + (1:length (vars))' * coeffs(2);

 badness2 = sum ((means - linear_fit) .^ 2 ./ vars) / 4;

 score = badness2 / (badness1 + badness2);
 
 ds = [1/badness1; 1; 1/badness2];
 ds = ds / sum (ds);
end

function [times, badness] = eval_edges (edges, days, match, st_en, all_pwr, my_pwr, cvw, overlaps, on, issolar)
  times = edges;
  non_missed = times;
  size_match = times;
  miss_slots = times;
  clear_step = times;
  if on
    on_off = 1;
  else
    on_off = -1;
  end

  for n = 1:length (edges)
    [times(n), wgts, ~, jumps] = get_weighted_edge (edges(n), cvw(:, days), ...
                                      all_pwr,  1, 0.2, on);
    wgts = wgts - 0.5;
    % How uniform are the trusted jumps?  Weighted variance.
    [wmean, wvar(n)] = hist_edge (jumps, wgts);
    %wmean = (wgts * jumps') / sum (wgts);
    %wvar(n) = sum ((wgts.^2 .* (jumps - wmean)) .^ 2) / sum (wgts.^4) / (wmean .^2);
    % xmean = diff (x(mod1 ([floor(t(n)-1), ceil(t(n))], size(cv,1))));
    % wmean = 0.5 * (wmean + xmean);
    size_match(n) = wabs (my_pwr - wmean * on_off, -0.3);
    if on
      if times(n) < match
        burst = ceil (times(n)):floor (match);
      else
        burst = [ceil(times(n)):size(cvw,1)/2, 1:floor(match)];
      end
    else
      if times(n) > match
        burst = ceil (match):floor (times(n));
      else
        burst = [ceil(match):size(cvw,1)/2, 1:floor(times(n))];
      end
    end

    if isempty (burst)
      burst = ceil(times(n));
    end
    enough = (cvw(burst, days) > min (my_pwr, all_pwr));
    no_miss = all (enough);
    first_hit = find (no_miss, 1);
    if isempty (first_hit) || first_hit > 0.25 * length (jumps)
      first_hit = 0.25 * length (jumps);
    end
    last_hit = find (no_miss, 1, 'last');
    if isempty (last_hit) || last_hit < 0.75 * length (jumps) + first_hit - 1
      last_hit = 0.75 * length (jumps) + first_hit - 1;
    end
    non_missed(n) = sum (no_miss) ...
                    + 0.5 * (first_hit - 1 + length (jumps) - last_hit);

    % These rolling mins should only be needed for solar.
    enough = rolling_min (-enough, 3);
    enough = -rolling_min (enough', 3)';
    miss_slots(n) = sum (~enough(:)) + ceil (length (enough(:))/100);

    clear_step(n) = not_ramp (times(n), on_off, days, cvw);
  end
  if isempty (overlaps)
    ov = 0;
  else
    ov = any (abs (bsxfun (@minus, overlaps, times)) < 1, 1);
  end
  % TODO: only apply this to rectangles crossing midday
  %       Is it better to nudge up powers instead?
  if issolar
    miss_slots = (miss_slots) .^ 0.75;
  end
  
  % TODO: Should eliminate negative badness, not just take max with 0
  badness = max (0, wvar .* (size_match / all_pwr + 0.25) .* (0.4 + ov) ...
                    .* miss_slots .* (1.1 - clear_step)...
                    ./ (0.5 + non_missed));
end

function [j, v] = hist_edge (jumps, wgts)
  % Detect quality / size of jump edge based on the histogram
  % Strip first and last quartiles
  %
  len = length (jumps);
  if len < 6
    j = mean (jumps);
    v = var (jumps);
    return
  end
  [h, edges] = histcounts (jumps, floor (len/2));
  cs = cumsum (h);
  Q1 = find (cs > len * 0.25, 1);
  Q3 = find (cs < len * 0.75, 1, 'last');
  [biggest, mode] = max (h);
  if isempty (Q3) ||  mode > Q3
    Q3 = mode;
  end
  Q3 = min (Q3+1, len);
  if isempty (Q1) || mode < Q1
    Q1 = mode;
  end

  idx = (jumps >= edges(Q1) & jumps <= edges(Q3));
  w = wgts(idx);
  d = jumps(idx);
  j = (w * d') / sum (w);
  v = sum ((w.^2 .* (d - j)) .^ 2) / sum (w.^4) / (j.^2);

  em = edges(mode);
  em1 = edges (min (mode + 1, length (edges)));
  [next_biggest, pos] = max (h([1:mode-1, mode+1:end]));
  % FIXME: This is sensitive to outliers.
  if em <= 0 && 0 <= em1
    j = 0;
  elseif j < em || j > em1
    if biggest > 2 * next_biggest
      j = edges(mode);
    end
  else
    % if next biggest part of the same peak, compare the peak to the third
    if abs (mode - pos) == 1
      h(mode) = 0;
      h(pos) = 0;
      next_biggest = max (h);
    end
    v = v * next_biggest / biggest;
  end

end

function rectangles = split_rectangle_days (rectangles, cv, orig_runs, issolar, valid_days)
  %show_runs (orig_runs, cv, valid_days);
  %figure(1); show_rectangles (rectangles, cv, []);
end

function rectangles = SNR_trust (rectangles, cv, cor, overlap)
  if isempty (rectangles)
    return
  end
  if nargin < 3
    cor = set_cor_pool (rectangles, cv', 1:size(cv, 2));
  end

  on_off = [rectangles.on_off];
  [pre_t, mid_t, post_t] = edge_neighbours (on_off(2,:), size (cv, 1));
  [pre_b, mid_b, post_b] = edge_neighbours (on_off(4,:), size (cv, 1));

  margin_h = 7;
  left_extent  = max (on_off(1,:) - margin_h, 1);
  right_extent = min (on_off(3,:) + margin_h, size (cv, 2));

  margin_v = 3;
  tops = max (pre_t - margin_v, 1);
  bots = min (post_b + margin_v, size(cv, 1));
  todo = pre_t < post_t & post_t <= pre_b & pre_b < post_b ...
         & tops > 0 & bots <= size (cv,1);
  wrap = (on_off(2,:) > on_off(4,:));

  for r = 1:length (rectangles)
    rectangles(r).ks_p = [NaN; NaN];
  end

  mask = (cor' == 0);

z = zeros (size(cv));
figure(11); show_rectangles (rectangles, z);

  for r = find (todo)
    outer_mask = mask (tops(r):bots(r), left_extent(r):right_extent(r));
    outer_cv   = cv   (tops(r):bots(r), left_extent(r):right_extent(r));
    outer = outer_cv(outer_mask);
    outer = outer(isfinite (outer));
    if isempty (outer)
      continue;
    end

%     z = zeros (size(cv));
%     z (tops(r):bots(r), left_extent(r):right_extent(r)) = outer_mask;
%     figure(11); show_rectangles (rectangles, z);

    %{
    a = (right_extent(r) > on_off(1,:) & left_extent(r) < on_off(3,:) ...
         & ((~wrap & bots(r) > on_off(2,:) & tops(r) < on_off(4,:)) ...
            | (wrap & (bots(r) > on_off(2,:) | tops(r) < on_off(4,:)))));
    b = 0;
    if sum (a) > 1
      %figure(11); show_on_off (on_off(:,r), cv);
      % Allow one other to be nearby
      a(r) = false;
      oth = find (a);
      oo(3:4) = min (on_doff(3:4, oth), [], 2);
      oo(1:2) = max (on_off(1:2, oth), [], 2);
      b(4) = right_extent(r) - oo(1);
      b(3) = oo(3) - left_extent(r);
      b(2) = bots(r) - oo(2);
      b(1) = oo(4) - tops(r);
      b(1:2) = b(1:2) * (b(3) + b(4));
      b(3:4) = b(3:4) * (b(1) + b(2));
      while (right_extent(r) > oo(1) && left_extent(r) < oo(3) ...
           && ((~wrap(r) && bots(r) > oo(2) && tops(r) < oo(4)) ...
                || (wrap(r) && (bots(r) > oo(2) || tops(r) < oo(4)))));
        if all (isinf (b))
          break;
        end
        [~, pos] = min (b);
        switch pos
          case 1
            tops(r) = ceil (oo(4));
            if oo(4) - tops(r) <= margin_v
              pret_t(r) = tops(r) - 1;
            end
          case 2
            bots(r) = floor (oo(2));
            if bots(r) - oo(2) <= margin_v
              bost_b(r) = bots(r) + 1;
            end
          case 3
            %if oo(3) - left_extent(r) <= margin_h
              left_extent(r) = oo(3);
            %end
          case 4
            %if right_extent(r) - oo(1) <= margin_h
              right_extent(r) = oo(1);
            %end
        end
        b(pos) = Inf;
      end
    end

    if all (isinf (b))
      continue;
    end

    %figure(11); show_on_off (on_off(:,r), cv);

    outer1 = cv([tops(r):pre_t(r),post_b(r):bots(r)], ...
                left_extent(r):right_extent(r));
    outer2 = cv(post_t(r):pre_b(r), ...
                [left_extent(r):min(on_off(1,r), right_extent(r))-1, ...
                 max(on_off(3,r), left_extent(r)):right_extent(r)]);
    outer = [outer1(:); outer2(:)];
%}
    inner = cv (round(on_off(2,r)):round (on_off(4,r))-1, on_off(1,r):on_off(3,r)-1);

    % test if these may have come from the same distribution
    % by the Kolmogrov-Smirniv two-sample test.
    [~, a] = kstest2 (outer(:), inner(:));
    %[~, b] = kstest2 (outer, inner(:)-rectangles(r).power);
    b = 0;
    rectangles(r).ks_p = [a; b];

    disp (rectangles(r).ks_p);
    1;
  end

  for r = 1:length (rectangles)
    if length (rectangles(r).burst) <= 2
    else
      days = on_off(1,r):on_off(3,r)-1;
      if isempty (days)
        continue;
      end

      % TODO: find better way of selecting how many to take from idx
      % TODO: consider mid_t.
      top_jump = cv(post_t, days) - cv(pre_t, days);
      [~, idx] = sort (cv(pre_t, days));
      top = top_jump(idx(1:floor (length (idx) / 2 + 1)));

      bottom_jump = cv(pre_b, days) - cv(post_b, days);
      [~, idx] = sort (cv(post_t, days));
      bottom = bottom_jump(idx(1:floor (length (idx) / 2 + 1)));

      left_margin = cor(on_off(1,r):-1:left_extent(r), rectangles(r).burst);
      left_margin = min (left_margin, [], 2);
      left_margin = [find(left_margin, 1)-1, length(left_margin)];
      left_margin = on_off(1,r) - left_margin;

    end

    % Match top
    % Match bottom
    % Match left
    % Match right

    if nargin == 4
      % Match overlap
    end
  end
end

function [ol, overlap_hour, overlap_fraction] = find_overlap (on_off, rectangles)
  overlap_day = min (bsxfun (@minus, on_off(3,:), on_off(1,:)'), ...
                     bsxfun (@minus, on_off(3,:)', on_off(1,:)));
  overlap_day = max (overlap_day, 0);

  overlap_hour = min (bsxfun (@minus, on_off(4,:), on_off(2,:)'), ...
                      bsxfun (@minus, on_off(4,:)', on_off(2,:)));
  cross_midnight = (on_off(4,:) < on_off(2,:));
  overlap_hour(cross_midnight,cross_midnight) = 1;
  % Handle rectangles wrapped over midnight intersecting those that don't
  if any (cross_midnight)
    ov = max (0, bsxfun (@minus, on_off(4,~cross_midnight), on_off(2, cross_midnight)')) ...
       + max (0,-bsxfun (@minus, on_off(2,~cross_midnight), on_off(4, cross_midnight)'));
    overlap_hour(~cross_midnight, cross_midnight) = ov';
    overlap_hour(cross_midnight, ~cross_midnight) = ov;
  end
  overlap_hour = max (overlap_hour, 0);

  overlap = triu (overlap_day .* overlap_hour, 1);
  if any (overlap(:))
    [ol(1,:), ol(2,:)] = find (overlap);
  else
    ol = zeros(2,0);
  end

  if nargout >= 3
    sizes = bsxfun (@max, on_off(3,:) - on_off(1,:), ...
                    cellfun (@length, {rectangles.burst})');

    overlap_fraction = overlap ./ sizes;
  end
end

function [on_off, rectangles]  = align_pair_edge (rectangles, cv, on_off, detail, detail_idx, transition)
  ool = size (on_off, 2);
  if ool <= 1
    return;
  end
  d1 = [0; diff(cv(:))];
  d1(1) = d1(1 + size (cv, 1));       % guess first jump
  d1  = reshape (d1, size (cv));      % half-hour differences

  days = size (cv, 2);
  pwr = mean ([rectangles.power]);

  [ol, overlap_hour, overlap_fraction] = find_overlap (on_off, rectangles);
  ol_frac = overlap_fraction (sub2ind (size (overlap_fraction), ol(1,:), ol(2,:)));
  [~, ol_idx] = sort (ol_frac);

  after = size (cv,2) + 1;
  on_off(:, size (on_off, 2)+1) = [after; 1; after; 1];    % dummy when rectangle is at an edge
  transition(after) = 2;

  for i = ol_idx(end:-1:1)
    % If one is deleted, or shifted so as not to overlap, skip this iteration
    if any (on_off (2, ol(:,i)) < 0) ...
        || all (diff ([on_off([2,4], ol(1,i)); on_off([2,4], ol(2,i))]) > 0) ...
        || all (diff ([on_off([2,4], ol(2,i)); on_off([2,4], ol(1,i))]) > 0)
      continue
    end
    %figure(11); show_on_off (on_off(:, ol(:,i)), cv);
    if on_off (3, ol(1,i)) ~= on_off (3, ol(2,i))

      %  If big overlap
      d = abs (bsxfun (@minus, on_off([1,3], ol(1,i)), on_off([1,3], ol(2,i))'));

      % TODO: Find better criterion for matching ends vs splitting
     % if min (d(1,1), d(2,2)) < min (d(1,2), d(2,1))
     if d(2,2) < min (d(1,2), d(2,1))

        mask_L = false (size (cv,1), 2);
        mask_L(rectangles (ol(1,i)).burst, 1) = true;
        mask_L(rectangles (ol(2,i)).burst, 2) = true;

        % If rectangles overlap enough in time, align their end dates
        if 3 * sum (mask_L(:,1) & mask_L(:,2)) > sum (mask_L(:))
          % find next few rectangles
          a = find (overlap_hour (ol(1,i),:));
          %a = a (abs (on_off(1, a) - on_off(3, ol(1,i))) < 20);
          idx = (on_off(1,a) > on_off(3, ol(1,i)) - 20) ...
                | (on_off(1,a) > on_off(3, ol(1,i)) - 40 ...
                    & on_off(3,a) > on_off(3, ol(1,i)));
          a = a(idx);
          a = setdiff (a, ol(:,i));
          if any (on_off(1,a) < max (on_off(3, ol(:,i))) + 20)
            a = a(on_off(1,a) < max (on_off(3, ol(:,i))) + 20);
          else
            [~, idx] = min (on_off(1,a));
            a = a(idx);
          end

          mask_R = false (size (cv,1), length (a));
          for j = 1:length (a)
            mask_R(rectangles(a(j)).burst,j) = true;    % find right element
          end

          % create dummy rectangle to the right if none exists
          % or if none is to the right of both matrices
          if isempty (a) || max (on_off(1,a)) < max (on_off(3, ol(:,i)))
            a(end+1) = ool+1;
            on_off(:,end) = [size(cv,2)+2; on_off(2, ol(1,i)); size(cv,2)+2; on_off(4, ol(1,i))];
            mask_R(:, end+1) = mask_L(:,1);
          end

          % check for cases
          % Left meets right (case 2: left moves forward, case 3: right move back)
          % Left - gap - right (case 1)
          steps = Inf (2, length(a), 3);


          allowed_keep = [0, 0];
          for n = 1:2
            j = ol(n,i);
            burst = rectangles(j).burst;
                % "sum" to force [] to 0
            mean_before = sum (mean (mean (cv(burst, max(on_off(3,j)-20,on_off(1,j)):on_off(3,j)-1))));
            % Find which days are "allowed",
            % i.e., do not include times the pump seems off.
            min_end = min (on_off(3, j));
            max_end = min (max ([min_end, on_off(1,a)]), size (cv,2) + 1);
             % These next two may be empty
            body = medfilt1 (min (cv (burst, min_end:max_end-1)), 5, 'omitNaN', 'truncate');
            % TODO denominator should be less for non-solar.
            allowed = find (body < rectangles(j).power / 2.5, 1) + min_end - 1;
            allowed_keep(n) = min ([allowed, size(cv,2)]);

            for k = 1:length (a)
              % case 1: gap
              if on_off(3, j) < on_off(1, a(k))
                if a(k) == length (rectangles) + 1
                  if on_off (3,j) == size (cv,2) + 1  % no "after"
                    % TODO: take mean over more than the last column
                    mean_after = mean (cv(rectangles(j).burst, end));
                  else
                    mean_after = mean (mean (cv(rectangles(j).burst, on_off(3,j):size (cv,2))));
                  end
                else
                  mean_after = mean (mean (cv(rectangles(j).burst, on_off(3,j):on_off(1,a(k))-1)));
                end
                %d = mean_before - mean_after;
                %if isempty (d)
                %  steps(j, k, 1) = 0;
                %else
                steps(n, k, 1) = mean_after - mean_before;

                % Bias towards times with existing transitions.
                if transition(on_off(3,j)) > 2
                  steps(n, k, 1) = steps(n, k, 1) - rectangles(j).power / 2;
                end

                %end
              end

              if on_off(1, a(k)) > allowed
                steps(n, k, 2:3) = inf;
                continue;
              end

              % cases 2 and 3: abut
              left_mask = mask_L (:, n);
              right_mask = mask_R (:, k);
              both = left_mask & right_mask;
              left_mask(both) = false;
              right_mask(both) = false;
              if a(k) > length (rectangles)
                me = mean (cv(both, end));
                ak = on_off(1,a(k)) - 2;
              else
                ak = on_off(1,a(k));
                everywhere = cv(both, on_off(3,j):ak);
                me = mean (everywhere (:));
              end

              % case 2: ???
              if any (right_mask)
                post = cv(right_mask, on_off(3,j):ak);
                steps (n, k, 2) = mean (post(:)) - me;
              elseif any (left_mask)
                steps (n, k, 2) = -me;
              end

              % case 3: ???
              if any (left_mask)
                pre = cv(left_mask, on_off(3,j):ak);
                steps (n, k, 3) = mean (pre(:)) - me;
              else
                steps (n, k, 3) = -me;
              end
            end
          end

          % Find best match to the jump power
          type = 3;
          while type == 3
            [s, m] = min (steps(:));    % TODO: handle inhomogeneous powers
            [n, k, type] = ind2sub (size (steps), m);
            if type == 3
              % TODO Before extending to end, we should check jumps.
              % Here, just check that the extension isn't too long
              if on_off(1, a(k)) - max (on_off(3, ol(:,i))) ...
                  > abs (on_off(3, ol(1,i)) - on_off(3, ol(2,i)))
                steps (:,:,3) = Inf;
              else
                break;
              end
            end
            if isinf (s)
              type = 4;
            end
          end

          % Set end of both rectangles to chosen end day.
          % Don't yet set start of matching rectangle, as it may overlap others.
          if type <= 2 && s < 0
            if  1.5 * length (intersect (rectangles(ol(1,i)).burst, ...
                                 rectangles(ol(2,i)).burst))...
                     > min (length (rectangles(ol(1,i)).burst), ...
                            length (rectangles(ol(2,i)).burst))
              en = find_end (on_off(3, ol(:,i)), burst, cv, n, ...
                             squeeze (detail(detail_idx(ol(1,i)), :, :)), ...
                             transition);
              % Only trucate longer if any "strong" edge is
              % already covered by another rectangle
              [mx, pos] = max (on_off(3, ol(:,i)));
              p = ol(pos, i);
              olap_d = (on_off(1,:) < mx & on_off(3,:) > en);
              ov1 = (abs (on_off(2, olap_d) - on_off(2, p)) < 0.5);
              ov2 = (abs (on_off(4, olap_d) - on_off(4, p)) < 0.5);
              % TODO: if transition(en) > 2, look for rectangle in
              %       truncated range
              if mx - en < 5 ...
                 || transition(en) > 2 ...
                 || ((mean (detail (detail_idx(p), en:mx-1, 1)) < 0.3 || sum (ov1) > 1)...
                   &&(mean (detail (detail_idx(p), en:mx-1, 2)) < 0.3 || sum (ov2) > 1))
                on_off(3, ol(pos,i)) = en;
              end
              on_off(3, ol(3-pos, i)) = en;
            end
          elseif type < 4
              % Make rectangle abut the following one
            oo = min (on_off(1, [a(k), ol(1,i), ol(2,i)]), days + 1);
            new_end = min (oo(1), allowed_keep+1);
            clear1 = 0;   % should we delete rectangle 1?
            clear2 = 0;
            
            % If one has a smaller "allowed keep"
            % (i.e., wider range, getting below required power)
            % TODO: look for clock drift.  Handle case of uneven starts.
            [ne,  first] = min (new_end);
            if abs (new_end(1) - new_end(2)) > 10 ...
              && on_off(1, ol(3-first, i)) > on_off(1, ol(first, i)) - 10
                on_off(1, ol(3 - first, i)) = new_end(first);
                ne = new_end;
                clear1 = -length (rectangles);  % don't allow it to become >0
                clear2 = clear1;
            end
            on_off(3,ol(:,i)) = ne;
              % push back the start of "after" rectangles we would overap
            for j = find (on_off(1,a) <= oo(1))
              if isequal (on_off([3,2,4], ol(1,i)), on_off([1,2,4], a(j)))
                if a(j) == length (rectangles) + 1
                  rectangles(ol(1,i)).on_off(3) = days + 1;
                  on_off(3,ol(1,i)) = days + 1;
                  clear1 = -length (rectangles);  % don't allow it to become >0
                elseif oo(2) < rectangles (a(j)).on_off(3)
                  on_off(1,a(j)) = oo(2);
                  clear1 = clear1 + 1;
                end
              elseif isequal (on_off([3,2,4], ol(2,i)), on_off ([1,2,4], a(j)))
                if a(j) == length (rectangles) + 1
                  on_off(1,ol(2,i)) = days + 1;
                  clear2 = -length (rectangles);
                elseif oo(3) < rectangles (a(j)).on_off(3)
                  on_off(1,a(j)) = oo(3);
                  clear2 = clear2 + 1;
                end
              elseif 1.5 * length (intersect (rectangles(a(j)).burst, ...
                                 rectangles(ol(1,i)).burst))...
                     > min (length (rectangles(a(j)).burst), ...
                            length (rectangles(ol(1,i)).burst)) ...
                     && oo(1) < on_off(3, a(j))
                on_off(1,a(j)) = oo(1);
              end
            end
            if clear1 > 0
              on_off(2,ol(1,i)) = -1;
            end
            if clear2 > 0
              on_off(2,ol(2,i)) = -1;
            end
          elseif all (on_off (2, ol(:,i)) >= 0)   % if neither deleted
            % Next rectangle tells us nothing.
            % See if both can be set to the same
            oo1 = rectangles(ol(1,i)).on_off; % can't use on_off, as it is rounded
            oo2 = rectangles(ol(2,i)).on_off;
            if all (abs (oo1([2,4]) - oo2([2,4])) <= 1)
              en = find_end (on_off(3, ol(:,i)), burst, cv, 0, ...
                             squeeze (detail(detail_idx(ol(1,i)), :, :)), transition);
              [mx, pos] = max (on_off(3,ol(:,i)));
              if mx > en + 5
                % If the mismatch looks like a separate rectangle,
                % make it one.
                if abs (diff (on_off(1, ol(:,i)))) < 5
                  j = ol(pos, i);
                  r = find_rect (cv, en, mx, rectangles(j).power, ...
                                 j, rectangles);
                  if ~isempty (r)
                    rectangles (j).on_off = r.on_off;
                    on_off(:, j) = round (r.on_off);
                    on_off(3, ol(3-pos, i)) = en;
%                  st = find_end ([en, mx], burst, cv, -1, ...
%                                  squeeze (detail(detail_idx(ol(1,i)), :, :)));
%                   if st < mx - 2
%                     on_off(1, ol(pos, i)) = st;
%                     on_off(3, ol(3-pos, i)) = en;
                  else
                    on_off(3,ol(:,i)) = en;
                  end
                end
              else
                on_off(3,ol(:,i)) = en;
              end
            end
          end

            % Mark rectangle for removal if they become equivalent
          if all (abs (on_off(:, ol(2,i)) - on_off(:, ol(1,i))) <= [0;1;0;1])
            on_off([1,3], ol(1,i)) = (on_off([1,3], ol(1,i)) ...
                                    + on_off([1,3], ol(2,i))) / 2;
            on_off(2,ol(2,i)) = -1;
          end
        end
      end

      % Shrink one or both rectangles to remove overlap
      if on_off (3, ol(1,i)) ~= on_off (3, ol(2,i)) ...
         && any (on_off([2,4],ol(1,i)) ~= on_off([2,4],ol(2,i))) ...
         && all (on_off(2,ol(:,i)) ~= -1)
        % Should we shrink horizontally, vertically or both?
        % For now, just shrink horizontally.
        a = ol(1,i);
        b = ol(2,i);
        c = [on_off(1,a), on_off(1,b), on_off(3,a), on_off(3,b)];
        if c(1) > c(2)
          [a, b] = deal (b, a);
          c = [on_off(1,a), on_off(1,b), on_off(3,a), on_off(3,b)];
        end
        d = diff (c);
        p = [rectangles([a,b]).power];
        if all (d >= 0)
          on_off = find_switch_day (on_off, a, b, cv, d1, p, pwr, rectangles);
        elseif on_off(1,a) >= on_off(3,b) || on_off(1,b) >= on_off(3,a)
          continue;
        elseif int_overlap (on_off(2,a), on_off(4,a), on_off(2,b), on_off(4,b))
          % one rectangle starts and ends within the other (in days)
            % trust "bigger" jump
          if p(1) / p(2) > 5 || on_off(1,b) == on_off(3,b)
            on_off(2, b) = -1;
          elseif p(1) / p(2) < 0.2 || on_off(1,a) == on_off(3,a)
            on_off(2, a) = -1;
          else
            if on_off(3,a) - on_off(1,a) > on_off(3,b) - on_off(1, b)
              [a, b] = deal (b, a);
            end
            %figure(8); show_on_off (on_off(:, [a,b]), cv);

            m = quantile (cv (:, on_off(1,a):on_off(3,a)-1), 0.25, 2);
            r = [];
            if ~isempty (ranges (find (m > pwr)))
              new_rect = find_rect (cv, on_off(1,a), on_off(3,a), ...
                                 pwr, a, rectangles);
              if ~isempty (new_rect)
                r(1:2,1) = new_rect.on_off([2,4]);
              end
            end
%{
            % This code should find the on/off times of the rectangle,
            % but over-estimates when the top three quantiles have noise.
            m = quantile (cv (:, on_off(1,a):on_off(3,a)-1), 0.25, 2);
            r = ranges (find (m > pwr));
            if isempty (r)
              continue;
            end
            if r(1) == 1 && r(end) == size(cv,1)
              r(1) = r(1,end);
              r = r(:, 1:end-1);
            end
            r = r(:, abs (r(1,:) - on_off(2,a)) <= 1 ...
                   | abs (r(2,:)+1 - on_off(4,a)) <= 1);
            if size(r, 2) > 1
              r = r(:, abs (r(1,:) - on_off(2,b)) > 1 ...
                     | abs (r(2,:)+1 - on_off(4,b)) > 1);
              if size(r,2) > 1
                r = r(:, 1);
              end
            end
%}
            if ~isempty (r) && max (abs (r - on_off([2,4],b))) > 1
              if r(1) <= r(2)
                rectangles(a).burst = r(1):r(2);
              else
                rectangles(a).burst = [r(1):size(cv,1), 1:r(2)];
              end
              rectangles(a).on_off = [on_off(1,a); r(1); on_off(3,a); mod1(r(2)+1, size(cv,1))];
              on_off(:,a) = round (rectangles(a).on_off);
              % If this rectangle is entirely within the other, delete it.
              if isempty (setdiff (rectangles(a).burst, rectangles(b).burst))
                on_off(2,a) = -1;
              end
            else
              on_off(2,a) = -1;
            end

            %figure(8); show_on_off (on_off(:, [a,b]), cv);
          end
        end
      end
    end
    if diff (on_off(3, ol(:,i))) == 0 ...
           && abs (diff (on_off(1, ol(:,i)))) > 4 ...
           && max (abs (diff (on_off([2,4], ol(:,i)), 2))) < 0.5
            % If the on/off times are equal, shrink day range to remove overlap
      [~, pos] = min (on_off(1, ol(:,i)));
      p = ol(pos,i);
      en = find_end ([floor(sum(on_off(1,ol(:,i)))/2), on_off(3, ol(1,i))], ...
                     rectangles(ol(1,i)).burst, cv, 0, ...
                     squeeze (detail(detail_idx(ol(1,i)), :, :)));
      if en > on_off(3, ol(3-pos,i)) - 4
        ngbr = (-1:1) + on_off(1, ol(3-pos, i));
        ngbr = max (min (ngbr, length (transition)), 1);
        if max (transition (ngbr)) < 3
          on_off(2, ol(3-pos, i)) = -1;
        end
      elseif en < on_off(1, ol(3-pos,i))
        on_off(3, ol(pos, i)) = en;
      end
    end
    %figure(11); show_on_off (on_off(:, ol(:,i)), cv);
    i;
  end
  on_off = on_off(:, 1:end-1);      % Remove dummy.
end

function on_off = find_switch_day (on_off, a, b, cv, d1, p, pwr, rectangles)
  % Reconcile the overlap between rectangles(a) and rectangles(b)
  % by estimating the day that the timer settings change.
  % rectangles(a) must start first
  %
  % TODO: Give preference to multiple simultaneous transitions
  days = size(cv, 2);
  e  = on_off(1,b);   % range over which to search for change-over
  e1 = min (size(cv,2), on_off(3,a));

  if p(1) / p(2) > 5                    % trust "bigger" jump
    best_pos = e1 - e + 1;
  elseif p(1) / p(2) < 0.2
    best_pos = 1;
  else
    % TODO: Handle case where one rectangle takes an edge from each
    %       of two other rectangles.
    mask_1 = false (size (cv,1), 1);    % mask of "before" rectangle
    mask_2 = mask_1;                    % mask of "after"  rectangle
    mask_1(rectangles (a).burst, 1) = true;
    mask_2(rectangles (b).burst, 1) = true;

    % copied from align_rectangles; search for mis_before_penalty.
    mask_diff = mask_2 - mask_1;
    idx_diff  = (mask_diff ~= 0);
        % idx_diff can be empty if =height rectangles separated by gap
    wgt_diff  = max(ceil(sum(idx_diff)/2), 1);
    int = max (1,min (days, (e:e1)-1)); % e:e1, clipped %%%%%%%%%%%

    %-% Include partial half-hours
    miss = min(cv(mask_1, int), [], 1) < 0.9*rectangles(a).power;
    if isempty (miss)
      miss = zeros (1, length (int));
    end
    h = floor (rectangles(a).on_off(2));
    frac = 1 - (rectangles(a).on_off(2) - h);
    if h == 0
      h = size(cv, 1);
    end
    miss = miss | cv(h, int) < 0.9 * rectangles(a).power * frac;
    h = ceil(rectangles(a).on_off(4));
    frac = rectangles(a).on_off(4) - h;
    if h > size(cv, 1)
      h = 1;
    end
    miss = miss | cv(h, int) < 0.9 * rectangles(a).power * frac;
    mis_before_penalty = [0,cumsum(miss)];

    %-% Include partial half-hours
    % TODO: Fix duplication with the above, and with cor_pool_precise
    miss = min(cv(mask_1, int), [], 1) < 0.9*rectangles(b).power;
    if isempty (miss)
      miss = zeros (1, length (int));
    end
    h = floor (rectangles(b).on_off(2));
    frac = 1 - (rectangles(b).on_off(2) - h);
    if h == 0
      h = size(cv, 1);
    end
    miss = miss | cv(h, int) < 0.9 * rectangles(b).power * frac;
    h = ceil(rectangles(b).on_off(4));
    frac = rectangles(b).on_off(4) - h;
    if h > size(cv, 1)
      h = 1;
    end
    miss = miss | cv(h, int) < 0.9 * rectangles(b).power * frac;
    mis_after_penalty = [0,cumsum(miss(end:-1:1))];
    mis_after_penalty = mis_after_penalty(end:-1:1);

    pre  = max (e-1:e1, 1);
    post = min (e:e1+1, days);

    oo = round (on_off([2,4], [a,b]));
    terms_a = abs (d1(oo(1), pre) - pwr) ...
            + abs (d1(oo(2), pre) + pwr);
    terms_b = abs (d1(oo(3), post) - pwr) ...
            + abs (d1(oo(4), post) + pwr);
    cost1_a = cumsum (terms_a);
    cost1_b = cumsum (terms_b(end:-1:1));
    cost1_b = cost1_b (end:-1:1);

    cost1 = cost1_a + cost1_b;
    cost2 = sum (abs(bsxfun (@minus, mask_diff(idx_diff), (cv(idx_diff,post)-cv(idx_diff,pre))/wgt_diff)), 1);
    cost3 = mis_before_penalty + mis_after_penalty;

%{
    costs = zeros (1+e1-e,1);
    for j = 1:1+e1-e
      pre = max(e+j-2, 1);                                %%%%%%%%%%%
      post = min(e+j-1, days);                            %%%%%%%%%%%
                  % horizontal edges
      cost = sum(abs(d1((on_off(2,a)),pre)-pwr) ...
                +abs(d1((on_off(4,a)),pre)+pwr)) ...
           + sum(abs(d1((on_off(2,b)),post)-pwr) ...
                +abs(d1((on_off(4,b)),post)+pwr));
                  % vertical edges
                  % TODO:  was cc.  Should it be??
      cost = cost + (sum(abs(mask_diff(idx_diff) - (cv(idx_diff,post)-cv(idx_diff,pre))))/wgt_diff);
      cost = cost + mis_before_penalty(j)+mis_after_penalty(j+1);
      costs(j) = cost;
    end
%}
    costs = cost1 + cost2 + cost3;
    [~, best_pos] = min (costs);
  end
  e = min (e + best_pos - 1, e1);
  on_off (3,a) = e;
  on_off (1,b) = e;
end

function en = find_end (first_second, burst, cv, n, detail, transition)
% Choose most likely end of rectangle, or start if n<0
% first_second is a pair specifying the range in which to end (any order)
% burst is the set of half-hours that the pump is on
  first  = min (first_second);
  second = max (first_second);

  if n >= 0
    mid = first:second-1;
  else
    mid = second-1:-1:first;
  end
  pmid  = mean (cv(burst, mid));

  mean_fwd = cumsum (pmid)           ./ (1:length (pmid));
  mean_bkd = cumsum (pmid(end:-1:1)) ./ (1:length (pmid));
  mean_bkd = mean_bkd(end:-1:1);

  ee = detail (mid,:);
  ee1 = cumsum (max (ee, [], 2) > 0.5)'  - (1:length (ee)) * 0.5;
  ee2 = cumsum (min (ee, [], 2) > -0.1)' - (1:length (ee)) * 0.5;

  total = (mean_fwd - mean_bkd) * length (pmid) + ee1 + ee2;
  if nargin >= 6
    total = total + max (0, max(total)/4 * (transition(mid)' - 2));
  end
  [mx, pos] = max (total);
  if n > 0
    p = max (first_second(n) - first, 1);
    if mx < 1 + total(p);
      pos = p;
    end
  end
  if n >= 0
    en = first + pos;
  else
    en = second - pos;
  end
end

function diffs = vert_diff (used, cv, row, st, en)
  row = round (row);
  if row == 0
    row = size (used, 1);
  end
  prev_row = row - 1;
  if prev_row == 0
    prev_row = size (used, 1);
  end
  diffs = used (row, st:en-1) ~= used (prev_row, st:en-1);
end

function diffs = hor_diff (used, cv, col, burst)
  prev_col = col - 1;
  if prev_col == 0
    prev_col = size (used, 1);
  end
  diffs = used (col, burst) ~= used (prev_col, burst);
end

function x = wabs (x, w)
 % "Weighted" absolute value.  w==-1 gives abs.
  a = (x < 0);
  x(a) = w * x(a);
end

function [detail, enough_power] = top_and_bottom (rectangles, cv, power)
  noisy = (size (cv, 2) > 120) && (mean (cv(:)) > mean (power));
  changed = false (1, length(rectangles));
  on_off = [rectangles.on_off];

 % Threshold based on:
 %    proximity to original rectangle
 %    compatibility with change times in other rectangles
 %    presence or trustworthiness of intervening rectangles
 %    length (or stand-along trustworthiness?)

      % If no next rectangle, or non-overlapping next rectangle, or
      % "unreliable" next rectangle, check if edges extend
  if nargin < 3 || isempty (power)
    power = [rectangles.power];
  end
  bg = zeros (size (cv));

  %-% to be used later: % new_power = zeros (2, length (rectangles));
  on_off_rounded = round (on_off([2,4],:));
  enough_power = ones (length (rectangles), size (cv,2), 2);

  for i = 1:2
%    a = (abs (on_off(2*i,:) - round (on_off(2*i,:))) > 0.2);
%    pre(~a) = on_off_rounded(i,~a) - 1;
%    mid(~a) = pre(~a);
%    pre(a) = floor (on_off(2*i, a)) - 1;
%    mid(a) = pre(a) + 1;
%    post = mid + 1;
%
%    pre  = mod1 (pre,  size(cv, 1));
%    mid  = mod1 (mid,  size(cv, 1));
%    post = mod1 (post, size(cv, 1));
    [pre, mid, post, a] = edge_neighbours (on_off(2*i,:), size (cv, 1));

    % TODO: Recalculate power based on jumps.
    %       Take err1 from original power, err2 from recalculated power
    %for j = 1:length (rectangles)
    %  new_power(i,j) = mean (cv(post(j),on_off(1,j):on_off(3,j)-1) ...
    %                       - cv(pre (j),on_off(1,j):on_off(3,j)-1));
    %end
    if noisy
      c_post = zeros (length (post), size (cv, 2));
      c_mid = c_post;
      c_pre = c_post;
      for j = 1:length (rectangles)
        if (i == 1)
          idx = find (cv(pre(j),:) < abs (power(j)));
        else
          idx = find (cv(post(j),:) < abs (power(j)));
        end
        if length (idx) < 2
          idx = 1:size(cv, 2);
        end

        % interpolate
        cc = cv([pre(j), mid(j), post(j)], idx);
        cc = (3 * cc + 1 * medfilt1 (cc, 3, [], 2)) / 4;
        c_pre (j,1:idx(1)) = cc(1, 1);
        c_mid (j,1:idx(1)) = cc(2, 1);
        c_post(j,1:idx(1)) = cc(3, 1);

        d = diff(idx);
        for k = 1:length (d)-1
          wghts = (0:d(k)-1) / d(k);
          x = idx(k):idx(k+1) - 1;
          c_pre (j,x) = cc(1, k) * (1-wghts) ...
                      + cc(1, k+1) *  wghts;
          c_mid (j,x) = cc(2, k) * (1-wghts) ...
                      + cc(2, k+1) *  wghts;
          c_post(j,x) = cc(3, k) * (1-wghts) ...
                      + cc(3, k+1) *  wghts;
        end
        c_pre (j,idx(end):end) = cc(1, end);
        c_mid (j,idx(end):end) = cc(2, end);
        c_post(j,idx(end):end) = cc(3, end);
      end
    else
      c_pre  = cv(pre, :);
      c_mid  = cv(mid, :);
      c_post = cv(post,:);
    end

    if isfield (rectangles, 'alt_days')
      has_alt = find (~cellfun (@isempty, { rectangles.alt_days }));
      if ~isempty (has_alt)         % alt_days flag may be out of date
        alt_rect = [rectangles(has_alt)];
        alt_hrs = cellfun(@(x)(x(:)), {alt_rect.alt_hrs}, 'UniformOutput', false);
        alt_hrs = [alt_hrs{:}];
        [alt_pre, alt_mid, alt_post] ...
            = edge_neighbours (alt_hrs(i,:), size (cv, 1));
        for r = 1:length (has_alt)
          days = [rectangles(has_alt(r)).alt_days(1):7:size(cv, 2), ...
                  rectangles(has_alt(r)).alt_days(2):7:size(cv, 2)];
          c_pre (has_alt(r), days) = cv(alt_pre (r), days);
          c_post(has_alt(r), days) = cv(alt_post(r), days);
          % TODO: if a(has_alt(r)), set c_mid based on on_off's fraction(a)
          c_mid (has_alt(r), days) = cv(alt_mid (r), days);
        end
      end
    end

    c_pre  = bsxfun (@rdivide, c_pre,  power(:));
    c_mid  = bsxfun (@rdivide, c_mid,  power(:));
    c_post = bsxfun (@rdivide, c_post, power(:));

    % Give less weight to jumps that are too big,
    % since power is often under-estimated,
    % and consistent big jumps are likely to be interesting.
    err = wabs (1 - medfilt1 (c_post - c_pre, 5, [], 2), -0.25);

    % Is the "on" power high enough?
    % Allow an odd exception, and weekends.
    ep = medfilt1 (max (abs (c_post), abs (c_pre)), 3, [], 2);
    enough_power(:,:,3-i) = logical (medfilt1 (single(ep > 0.7), 7, [], 2));
%    enough_power(:,:,3-i) = (medfilt1 (max (abs (c_post), abs (c_pre)), 7, [], 2) > 0.8);

    fraction = on_off(2*i,:) -  floor (on_off (2*i,:));
    if any (a)
      err_pos(a,:) = wabs (bsxfun (@minus, 1-fraction(a)', medfilt1 (c_mid(a,:) - c_pre(a,:), 5, [], 2)), -0.25);
      err(a,:) = (err(a,:) + err_pos(a,:)) / 2;
    end

    smoothed = medfilt1 (-rolling_min (-rolling_min (err', 9), 9))';
    distrust = 1 ./ (abs (c_post + c_pre) + 0.5);
    detail (:,:,3-i) = (0.5 - smoothed) ./ distrust;

    bg (pre+1,:) = detail (:,:,3-i);
    power = -power;
  end
end

function [rectangles, detail, new_rectangles, new_details] = extend_flat(rectangles, cv, power, issolar, all_on_off)
 % Top/bottom edges:
 %   Find jumps for whole image
 %   for top and bottom edges of each rectangle,
 %      find exp (-error between observed difference and edge size)
 %   Subtract a middle-ish value
 %   Weight by exp (-local variance) -- over what interval?
 %                                   -- average of horiz variances?
 %   (This gives big negative weight to confident mismatches, instead of
 %    giving them no weight.)
 %   Take average of top and bottom edge
 %   Median filter
 %
 % Side edge:
 %   Extract code from align edges
 %   Vectorize it
 %
 % This is currently too sensitive to the estimated power.
  if isempty (rectangles)
    detail = [];
    new_rectangles = [];
    new_details = [];
    return
  end

  on_off = [rectangles.on_off];

  cvw = [cv; cv(:,2:end), cv(:,end)];   % wrap
  cv_rev = cv(:, end:-1:1);             % Reverse, to reuse code for left/right
  cvw_rev = cvw(:, end:-1:1);
  on_off_rev([2,4],:) = on_off([2,4],:);
  on_off_rev([1,3],:) = size(cv, 2) + 2 - on_off([3,1],:);
  d1 = reshape (diff ([cv(end); cv(:)]), size (cv));

  new_rectangles = [];
    if nargin < 3 || isempty (power)
    power = [rectangles.power];
  end
  pwr = mean (power);
  if isscalar (power)
    power = repmat (power, size(rectangles));
  end

  [detail, enough_power] = top_and_bottom (rectangles, cv, power);

  % Candidate points:
  % - When min (detail, [], 3) changes sign
  % - Big jumps in the above
  % - When max (detail, [], 3) changes sign

  % Compare with jumps in (medfilt1'd) "on" value
  % Compare with start/end of other rectangles,
  %      especially those overlapping in time.

  % Generate new "candidate" regions
  % Discard if there are already reliable rectangles there

  mn = min (detail, [], 3) .* (enough_power(:,:, 1) & enough_power(:,:,2));
  mx = max (detail, [], 3) .* (enough_power(:,:, 1) & enough_power(:,:,2));
  mn1 = 0.9*mn + 0.1*mx;
  mx1 = 0.9*mx + 0.1*mn;
  start_missed = 0;
  end_missed = 0;
  updated = 0;
  for i = 1:length (rectangles)
    st = -1;
    if all (mx(i, on_off(1,i):on_off(3,i)-1) > 0)
      st1 = find ([0, mn(i, 1:on_off(1,i))]   <= 0.03, 1, 'last');
      st2 = find ([0, mx(i, 1:on_off(1,i))]   <= 0.03, 1, 'last');
      st3 = find ([0,mn1(i, 1:on_off(1,i))]   <= 0.03, 1, 'last');
      st4 = find ([0,mx1(i, 1:on_off(1,i))]   <= 0.03, 1, 'last');
      en1 = find ([ mn(i, on_off(3,i):end),0] <= 0.03, 1) + on_off(3,i)-1;
      en2 = find ([ mx(i, on_off(3,i):end),0] <= 0.03, 1) + on_off(3,i)-1;
      en3 = find ([mn1(i, on_off(3,i):end),0] <= 0.03, 1) + on_off(3,i)-1;
      en4 = find ([mx1(i, on_off(3,i):end),0] <= 0.03, 1) + on_off(3,i)-1;

      %%%%%
      % If st and/or st2 differs greatly from current start,
      %   Recalculate both using new_power, giving five possibilities?
      %   Look for shared big jump before and after, rather than min/max
      %       of jumps?
      %   Choose to corresponding to the same power level
      %   Update power level to match.
      %   Do we consider whether these jumps are "close" to the original?
      %   Do we consider proximity to other rectangles?

      % Hack -- should deal with discrepancy more carefully
      st = floor (sort ([on_off(1,i), st1, st2, st3, st4]));
      en = floor (sort ([on_off(3,i), en1, en2, en3, en4]));

      % Currenly, just check that if we are extending the range,
      % then we're not extending beyond an obvious "edge"
      if length (en) > 2 && en(ceil (end/2)) >= on_off(3,i) && en(end) - en(1) > 5
        burst = rectangles(i).burst;
        if en(1) <= 2   % avoid illegal index in the next statement.
          en(1) = 3;
        end
        mins = min (cv(burst, en(1)-2:en(end)-1), [], 1);
        mins = medfilt1 ([mins(1), mins(1), mins, mins(end), mins(end)]);
        mins = mins(3:end-2);
        jumps = wabs (mins (en - en(1) + 2) - mins (en - en(1) + 1) - power(i), -0.3);
        % TODO Also get jumps in rows just above and just below rectangle.
        above = mod1 (burst(1) - 1, size (cv, 1));
        below = mod1 (burst(end) + 1, size (cv, 1));
        edge = (cv([above,below], en(1)-2:en(end)-1));
        edge = medfilt1 ([edge(:,1), edge(:,1), edge, edge(:,end), edge(:,end)], 3, [], 2);
        edge = edge(:,3:end-2);
        jumps_edge = wabs (edge (:,en - en(1) + 2) - edge (:,en - en(1) + 1) - power(i), -0.3);
        jumps_edge = min (jumps_edge);

        weights = min (1:length (en), length (en):-1:1);
        weights = 1 + weights(ceil (end/2)) - weights;
        [~, idx1] = min ((0.025 + min (jumps, jumps_edge)) .* weights);
        [~, idx2] = min ((0.1   + min (jumps, jumps_edge)) .* weights);
        if en(idx1) ~= en(idx2)
          s = min (idx1, idx2);
          e = max (idx1, idx2);
          if (mean (mn(en(s):en(e)) + mx(en(s):en(e))) < 0)
            idx1 = s;
          else
            idx1 = e;
          end
        end
        en = en(idx1);

%{
        % Check for "trustworthy" rectangles that the extension would overlap.
        % Reduce the degree of extension to avoid such overlap.
        if en > on_off(3,i) + 5
          overlap_days = (on_off(1,:) <= en & on_off(1,:) >= on_off(3,i));
          overlap_days(i) = false;
          overlap_days = find (overlap_days);

          % Process from closest to furthest
          [~, idx] = sort (on_off(1, overlap_days));
          overlap_days = overlap_days(idx);

          for n = 1:length (overlap_days)
            r = overlap_days(n);
            mid1 = any (rectangles(i).burst == ceil (on_off(2, r)));
            mid2 = any (rectangles(i).burst == floor (on_off(4,r)));
            if mid1 || mid2         % if overlap hours
              % Compare means of four rectangles,
              % delimited by ...
              en1 = min (on_off(3,r), en);
              st1 = max (2 * en - en1, on_off(1,i));
              if mid1
                mid = ceil (on_off(2, r));
              else
                mid = floor (on_off(4,r));
              end
              if on_off(2,i) < on_off(4,i)
                burst1 = ceil (on_off(2,i)):mid-1;
                burst2 = mid:on_off(4,i)-1;
              else
                if on_off(2,i) < mid
                  burst1 = ceil(on_off(2,i)):mid-1;
                  burst2 = [mid:size(cv,1), 1:on_off(4,i)-1];
                else
                  burst1 = [on_off(2,i):size(cv,1), 1:mid-1];
                  burst2 = mid:on_off(4,i)-1;
                end
              end
              % Check edges
              [~, b1] = eval_edges (on_off(2,[i,r]), on_off(4,i), [on_off(1,r), en1], power(i), rectangles(i).power, cvw, [], true);
              [~, b2] = eval_edges (on_off(4,[i,r]), on_off(2,i), [on_off(1,r), en1], power(i), rectangles(i).power, cvw, [], false);
              badness = b1 + b2;
              if badness(2) < badness(1)
                if mid2
                  [burst1, burst2] = deal (burst2, burst1);
                end
                % Check means
                mis_left = cv(burst1, st1:on_off(1,r)-1);
                hit_left = cv(burst2, st1:on_off(1,r)-1);
                mis_rght = cv(burst1, on_off(1,r):en1);
                hit_rght = cv(burst2, on_off(1,r):en1);
                if isempty (burst1) || isempty (burst2) || mean (top_rght(:)) < min ([mean(top_left(:)), mean(bot_left(:)), mean(bot_rght(:))])
                  en = on_off(1,r);
                  break
                end
              end
            end
          end
        end
%}
%{
        burst = rectangles(i).burst;
        region1 = cv (burst, st(end):en(1)-1);
        region2 = cv (burst, en(1):en(2)-1);
        region3 = cv (burst, en(2):en(3)-1);

        region1 = mean (region1(:));
        region2 = mean (region2(:));
        region3 = mean (region3(:));

        if region2 < 0.9 * power(i) ...
           || abs   (power(i) + region2 - region1) < ...
              0.3 * (power(i) + region3 - region2)
          en = en(1);
        else
          en = en(2);
        end
%}
      else
        en = en(ceil (end/2));
      end

      start = on_off(3,i);
        % If we end inside another rectangle, possibly reduce start
      if on_off(2,i) < on_off(4,i)
        overlap = on_off(2,:) < on_off(4, i) & on_off(4,:) > on_off(2,i);
      else
        overlap = on_off(2,:) > on_off(4,:) | on_off(2,:) < on_off(4,i) | on_off(4,:) > on_off(2,i);
      end
      overlap = overlap & on_off(1,:) <= en & on_off(3,:) > en;
      overlap(i) = false;
      if any (overlap)
        start = max (on_off(1,i), min (on_off(1,overlap)));
      end

      r = split_at_overlap ([on_off(1:2,i); en; on_off(4,i)], on_off, cv, ...
                            cvw, start, i, rectangles, true, issolar);
      if length (r) >= 3
        en = r(3,1);
      end

      % If we overlap the next rectangle in our chain, find the best
      % crossover day.

%       if chain_next(i) && on_off(1, chain_next(i)) < en
%         p = [rectangles([i, chain_next(i)]).power];
%         oo = find_switch_day (on_off, i, chain_next(i), cv, d1, p, pwr, rectangles);
%         en = oo(3, i);
%       end
      en = keep_matched_switches (en, on_off, cv, cvw, i, rectangles, issolar);

      % Repeat for starting edge.
      if length (st) > 2 && st(ceil (end/2)) <= on_off(1,i) && st(end) - st(1) > 5
        burst = rectangles(i).burst;
        st(st >= size (cv, 2)) = size (cv, 2) - 1;
        mins = min (cv(burst, st(1):st(end)+1), [], 1);
        mins = medfilt1 ([mins(1), mins(1), mins, mins(end), mins(end)]);
        mins = mins (3:end-2);
        jumps = abs (mins (st - st(1) + 2) - mins (st - st(1) + 1) - power(i));
        % TODO Also get jumps in rows just above and just below rectangle.
        above = mod1 (burst(1) - 1, size (cv, 1));
        below = mod1 (burst(end) + 1, size (cv, 1));
        edge = (cv([above,below], st(1):st(end)+1));
        edge = medfilt1 ([edge(:,1), edge(:,1), edge, edge(:,end), edge(:,end)], 3, [], 2);
        edge = edge(:,3:end-2);
        jumps_edge = wabs (edge (:,st - st(1) + 2) - edge (:,st - st(1) + 1) - power(i), -0.3);
        jumps_edge = min (jumps_edge);

        weights = min (1:length (st), length (st):-1:1);
        weights = 1 + weights(ceil (end/2)) - weights;
        [~, idx1] = min ((0.025 + min (jumps, jumps_edge)) .* weights);
        [~, idx2] = min ((0.1   + min (jumps, jumps_edge)) .* weights);
        if st(idx1) ~= st(idx2)
          s = min (idx1, idx2);
          e = max (idx1, idx2);
          if (mean (mn(st(s):st(e)) + mx(st(s):st(e))) < 0)
            idx1 = s;
          else
            idx1 = e;
          end
        end
        st = st(idx1);


%{
        region1 = cv (burst, st(1):st(2)-1);
        region2 = cv (burst, st(2):st(3)-1);
        region3 = cv (burst, st(3):en(1)-1);

        region1 = mean (region1(:));
        region2 = mean (region2(:));
        region3 = mean (region3(:));

        if region2 < 0.9 * power(i) ...
           || abs   (power(i) + region2 - region3) < ...
              0.3 * (power(i) + region1 - region2)
          st = st(3);
        else
          st = st(2);
        end
%}
      else
        st = st(ceil (end/2));
      end

      % TODO: consider overlap existing before the extension,
      %       as done for en.
      r = split_at_overlap ([on_off_rev(1:2,i); size(cv,2)+2-st; on_off_rev(4,i)], ...
                            on_off_rev, cv_rev, cvw_rev, on_off_rev(3,i), ...
                            i, rectangles, true, issolar);
      if length (r) >= 3
        st = size(cv,2)+2 - r(3,1);
      end
    else      % if mx isn't positive for the whole range
      a = ranges (find (mn (i, on_off(1,i):on_off(3,i)-1) > 0));
      if length (a) == 2 && a(2) - a(1) > 0.6 * (on_off(3,i)-on_off(1,i));
        fprintf ('st %d + %d  en %d + %d\n', on_off(1,i), on_off(1,i)+a(1)-1, on_off(3,i), on_off(1,i)+a(2));
        st = on_off(1,i) + a(1) - 1;
        en = on_off(1,i) + a(2);
      else
        fprintf ('smaller\n');
      end
    end
    if st > 0 && (st ~= on_off(1,i) || en ~= on_off(3,i))
%fprintf ('st %g (was %g) en %g (was %g)\n', st, on_off(1,i), en, on_off(3,i));
      if en < on_off(3,i) - 5
        % Look for "real" rectangle in truncated section
        new_rect_s = find_rect (cv, en, on_off(3,i), ...
                                rectangles(i).power, i, rectangles);
        if ~isempty (new_rect_s)
          % don't count it if it is part of an existing rectangle
          if any (abs (on_off(2,:) - new_rect_s.on_off(2)) < 1 ...
                & abs (on_off(4,:) - new_rect_s.on_off(4)) < 1 ...
                & on_off(1,:) <= new_rect_s.on_off(1) ...
                & on_off(3,:) >= new_rect_s.on_off(3))
            new_rect_s = [];
          else
            new_rect_s.colour = [1; 1; 0];
          end
        end
        new_rectangles = [new_rectangles, new_rect_s];
      end
      on_off([1,3], i) = [st, en];
      on_off_rev([3,1], i) = size(cv, 2) + 2 - [st, en];
      rectangles(i).on_off(1) = st;
      rectangles(i).on_off(3) = en;
      rectangles(i).colour = [0 1 0]';

      slice = cv(rectangles(i).burst(1:end-1),st:en-1);
      rectangles(i).missed = st - 1 + find (min(slice,[],1) < 0.8 * rectangles(i).power);

      changed(i) = true;
      updated = updated + 1;
    end
    %figure(100); show_rectangles (rectangles, cv);
  end
%fprintf ('updated %d  start_missed %d  end_missed %d\n', updated, start_missed, end_missed);fprintf('length (rectangles) = %d\n', length (rectangles));
  if ~isempty (new_rectangles)
    new_details = top_and_bottom (new_rectangles, cv, [new_rectangles.power]);
  else
    new_details = [];
  end
  figure (9); show_rectangles ([rectangles, new_rectangles], cv, []);
  if isfield (rectangles, 'colour')
    rectangles = rmfield (rectangles, 'colour');
  end
  if ~isempty (new_rectangles)
    new_rectangles = rmfield (new_rectangles, 'colour');
  end
end

function en = keep_matched_switches (en, on_off, cv, cvw, i, rectangles, issolar)
  if rectangles(i).pairings > 1
    en = on_off(3, i);
  end
end

function [new_rect adjusted_rect] = split_at_overlap (r_on_off, on_off, cv, cvw, st, i, rectangles, abort, issolar)
  % Split the rectangle described by 4x1 r_on_off
  % at days where it overlaps with rectangles  that seem to fit better,
  % starting from day st.
  % Rectangles to check against are in on_off and rectangles.
  % r_on_off is assumed to have the same power and burst as rectangles(i).
  % If abort is true, stop after the first split
  %   (for use extending rectangles, rather than creating new ones).
  new_rect = r_on_off;
  adjusted_rect = zeros (5,1);

  if r_on_off(3) > st
    en = r_on_off(3);
    overlap_days = (on_off(1,:) <= r_on_off(3) & on_off(3,:) >= st);
    overlap_days(i) = false;
    overlap_days = find (overlap_days);

    % Process from closest to furthest
    [~, idx] = sort (on_off(1, overlap_days));
    overlap_days = overlap_days(idx);

    j = 1;
    for n = 1:length (overlap_days)
      r = overlap_days(n);
      mid1 = any (rectangles(i).burst == ceil  (on_off(2,r)));
      mid2 = any (rectangles(i).burst == floor (on_off(4,r)));
      mid3 = any (rectangles(r).burst == rectangles(i).burst(1)) ...
             || any (rectangles(r).burst == rectangles(i).burst(end));
      if (mid1 || mid2 || mid3) ...                  % if hours overlap but are
          && max (abs (on_off([2,4], r) - r_on_off([2,4]))) > 0.5  % not the same
        % Compare means of four rectangles,
        % delimited by ...
        en1 = min (on_off(3,r), en);
        st1 = max (2*on_off(1,r)-en1, on_off(1,i));
        %st1 = max (on_off(1,r), on_off(1,i));
        if mid1
          mid = ceil (on_off(2, r));
          if mid > size (cv,1)
            mid = 1;
          end
        else
          mid = floor (on_off(4,r));
          if mid < 1
            mid = size (cv, 1);
          end
        end
        if on_off(2,i) < on_off(4,i)
          burst1 = ceil (on_off(2,i)):mid-1;
          burst2 = mid:on_off(4,i)-1;
        else
          if on_off(2,i) < mid
            burst1 = ceil(on_off(2,i)):mid-1;
            burst2 = [mid:size(cv,1), 1:on_off(4,i)-1];
          else
            burst1 = [ceil(on_off(2,i)):size(cv,1), 1:mid-1];
            burst2 = mid:on_off(4,i)-1;
          end
        end

        % skip four-corner test if the overlap rectangle is taller.
        if (~mid1 && ~mid2)
          burst1 = [];
        end

        % Tentatively find change-over point in timer settings
        %   Maybe pass subset of rectangles.
        %   Maybe precompute d1 in caller
        %   If  b  is marked as deleted, what do we do?
        %   If  a  is marked as deleted, what do we do?
        % [on_off, rectangles] = find_switch_day (on_off, a, b, cv, d1, pwr, rectangles)


        % Check edges
        days = get_days (rectangles(r), on_off(:,r), false);
        [~, b1] = eval_edges (on_off(2,[i,r]), days, on_off(4,i), ...
                              [on_off(1,r), en1], rectangles(i).power, ...
                              rectangles(i).power, cvw, [], true, issolar);
        [~, b2] = eval_edges (on_off(4,[i,r]), days, on_off(2,i), ...
                              [on_off(1,r), en1], rectangles(i).power, ...
                              rectangles(i).power, cvw, [], false, issolar);
        badness = b1 + b2;
        if badness(2) < badness(1)
          if mid2
            if ~mid1
              [burst1, burst2] = deal (burst2, burst1);
            else
              % TODO: shrink burst 2
            end
          end
          % Check means
          mis_left = cv(burst1, st1:on_off(1,r)-1);
          hit_left = cv(burst2, st1:on_off(1,r)-1);
          mis_rght = cv(burst1, on_off(1,r):en1-1);
          hit_rght = cv(burst2, on_off(1,r):en1-1);
          if isempty (burst1) || isempty (burst2) ...
             || mean (mis_rght(:)) < min ([mean(mis_left(:)), ...
                                           mean(hit_left(:)), ...
                                           mean(hit_rght(:))])
            if on_off(1,r) > new_rect(1,j) + 2
%}
              new_rect(3,j+1) = new_rect(3,j);
              new_rect(3, j) = on_off(1,r);
              j = j + 1;
            end
            new_rect(1,j) = on_off(3,r);
            if abort || new_rect(3,j) < on_off(3,r)
              break
            end
          end
        end
      end
    end
    if new_rect(1,j) > new_rect(3,j)
      new_rect = new_rect(:,1:end-1);
    end
  end
  % Omit empty rectangles
  new_rect = new_rect(:, new_rect(3,:) > new_rect(1,:));
end

function days = get_days (rect, on_off, alt)
  % Find days that should match an edge of rectangle.
  % If alt is true, use alt_days.
  days = on_off(1):on_off(3)-1;
  if isfield (rect, 'alt_days') && ~isempty (rect.alt_days)
    d = false (1, days(end));
    d(days) = true;
    if all (rect.on_off([1,3]) == on_off([1,3]))
      alt_days = rect.alt_days;
    else
      alt_days = mod1 (366 - rect.alt_days, 7);    % TODO: Don't assume 365.
    end
    skip = [alt_days(1):7:days(end), alt_days(2):7:days(end)];
    if nargin < 3 || ~alt
      d(skip) = false;
    else
      d(~skip) = false;
    end
    days = d;
  end
end

function [new_rect_s] = find_rect (cv, st, en, power, i, rectangles)
  % If most but not all of st:en is a potential rectangle extension,
  % look for a rectangle within the restricted set of days.
  new_rect_s = [];
  len = length (rectangles(i).burst);
  if len >= 5
    r = 3:len - 2;
  elseif len >= 3
    r = 2:len - 1;
  else
    r = 1:len;
  end
  m = min (cv (rectangles(i).burst(r), st:en-1));
  r = ranges (find (m > power));
  [m, pos] = max (diff (r, 1));
  if m >= 0.7 * (en-st)
    en = st + r(2,pos);
    st = st + r(1,pos) - 1;
  end

  % Step 1:
  %  Min over width
  %  Look for jumps of size  power
  %  Look for jumps spaced by on_off(4,i) - on_off(2,i)
  %  Bonus if location is close to on_off(2,i):on_off(4,i)
  %  Penalty if overlaps existing rectangle in on_off, except i.
  m = min (cv (:, st:en-1), [], 2);
  if 0
    [u_jp, times, jp] = find_jumps (m, 0.01, 1);
    times = times + 1;        % TODO: is this an initial offset issue?
  else
    mm = m([end-2:end, 1:end, 1:3]);
    r = ranges (find (mm > power));
    if isempty (r)
      return
    end
    short = r(:, r(1,:)+2 > r(2,:));
    r = r(:, r(1,:)+2 <= r(2,:));   % only consider rectangles of length > 2

    r(r <= 3) = r(r <= 3) + size (cv, 1);
    idx = (r > size (cv,1) + 3);
    r(idx) = r(idx) - size (cv, 1);
    if size(r,2) > 1 && all (r(:, end) == r(:,1))
      r = r(:, 1:end-1);
    end

    jp_up = max ([mm(r(1,:)), mm(r(1,:)+1), mm(r(1,:)+2)]') - min ([mm(r(1,:)-1), mm(r(1,:)-2)]');
    jp_dn = min ([mm(r(2,:)+1), mm(r(2,:)+2)]') - max ([mm(r(2,:)), mm(r(2,:)+1), mm(r(2,:)+2)]');
    jp_up = max (jp_up, 1e-10);
    jp_dn = min (jp_dn, 0);
    times = [r(1,:), r(2,:)] - 3;
    jp = [jp_up, jp_dn];
    u_jp = abs (jp);
  end

  score = power ./ (0.2 * power + wabs (u_jp - power, -0.2));
  up = (jp > 0);
  score(up)  = score(up)  + 10 ./ (10 + abs (times(up)  - rectangles(i).on_off(2)));
  score(~up) = score(~up) + 10 ./ (10 + abs (times(~up) - rectangles(i).on_off(4)));
  pairs = bsxfun (@plus, score(up), score(~up)');
  pairs = pairs + 0.5 ./ (0.5 + bsxfun (@minus, times(~up)', times(up) + rectangles(i).on_off(4)-rectangles(i).on_off(2)));
  pairs = pairs + power ./ (0.2 * power + abs (bsxfun (@minus, u_jp(~up)', u_jp(up))));
  down = find (~up);
  up = find (up);

  oo = [rectangles.on_off];
  olap = [oo(:, oo(3,:) >= en-2 & oo(1,:) <= st-2 & i ~= 1:size(oo,2))];
  for n = 1:length (up)
    for j = 1:length (down)
      burst = ceil (times(up(n))):floor (times(down(j)));
      if isempty (burst)
        burst = [ceil(times(up(n))):size(cv,1), 1:floor(times(down(j)))];
      end
      if any (m(burst) < power)
        pairs(n,j) = 0;
      end
      % Skip if we substantially overlap an existing rectangle
      if any (max (abs ([burst(1)-olap(2,:) burst(end)-olap(4,:)+1])) <= 1)
        pairs(n,j) = 0;
      end
    end
  end

  % For now, just take the best rectangle, if it is good enough
  % Should sometimes take two.
  [~, pos] = max (pairs(:));
  [u, d] = ind2sub (size(pairs), pos);
  if pairs (u, d) > 5
    new_rect_s = rectangles(i);
    on_off = [st; times(up(u)); en; mod1(times(down(d))+1, size(cv,1))];
    new_rect_s.on_off = on_off;
    burst = ceil (on_off(2)):floor (on_off(4));
    if isempty (burst)
      burst = [ceil(on_off(2)):size(cv,1), 1:floor(on_off(2))];
    end
    new_rect_s.burst = burst;
  end
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

function rectangles = extend_irregular (rectangles, cv, power)
end

function [timer, rectangles] = timer_settings(rectangles, cv, issolar)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Convert a list of "rectangles" into a list of dates at which timer settings
 % change and the settings during those intervals.
 % Fields:
 %    start time    - day (or 30 min?) at which settings were changed
 %    rects         - list of indices into 'rectangles' structure after this time
 %    alt_days      - days of the week using non-standard timer settings
 %    alt_rects     - list of indices into 'rectangles' for non-standard settings
 %
 %%%    start time    - day (or 30 min?) at which settings were changed
 %%%    on/off times  - array of times (in fraction of a 30-min interval) for
 %%%                    switching on and off
 %%%    power         - some users seem to have double-power or half-power modes.
 %%%                    record which intervals are on unusual power
 %%%    days of week  - some users have different settings on different days of
 %%%                    the week.  They seem to have a choice of two settings
 %%%  structure of arrays
  on_off = [rectangles.on_off];
  evDays = [1, union(on_off(1,:), on_off(3,:))];
  T = length(evDays);
  power = common_power(rectangles);
  timer.start = zeros(T+1,1);
  timer.power = power * ones(T,1);
  timer.rects = cell(T,1);
  timer.alt_days = zeros(T,0);
  timer.alt_rects = zeros(T,0);
  for e = 1:T
    timer.start(e) = evDays(e);
                        % Comment was: will be <= and > after fixing counting in runs/on_off
    timer.rects{e} = find (on_off(1,:) <= evDays(e) & on_off(3,:) > evDays(e));
  end

  % to find alt_days, consider three possibilities:
  % 1. all days are the same
  % 2. A group of days have different settings throughout the year
  % 3. The set of days having different settings changes from time to time.
  % Heuristics include:
  % - Days with alternate settings are usually contiguous
  % - if some rectangles have alternate settings, many probably do
  % - alternate settings usually don't differ too much from main settings

  [~, idx] = sort(-(on_off(3,:) - on_off(1,:)));        % most days first
  %days = 1:7;
  L = length (idx);
  d_off = cell (L, 1);
  d_on  = cell (L, 1);
  t_off = zeros(L, 1);
  t_on  = zeros(L, 1);
  c_off = zeros(L, 1);
  c_on  = zeros(L, 1);
  if isfield (rectangles, 'alt_days')
    rectangles = rmfield (rectangles, {'alt_days', 'alt_hrs'});
  end
  for i = idx
    % TODO: allow rect of < 2 weeks if it has the same on or off time
    %       as its predecessor or successor.
    if on_off(3,i) - on_off(1,i) >= 14   % need two weeks to see pattern
      if on_off(4,i) > on_off (2,i)
        mid = (on_off(2,i) + on_off(4, i)) / 2;
      else
        mid = mod1 ((on_off(2,i) + on_off(4, i) + size (cv,1)) / 2, size (cv,1));
      end
      mid = mod1 (round (mid), size (cv, 1));
      range = on_off(1,i):on_off(3,i)-1;

      off_time = mod1 (round(on_off(4,i)), size(cv, 1));
      on_time  = mod1 (round(on_off(2,i)), size(cv, 1));
      before_off = mod1 (off_time-1, size(cv, 1));
      before_on  = mod1 (on_time -1, size(cv, 1));

      % jump in the expected direction on each day
      steps_off = cv(before_off,range) - cv(off_time, range);
      steps_on = cv(on_time, range) - cv(before_on, range);

      % mean jump on each day-of-week
      DoW_off = mean7 (steps_off);
      DoW_on  = mean7 (steps_on);

      % Seems to be no jump.  (Choose threshold more carefully!  Should depend
      % on noise, whether other rectangles have day-of-week dependence etc.
      % Sufficient statistic should also include more times either side of jump.)
      no_jump_off = (DoW_off < power/3);
      no_jump_on  = (DoW_on  < power/3);

%      % Alternative way of estimating alt_days.
%      % This identifies more (too many?) rectangles.
%      % It is probably better to consider turn-on and turn-off together.
%      DoW_pairs = DoW_on + DoW_on([2:7, 1]);
%      [m, pair] = min (DoW_pairs);
%      if m < power/2
%        days1 = [pair, 1+mod(pair, 7)];
%        no_jump_on = false (size (no_jump_on));
%        no_jump_on(days1) = true;
%      end
%      DoW_pairs = DoW_off + DoW_off([2:7, 1]);
%      [m, pair] = min (DoW_pairs);
%      if m < power/2
%        days2 = [pair, 1+mod(pair, 7)];
%        no_jump_off = false (size (no_jump_off));
%        no_jump_off(days2) = true;
%      end

      days1 = find (no_jump_on);
      days2 = find (no_jump_off);
      burst = rectangles(i).burst;
      alt_is_off = false;   % alt days have zero power (during burst)
      if length (days1) == 2 && length (days2) == 2 && all (days1 == days2)
        x = cv (burst, range([days1(1):7:end, days1(2):7:end]));
            % TODO: should look at the *diff* between alt/no-alt days
        if all (mean (x, 2) < 0.5 * power)
          days_off = mod1 (days1 + range(1) - 1, 7);
          time_off = mid;
          conf_off = 2;
          days_on  = days_off;
          time_on  = time_off;
          conf_on  = conf_off;
          alt_is_off = true;
        %else
          %keyboard;
        end
      end
      if ~alt_is_off
        [days_off, time_off, conf_off] = alt_jump (no_jump_off, DoW_off, on_off (4, i), mid,  power, range, cv);
        [days_on,  time_on,  conf_on]  = alt_jump (no_jump_on,  DoW_on,  on_off (2, i), mid, -power, range, cv);
        if ~intseteq (days_on, days_off)
          if conf_off > conf_on
            days_on = [];
          else
            days_off = [];
            if conf_off + conf_on == 0
              days_on = [];
            end
          end
        end
      end

      d_off{i} = days_off;
      d_on {i} = days_on;
      t_off(i) = time_off;
      t_on (i) = time_on;
      c_off(i) = conf_off;
      c_on (i) = conf_on;
    end
  end

  noise = (abs (t_on' - on_off(2,:)) < 1) ...
          | (abs (t_off' - on_off(4,:)) < 1);
  tmp_off = t_off;
  tmp_on  = t_on;
  tmp_off(t_off == 0) = on_off(4, t_off == 0);
  tmp_on (t_on  == 0) = on_off(2, t_on  == 0);
  noise = noise | ((tmp_off - tmp_on)' > abs (on_off(4,:) - on_off(2,:)) .^ 1.5);

  c_on (noise) = 0;
  c_off(noise) = 0;
  d_on (noise) = {[]};
  d_off(noise) = {[]};

  % Check for consistency of days-of-week for overlapping rectangles
  matches = 0;    % intervals on which alt times are present and consistent
  m_days = 0;     % sum of lengths of intervals during which alt times present
  half_matches = 0;  % intervals on which one rectangle has alternative times
  hm_days = 0;
  mismatches = 0;
  mm_days = 0;
  none = 0;
  n_days = 0;
  %fprintf ('\n');
  for e = 1:T-1
    days = timer.start(e+1) - timer.start(e);
    r = timer.rects{e};
    %fprintf ('Interval %g to %g: ', timer.start(e), timer.start(e+1));
    % TODO: If overlap, truncate alt times (?)
    if length (r) == 2
      if intseteq ([d_off{r(1)}, d_on{r(1)}], [d_off{r(2)}, d_on{r(2)}])
        if ~isempty ([d_off{r(1)}, d_on{r(1)}])
          matches = matches+max([c_on(r(1:2)); c_off(r(1:2))]);
          m_days = m_days + days;
          %fprintf ('*** %d and %d have days ', r(1), r(2));
          %fprintf (' %d', [d_off{r(1)}, d_on{r(1)}]);
        else
          none = none + 1;
          n_days = n_days + days;
          %fprintf ('%d and %d are both empty', r(1), r(2));
        end
      elseif isempty ([d_off{r(1)}, d_on{r(1)}])
        half_matches = half_matches + max(c_on(r(1)), c_off(r(1)));
        hm_days = hm_days + days;
        %fprintf ('* %d has ', r(2));
        %fprintf (' %d', [d_off{r(2)}, d_on{r(2)}]);
      elseif isempty ([d_off{r(2)}, d_on{r(2)}])
        half_matches = half_matches + max(c_on(r(2)), c_off(r(2)));
        hm_days = hm_days + days;
        %fprintf ('* %d has ', r(1));
        %fprintf (' %d', [d_off{r(1)}, d_on{r(1)}]);
      else
        %fprintf ('%d has ', r(1));
        %fprintf (' %d', [d_off{r(1)}, d_on{r(1)}]);
        %fprintf (' but %d has ', r(2));
        %fprintf (' %d', [d_off{r(2)}, d_on{r(2)}]);
        mismatches = mismatches+1;
        mm_days = mm_days + days;
      end
    elseif length (r) == 1
      if ~isempty ([d_off{r(1)}, d_on{r(1)}])
        matches = matches + 1;
        m_days = m_days + days;
      else
        none = none + 1;
        n_days = n_days + days;
      end
    else
      %fprintf ('has %d rectangles', length (r));
    end
    %fprintf ('\n');
  end
  %fprintf ('matches %d half_matches %d mismatches %d none %d\n', ...
  %          matches, half_matches, mismatches, none);
  % If sufficiently many matches
  %   Create "alt" rectangles (or turn rectangles into alt rectangles)
  %        alt field with start/end time and days
  %   Recalculate reliability of each, and overall power.

  if matches + half_matches > 2 * (mismatches + none) ...
     || m_days + hm_days > 2 * (mm_days + n_days) ...
     || (matches > 1.5 * (mismatches+1) && m_days > 1.5 * mm_days && m_days > 60)    % Deem matches to be "real"
    for i = idx
      if on_off(3,i) - on_off(1,i) > 14   % pattern looked for explicitly
        if ~isempty (d_off{i}) || ~isempty (d_on{i})
          if ~isempty (d_off{i})
            rectangles(i).alt_days(1:length (d_off{i}), 1) = d_off{i};
            rectangles(i).alt_hrs(2) = t_off(i);
          else
            rectangles(i).alt_hrs(2) = on_off(4, i);
          end
          if ~isempty (d_on{i})
            rectangles(i).alt_days(1:length (d_on{i}), 1) = d_on{i};
            rectangles(i).alt_hrs(1) = t_on(i);
          else
            rectangles(i).alt_hrs(1) = on_off(2, i);
          end
        end
      else
        % Eventually: If this overlaps with a rectangle with alt times,
        % check those days to see if this rectangle does too.
      end
    end
  elseif any (c_on > 1)     % Hack: if we're "confident", force these
    for i = find (c_on')    % rectangles to have alt_days. Currently only set
      rectangles(i).alt_days(1:length (d_on{i}), 1) = d_on{i}; % when
      rectangles(i).alt_hrs = [t_on(i); t_off(i)]; % alt days are always off.
    end
  end

  % Check for alt_days in rectangles overlapping those with alt_days
  if isfield (rectangles, 'alt_days')
    % For each rectangle, record which alt_days have been tested
    checked = false (length (rectangles), 7);

    has_alt_day = double (~cellfun (@isempty, { rectangles.alt_days }));
    has_alt_len = sum (has_alt_day);
    has_alt_day(1:has_alt_len) = find (has_alt_day);
    i = 1;
    while i < has_alt_len
      ii = has_alt_day(i);
      overlap = find (on_off(3,:) > on_off(1,ii) + 21 ...
                    & on_off(1,:) < on_off(3,ii) - 21 ...
                    & on_off(3,:) - on_off(1,:) > 21 ...
                    & 1:size(on_off,2) ~= ii);
      alt_days = zeros (1, size(cv,2));
      alt_days ([rectangles(ii).alt_days(1):7:size(cv,2), ...
                 rectangles(ii).alt_days(2):7:size(cv,2)]) = true;

      ro = rectangles(overlap);
      no_alt = cellfun (@isempty, { ro.alt_days });
      if any (no_alt)
        o_no = overlap(no_alt);
        r = rectangles([o_no, ii]);
        roo = [r.on_off];
        st = min (roo(1,:));
        en = max (roo(3,:));
        days = zeros (1, size(cv, 2));
        days(st:en-1) = true;
        days = find (days & alt_days);
        [~, new_st] = max (bsxfun (@ge, days,           roo(1,:)')');
        [~, new_en] = max (bsxfun (@le, days(end:-1:1), roo(3,:)' - 1)');
        new_en = length (days) - new_en + 2;
        new_oo = num2cell ([new_st(:), roo(2,:)', new_en(:), roo(4,:)']', 1);
        [r.on_off] = deal (new_oo{:});
        r(end).on_off([2,4]) = r(end).alt_hrs;
        [r, kept] = tweak_on_off_pwr (r, cv(:, days), r(1).power, ...
                                      issolar, 1:length (days));
        % if tweak_on_off_pwr deleted rectangles, delete from o_no too
        o_no = o_no(kept(kept <= length (o_no)));

        for j = 1:length (o_no)
          jj = o_no(j);
          x = (abs (r(j).on_off([2,4]) - rectangles(jj).on_off([2,4])) > 1);
          if any (x)
            has_alt_len = has_alt_len + 1;
            has_alt_day(has_alt_len) = jj;
            rectangles(jj).alt_days = rectangles(ii).alt_days;
            rectangles(jj).alt_hrs = rectangles(jj).on_off([2,4]);
            if x(1)
              rectangles(jj).alt_hrs(1) = r(j).on_off(2);
            end
            if x(2)
              rectangles(jj).alt_hrs(2) = r(j).on_off(4);
            end
          end
        end
%          TODO:
%          if overlap is also alt_day, with incompatible days
%            if overlap is small, truncate.
%            if overlap is large, choose which alt_day is better
      end
      i = i + 1;
    end
  end
end

function m = mean7 (vec)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Find the seven means 1:7:end, 2:7:end, 3:7:end etc.

 %  pad = 7 - mod (length (vec), 7);
 %  if pad == 7
 %    pad = 0;
 %  end
 %  rsh = reshape (vec, [7, length(vec)+pad]);
 %  m(1:(7-pad)) = mean (rsh(1:(7-pad),:), 2);
 %  m(8-pad:7) = mean (rsh(8-pad:7, 1:end-1), 2);
  if (length (vec) < 7)
    m = vec;
  else
    m(7) = mean (vec(7:7:end));
    m(6) = mean (vec(6:7:end));
    m(5) = mean (vec(5:7:end));
    m(4) = mean (vec(4:7:end));
    m(3) = mean (vec(3:7:end));
    m(2) = mean (vec(2:7:end));
    m(1) = mean (vec(1:7:end));
  end
end

function [days, time, confidence] = alt_jump (no_jump, DoW, main_time, known_on, power, range, cv)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Find turn on/off time on days that don't match the standard time
 % no_jump = logical vector of which days don't have a big step at main_time
 % DoW = vector of mean jump size on each day of week
 % main_time = time of current estimate of jump
 % known_on = a time that the device is "known" to be on
 % range = list of days in this rectangle
 % Days are 1 = 1st,8th... day of valid_days, 7 = 7th,14th... day of valid_days
  confidence = 0;
  reduce_confidence_factor = 1; % If we "nudge" a choice, reduce this to (0,1)
  time = 0;
  days = find (no_jump);
  n_alt = length (days);
  if n_alt > 0
    if n_alt == 1
      tmp = DoW;
      tmp(no_jump) = max (DoW);
      [~, next_smallest] = min (tmp);
      if hr_diff (days, next_smallest, 7) == 1
        reduce_confidence_factor = 0.7;
        days = [days, next_smallest];
        n_alt = 2;
      end
    end
    match = zeros (size (days));
    for j = 1:length (days)
      x = mean (cv (:, range(days(j):7:end)), 2);
      match(j) = match_jump (x, -known_on, power);
    end
    % if 2 days with no jump at off_time have "better" jumps at similar times
    % to each other, but different from off_time, be probably alt timer
    if n_alt == 2 && any (diff (days) == [-1, -6, 1, 6])
      if hr_diff (match(1), match (2), size(cv,1)) < 2 ...
         &&  hr_diff ((match(1)+match(2))/2, main_time, size(cv,1)) > 1
          confidence = 1;
      else
        x = mean(cv(:,range([days(1):7:end,days(2):7:end])), 2);
        both_match = match_jump (x, -known_on, power);
        if min (abs (match - both_match)) < 1
          confidence = 0.5;
        end
      end
    % What other cases suggest alternate timer?
    else
      [~, most] = min (DoW + DoW ([2:end, 1]));
      diffs = hr_diff (match, main_time, size(cv,1));
      if sum (diffs > 1.5) == 2
        days = days(diffs > 1.5);
        if all (diff (days) ~= [-1, -6, 1, 6])
          days = [];
        elseif any (most == days) && any (mod1 (most+1, 7) == days)
          confidence = 0.9;
        else
          confidence = 0.5;
        end
      else
        m1 = find (days == most);
        m2 = find (days == mod1 (most+1, 7));
        if ~isempty (m1) && ~isempty (m2)
          if hr_diff (match (m1), match (m2), size (cv,1)) < 2
            days = [most, mod1(most+1,7)];
            if hr_diff ((match(1)+match(2))/2, main_time, size(cv,1)) > 1
              confidence = 0.9;
              % TODO -- find which other days also match
            end
          else
            %keyboard
            days = [];
          end
        else
          days = [];
          %keyboard
        end
      end
    end

    if ~isempty (days)
      % Match difference between normal and alt days
      % 1. difference of means
      y = 0;
      x = 0;
      for d = 1:7
        if any (d == days)
          x = x + mean (cv (:, range(d:7:end)), 2);
        else
          y = y + mean (cv (:, range(d:7:end)), 2);
        end
      end
      x = x / length (days);
      y = y / (7- length (days));

      % 2. First guess, to see whether alt times are earlier or later
      time = match_jump (x, -known_on, power);

      % 3. Match, if time seems to be different
      if hr_diff (time, main_time, size(cv,1)) > 1
        if (power < 0 && inorder (time, main_time, known_on)) ...
            || (power > 0 && inorder (known_on, main_time, time))
          time = match_jump (x - y, round (main_time), power);
        else
          time = match_jump (y - x, round (main_time), -power);
        end
        % TODO: If large difference between nominal and alt times, reduce conf
        confidence = confidence * reduce_confidence_factor;
      end
      days = mod1 (days + range(1) - 1, 7);
    end
  end
end


function yes = inorder (st, mid, en)
  % True if time  mid  is in the interval  st:en, mod  x  for some x>max(en,st)
  if st < en
    yes = st < mid && mid < en;
  else
    yes = st < mid || mid < en;
  end
end

function yes = between (A, mid, B, modulo_2)
  % Consider the shorter of the two "intervals" modulo 2*modulo_2
  % bounded by A and B.
  % Return true if mid in in that interval.
  if (abs (A - B) < modulo_2)     %usual case, not crossing midnight
    if A > B
      yes = (B < mid && mid < A);
    else
      yes = (A < mid && mid < B);
    end
  else
    if A < B
      yes = (B < mid || mid < A);
    else
      yes = (A < mid || mid < B);
    end
  end
end

function yes = int_overlap (st1, en1, st2, en2)
  % Return true if intervals st1:en1 and st2:en2 overlap (modulo size(cv,1)).
  if max (st1, st2) < min (en1, en2)
    yes = true;
  elseif st1 > en1
    yes = (max (st2, en2) >= st1 || min (st2, en2) < en1);
  elseif st2 > en2
    yes = (en1 >= st2 || st1 < en2);
  else
    yes = false;
  end
end
function interval = wrap (st, en, last)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Logical interval of time-of-day starting at st, ending at en, with
 % last  samples per day.
  interval = false (1,last);
  if st < en
    interval (st:en) = 1;
  else
    interval ([1:en, st:last]) = 1;
  end
end

function m = mod1 (value, modulus)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Modulo operator, returning values in [1, modulus], not [0, modulus-1].
  m = mod (value-1, modulus) + 1;
end

function d = hr_diff (x, y, modulus)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % hr_diff (x, y, modulus) is the difference, modulo modulus, between x and y
  d = abs (x - y);
  d = min (d, modulus - d);
end

function eq = intseteq (x, y)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % % True if sets of integers x and y are equal, false otherwise.
 % % Requires memory equal to the maximum value in x and y minus the minimum.
 %  eq = false;
 %  if length (x) == length (y) && max (x) == max (y)
 %    m = min ([x(:); y(:)]) - 1;
 %    a (x(:) - m) = true;
 %    eq = all (a(y(:) - m));
 %  end
  eq = isempty (setxor (x, y));
end

function rectangles = powers_reliabilities_DoW (rectangles, cv, cor)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % On/off:  Check lowest off times.
  %          If run of > "10" are "low", and step is "const",
  %             test hypothesis that this (2-slot?) step is the true step
  %          If no run or test fails
  %             take (weighted?) average over all junks
  % start/end days:
  %          If "constant" for three days prior and three days after
  %               and cor_pool predicts a jump
  %            take (raw) average over diff in days
  %          Take weighted average over all edges
  for r = 1:length (rectangles)
    on_off = rectangles(r).on_off;
    on = round (on_off(2));
    off = round (on_off(4));
    st = on_off (1);
    en = on_off (3);
    all_days = st:en-1;
    days = all_days;

    alt_days = [];
    alt_on = on;
    alt_off = off;

    if isfield (rectangles(r), 'alt_days') && sum (rectangles(r).alt_days)
      idx = ismember (mod1 (days, 7), rectangles(r).alt_days);
      alt_days = all_days(idx);
      days = all_days(~idx);
      alt_on  = round (rectangles(r).alt_hrs(1));
      alt_off = round (rectangles(r).alt_hrs(2));
    end

    [~, idx] = sort ([days, alt_days]);
    before = [index(cv, on - 1, days), index(cv, alt_on - 1, alt_days)];
    before = before(idx);
    after = [index(cv, on + 1, days), index(cv, alt_on + 1, alt_days)];
    after = after(idx);
    step = after - before;

    %%TODO: This bit needs work!
  end
end

function rectangles = set_missed (rectangles, power, cv)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Find days on which it seems the pump was off, considering alternate times
  for r = 1:length (rectangles)
    on_off = rectangles(r).on_off;
    st = on_off(1);
    en = on_off(3);
    all_days = st:en-1;
    days = all_days;
    alt_days = [];
    if isfield (rectangles(r), 'alt_days') && sum (rectangles(r).alt_days)
      idx = ismember (mod1 (days, 7), rectangles(r).alt_days);
      alt_days = all_days(idx);
      days = all_days(~idx);
    end

    ht = power;
    missed = false (size (all_days));

    if length (rectangles(r).burst) > 1
      b = rectangles(r).burst(1:end-1);
    else
      b = ceil (rectangles(r).on_off(2));
      if b > size (cv, 1)
        b = b - size (cv, 1);
      end
    end
    slice = cv(b,days);
    missed(days - st + 1) = (min (slice,[],1) < 0.9 * ht);

    slice = cv(b,alt_days);
    missed(alt_days - st + 1) = (min (slice,[],1) < 0.9 * ht);
    rectangles(r).missed = missed;
  end
end

function [rectangles, others] = set_missed_precise (rectangles, power, cv, cc, runs, issolar)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Find days on which it seems the pump was off, considering alternate times,
 % taking note of days that seem to have the wrong on/off times instead.
 others = [];
 on_off = [rectangles.on_off];
 [pre_on,  mid_on,  post_on,  a_on]  = edge_neighbours (on_off(2,:), size (cv, 1), true (1,size(on_off,2)));
 [pre_off, mid_off, post_off, a_off] = edge_neighbours (on_off(4,:), size (cv, 1), true (1,size(on_off,2)));

  frac_on  = on_off(2,:)   - mid_on;
  frac_off = on_off(4,:)-1 - mid_off;

  for r = 1:length (rectangles)
    % on_off = rectangles(r).on_off;
    st = on_off(1, r);
    en = on_off(3, r);
    all_days = st:en-1;
    days = all_days;
    alt_days = [];
    if isfield (rectangles(r), 'alt_days') && sum (rectangles(r).alt_days)
      idx = ismember (mod1 (days, 7), rectangles(r).alt_days);
      alt_days = all_days(idx);
      days = all_days(~idx);
    end

    ht = rectangles(r).power;
    missed = false (size (all_days));
    below = missed;
    no_step = missed;

    if length (rectangles(r).burst) > 1
      b = rectangles(r).burst(1:end-1);
    else
      b = ceil (rectangles(r).on_off(2));
      if b > size (cv, 1)
        b = b - size (cv, 1);
      end
    end
    slice  = cv(b,days);
    slice2 = cc(b,days);
    dy = days - st + 1;
    below(dy)  = (min (slice, [],1) < 0.9 * ht);
    below2(dy) = (min (slice2,[],1) < 0.9 * ht);
    no_step_on  = below;
    no_step_off = below;

    if any (below(dy))
      % no_step = ...;    % How can we estimate the lack of step?
      % If run of >4 missing, look for sub-burst.
      % If no_step, but high on both sides, look for extended burst
      % Look for matched sub-burst and extended burst
      % Look for being ~ht below nearest non-missed neighbour
      % Look for violation of frational step

      raw_on = cv(post_on(r), st:en-1) - cv(pre_on(r),  st:en-1);
      mf_on  = medfilt1 (raw_on, 5);
      no_step_on = (mf_on < 0.5 * ht);
      no_raw_on = (raw_on < 0.5 * ht);
      raw_off = cv(pre_off(r), st:en-1) - cv(post_off(r),st:en-1);
      mf_off = medfilt1 (raw_off, 5);
      no_step_off = (mf_off < 0.5 * ht);
      no_raw_off = (raw_off < 0.5 * ht);
      % process mid
      if frac_on(r) < 0.5
        raw_on = cv(mid_on(r), st:en-1) - cv(pre_on(r), st:en-1);
      else
        raw_on = cv(post_on(r),  st:en-1) - cv(mid_on(r), st:en-1);
        frac_on = 1-frac_on;
      end
      mf_on = medfilt1 (raw_on, 5);
      no_mid_on = (mf_on < 0.5 * (1 - frac_on(r)) * ht);
      no_on = no_step_on | no_mid_on;
      no_raw_on = no_raw_on | (raw_on < 0.5 * (1 - frac_on(r)) * ht);

      if frac_off(r) > -0.5
        raw_off = cv(mid_off(r), st:en-1) - cv(post_off(r), st:en-1);
      else
        raw_off = cv(pre_off(r),  st:en-1) - cv(mid_off(r), st:en-1);
        frac_off(r) = -1-frac_off(r);
      end
      mf_off = medfilt1 (raw_off, 5);
      no_mid_off = (mf_off < 0.5 * (1 + frac_off(r)) * ht);
      no_off = no_step_off | no_mid_off;
      no_raw_off = no_raw_off | (raw_off < 0.5 * (1 - frac_off(r)) * ht);

      missed = find (median (double ([below(dy);
                                      no_raw_on(dy);
                                      no_on(dy);
                                      no_raw_off(dy);
                                      no_off(dy)])));
      miss_runs= dy(ranges (missed)) + st - 1;
      if isempty (miss_runs)
        continue
      end
      if size(miss_runs, 1) == 1
        miss_runs = miss_runs';
      end
      miss_runs = miss_runs(:, miss_runs(2,:) - miss_runs(1,:) >= 3);

      % TODO: Compare mean of min (width or run, 5) vs
      %               mean of min (width of non-run, 5)
      %low_left = slice(:, 2:end) - slice(:, 1:end-1);
      %low_right = sum (low_left >  0.7 * ht);
      %low_left  = sum (low_left < -0.7 * ht);
      if issolar
        % TODO: be more tolerant of not meeting ht.
      else
        % take median of miss measures?
      end


      for i = 1:size(miss_runs, 2);
        overlap_runs = runs(:,1) <= miss_runs(2,i) & runs(:,2) > miss_runs(1,i);
        if on_off(4, r) > on_off(2, r)
          overlap_runs = overlap_runs & runs(:,3) > on_off(2, r) & runs(:,3) < on_off(4, r);
        else
          overlap_runs = overlap_runs & (runs(:,3) < on_off(2, r) | runs(:,3) < on_off(4,r));
        end
        overlap_runs = find (overlap_runs);
        for j = 1:length (overlap_runs)
          run = runs(overlap_runs(j),:);
          s = max (run(1), miss_runs(1,i));
          e = min (run(2), miss_runs(2,i) + 1);
          if (e - s >= 3)
            if (e - s > 4)    % Allow for mismatch in timer day
              e = e - 1;
              s = s + 1;
            end
            new_r = find_rect (cv, s, e, rectangles(r).power, r, rectangles);
            % TODO: find quality of new:
            %       - quality of jumps
            %       - jumps down relative to day before and day after run
            %       - match length of rectangle we're interrupting,
            %         or on or off.
            %others = [others, new_r];
          end
        end
      end
    end
    % TODO: compare missed days with those in overlap.
    % TODO: Trim overlapping  new_r

    slice = cv(b,alt_days);
    below(alt_days - st + 1) = (min (slice,[],1) < 0.9 * ht);

    rectangles(r).missed = missed;

    if ~isempty (others)
      figure(8); show_rectangles ([rectangles, others], cv);
    end
  end
end

function vals = index (cv, x, y)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Emulate vals = cv (x, y), where x is hours and y is days, so that cv is
  % really a linear array.  For x < 1 or x > size (cv, 1), wrap x and shift y.
  if min (length (x), length (y)) == 1
    idx = x + y * size (cv, 1);
    if any (idx < 1 | idx > length (cv (:)))
      idx (idx < 1) = mod1 (idx (idx < 1), size (cv, 1));
      idx (idx > length (cv (:))) = mod1 (idx (idx > length (cv (:))), size (cv, 1)) ...
                                   + (size (cv,2) - 1) * size (cv, 1);
    end
    vals = cv (idx);
  else
    after  = x > size (cv, 1);
    before = x < 1;
    others = ~before & ~after;
    vals (others,:) = cv(x(others), y);
    vals (before,:) = cv(x(before) + size (cv, 1), max (y-1, 1));
    vals (after,:)  = cv(x(after)  - size (cv, 1), min (y+1, size(cv, 2)));
  end
end

function vals = index1 (cv, x1, x2, y)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Emulate vals = cv (x1:x2, y), where x is hours and y is days, so that cv is
  % really a linear array.  For x < 1 or x > size (cv, 1), wrap x and shift y.
  if x2 > x1
    x = x1:x2;
  else
    x = x1:(x2 + size (cv, 1));
  end

  after  = x > size (cv, 1);
  before = x < 1;
  others = ~before & ~after;
  vals (others,:) = cv(x(others), y);
  vals (before,:) = cv(x(before) + size (cv, 1), max (y-1, 1));
  vals (after,:)  = cv(x(after)  - size (cv, 1), min (y+1, size(cv, 2)));
end

%function [timer, rectangles] = ground_truth (cv, timer, rectangles)
%
%end

function [timer, rectangles, ok] = pool_ground_truth_old (cv, timer, rectangles)
  dims = get (0, 'ScreenSize');
  if isempty (timer)
    dims = dims/2;
  end
  f = figure ('Visible', 'off', 'Position', [0, 0.5*dims(4), 0.8*dims(3), 0.8*dims(4)]);

  ha = axes ('Units', 'pixels', 'Units', 'normalized', ...
             'Position', [0.05, 0.05, 0.7, 0.9]);

  imagesc (cv);

  ok = true;
  undo_stack = {};
  days_to_hide = [];
  rectangles_visible = true;
  draw_rectangles;     % will set rectangles_visible = y

  %showHide =
  uicontrol ('Style', 'pushbutton', 'String', 'Show rectangles', ...
             'Units', 'normalized', ...
             'Position', [0.8, 0.9, 0.15, 0.04], ...
             'Callback', {@draw_or_clear_rectangles});

  %del_button =
  uicontrol ('Style', 'pushbutton', 'String', 'Delete rectangle',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.8, 0.15, 0.04], ...
             'Callback', {@del});

  %add_button =
  uicontrol ('Style', 'pushbutton', 'String', 'Add rectangle',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.7, 0.15, 0.04], ...
             'Callback', {@add});

  uicontrol ('Style', 'pushbutton', 'String', 'Merge',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.6, 0.15, 0.04], ...
             'Callback', {@merge_selected});

  uicontrol ('Style', 'pushbutton', 'String', 'Split rectangle',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.5, 0.15, 0.04], ...
             'Callback', {@set_splitting});

  uicontrol ('Style', 'pushbutton', 'String', 'Alt Days',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.4, 0.15, 0.04], ...
             'Callback', {@toggle_alt_days});

  uicontrol ('Style', 'pushbutton', 'String', 'Undo',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.25, 0.15, 0.04], ...
             'Callback', {@undo});
  %done_button =
  uicontrol ('Style', 'pushbutton', 'String', 'Done',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.10, 0.15, 0.04], ...
             'Callback', {@done});

  uicontrol ('Style', 'pushbutton', 'String', 'Cancel',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.02, 0.15, 0.04], ...
             'Callback', {@cancel});

  set (f, 'WindowButtonDownFcn', @start_select_or_move);
  set (f, 'WindowButtonUpFcn',   @final_select_or_move);

  % initialize for use below
  selected_rect = [];
  downtime = 0;
  coordinates = 0;

  if (~exist ('OCTAVE_VERSION', 'builtin'))
    movegui (f, 'center');
  end
  set (f, 'Visible', 'on');

  uiwait (f);

  function draw_or_clear_rectangles (~, ~)
    if ~rectangles_visible
      draw_rectangles;
    else
      %plot_handle =
      tmp = cv;
      for i = days_to_hide
        tmp(:,i:7:end) = 0;
      end
      imagesc (tmp);
      rectangles_visible = false;
    end
  end

  function draw_rectangles (~, ~)
    show_rectangles (rectangles, cv, days_to_hide);
    rectangles_visible = true;
  end

  function set_splitting (~, ~)
    set (f, 'WindowButtonUpFcn', @split_rectangle);
  end

  function split_rectangle (~,~)
    % split a rectangle into two, leaving a gap of size determined by
    % the mouse movement.  The split is horizontal if most mouse movement
    % (as a fraction of the size of the rectangle) is vertical and vice
    % versa.
    set (f, 'WindowButtonUpFcn', @final_select_or_move); % restore default

    old_coords = coordinates;
    coordinates_now = get (ha, 'CurrentPoint');
    coordinates_now = coordinates_now (1, 1:2);
    move = coordinates_now - old_coords;
    [r, found] = get_rectangle_from_coords ();  % clobbers 'coordinates'
    if ~found
      return
    end
    rectangles(end+1) = rectangles (r);
    add_undo ({ 'replace', r, rectangles(r) });
    add_undo ({ 'delete', length(rectangles), [] });

    width = rectangles(r).on_off(3) - rectangles(r).on_off(1);
    height = rectangles(r).on_off(4) - rectangles(r).on_off(2);
    if move(1) / width >  move(2) / height
      % vertical split
      rectangles(end).on_off(4) = old_coords (2);
      rectangles(r).on_off(2) = coordinates_now (2);
    else
      % horizontal split
      [en, st] = align_with_existing (old_coords, coordinates_now);
      rectangles(end).on_off(3) = round (en);
      rectangles(r).on_off(1) = round (st);
    end
    draw_rectangles;
  end

  function del (~, ~)
    sr = -sort (-selected_rect);
    e = length (undo_stack);
    undo_stack{e+length (sr)} = [];   % pre-allocate to avoid Matlab warning
    for i = 1:length (sr)
      r = sr(i);
      undo_action = {'insert', r, rectangles(r)};
      undo_stack{e + i} = undo_action;
      rectangles = rectangles ([1:r-1, r+1:end]);
    end
    selected_rect = [];

    draw_rectangles;
  end

  function yes = inorder (st, mid, en)
    % Copied from has_pool_pump
    % True if time  mid  is in the interval  st:en, mod  x  for some x>max(en,st)
    mask = (st < en);
    if length (mid) == 1 && length (mask) > 1
      mid = repmat (mid, size (mask));
    end
    yes(mask)  = st(mask)  < mid(mask)  & mid(mask)  < en(mask);
    yes(~mask) = st(~mask) < mid(~mask) | mid(~mask) < en(~mask);
  end

  function select_rectangle (object_handle, ~)
    [r, found] = get_rectangle_from_coords;

    if found
      modifiers = get_modifiers (object_handle);
      if ismember('shift',modifiers) || ismember('control', modifiers)
        selected_rect(end+1) = r;
      else
        draw_rectangles;
        selected_rect = r;
      end
      draw_as_selected (rectangles(r).on_off);
    else
      selected_rect = [];
      draw_rectangles;
    end
  end

  function draw_as_selected (oo)
    if isempty (oo)
      return
    end
    oo([1,3],:) = oo([1,3],:) - 0.5;

    s = (oo(4,:) > oo(2,:));   % "simple" case
    if any (s)
      x = [oo(1,s); oo(1,s); oo(3,s); oo(3,s)];
      y = [oo(2,s); oo(4,s); oo(4,s); oo(2,s)]-0.5;
      patch (x, y, 'w');
    end

    s = ~s;   % "split over days" case
    if any (s)
      one = ones (1, sum (s));
      x = [oo(1,s); oo(1,s); oo(3,s); oo(3,s)];
      y1 = [(size(cv, 1)+1)*one; oo(2,s); oo(2,s); (size(cv, 1)+1)*one];
      y2 = [one; oo(4,s); oo(4,s); one];
      patch ([x, x+1], [y1, y2] - 0.5, 'w');
    end

  end

  function [r, found] = get_rectangle_from_coords
    coordinates = get (ha, 'CurrentPoint');
    x = coordinates (1, 1);
    y = coordinates (1, 2);
    found = false;
    for r = 1:length (rectangles)
      if rectangles(r).on_off(1) <= x && rectangles(r).on_off(3) >= x ...
          && inorder (rectangles(r).on_off(2), y, rectangles(r).on_off(4))
        found = true;
        break
      end
    end
  end

  function merge_selected (~, ~)
    % Form one rectangle covering the extent covered by the highlighted
    % rectangles.
    % TODO: allow merging more than two at a time.  Need to decided
    %       horizontal or vertical merge, and sort into geometric order.
    if length (selected_rect) >= 2
      s = sort (selected_rect);
      for i = length(s)-1:-1:1
        merge_rectangles (s(i), s(i+1));
      end
    end
    selected_rect = [];
    draw_rectangles;
  end

  function merge_rectangles (r1, r2)
    add_undo ({ 'replace', r1, rectangles(r1) });
    add_undo ({ 'insert', r2, rectangles(r2) });

    on_off1 = rectangles(r1).on_off;
    on_off2 = rectangles(r2).on_off;
    rectangles(r2) = [];

    if on_off1(2) < on_off1(4)    % doesn't span midnight
      rectangles(r1).on_off(2) = min (on_off1(2), on_off2(2));
      rectangles(r1).on_off(4) = max (on_off1(4), on_off2(4));
    else
      rectangles(r1).on_off(2) = max (on_off1(2), on_off2(2));
      rectangles(r1).on_off(4) = min (on_off1(4), on_off2(4));
    end
    rectangles(r1).on_off(1) = min (on_off1(1), on_off2(1));
    rectangles(r1).on_off(3) = max (on_off1(3), on_off2(3));
  end

  function start_select_or_move (~, ~)
    coordinates = get (ha, 'CurrentPoint');
    coordinates = coordinates (1, 1:2);
    downtime = toc;
  end

  function modifiers = get_modifiers (h)
    s = get (h, 'selectiontype');
    if (strcmp (s, 'extend'))
      modifiers = {'shift'};
    elseif (strcmp (s, 'alt'))
      modifiers = {'control'};
    else
      modifiers = {};
    end
  end

  function final_select_or_move (object_handle, event_data)
    if isempty (rectangles)
      return
    end
    modifiers = get_modifiers (object_handle);

    coordinates_now = get (ha, 'CurrentPoint');
    coordinates_now = coordinates_now (1, 1:2);
    move = coordinates_now - coordinates;
    move(1) = round (move(1));
    scaled = [move(1), 5*move(2)];

    on_off = [rectangles(:).on_off];
    mn = min (coordinates_now, coordinates);
    mx = max (coordinates_now, coordinates);

    % Rectangles entirely within dragged region
    r = on_off(1,:) > mn(1) & on_off(2,:) > mn (2) ...
        & on_off(3,:) < mx(1) & on_off(4,:) < mx(2) ...
        & on_off(4,:) > on_off(2,:);

    if toc - downtime < 0.5 && norm (scaled) < 3
      % assume a click, not a drag
      select_rectangle (object_handle, event_data);
    elseif any (r)
      % select all rectangles in the range
      if ismember ('shift', modifiers)
        selected_rect = unique ([selected_rect, find(r)]);
      else
        selected_rect = find (r);
      end
      draw_as_selected ([rectangles(selected_rect).on_off]);
    else
      % Assume a drag.
      % Use both direction of drag and proximity to edges to work out
      % which edge was being dragged.

      % vertical move => horizontal edge
      seems_horiz = (abs (scaled(1)) < abs (scaled(2)));

      vert_candidates = find (inorder(on_off(2,:) - 0.5, coordinates(2), ...
                                      on_off(4,:) - 0.5));
      edges = [on_off(1, vert_candidates), ...
               on_off(3, vert_candidates)];
      [~, vert_closest] = min (abs (edges - coordinates(1)));
      if vert_closest > length (vert_candidates)
        vert_closest = vert_candidates(vert_closest ...
                                       - length (vert_candidates));
        vert_edge = 3;
      else
        vert_closest = vert_candidates(vert_closest);
        vert_edge = 1;
      end

      horiz_candidates = find (on_off(1,:) < coordinates(1) ...
                               & coordinates(1) < on_off(3, :));
      edges = [on_off(2, horiz_candidates), ...
               on_off(4, horiz_candidates)];
      if isfield (rectangles(1), 'alt_days')
        alt_cand = ~cellfun (@isempty, {rectangles(horiz_candidates).alt_days});
        hrs = [rectangles(horiz_candidates(alt_cand)).alt_hrs];
        edges = [edges, hrs];
      end

      [~, horiz_closest] = min (abs (edges - coordinates(2)));
      if horiz_closest > length (horiz_candidates)
        if horiz_closest <= 2*length (horiz_candidates)
          horiz_closest = horiz_candidates(horiz_closest ...
                                           - length (horiz_candidates));
          horiz_edge = 4;
        else
          horiz_closest = horiz_closest - 2*length (horiz_candidates);
          if horiz_closest > length (alt_candidates)
            horiz_closest = alt_candidates (horiz_closest...
                                            - length (alt_candidates));
            horiz_edge = -2;
          else
            horiz_closest = alt_candidates (horiz_closest);
            horiz_edge = -1;
          end
        end
      else
        horiz_closest = horiz_candidates(horiz_closest);
        horiz_edge = 2;
      end

      % Could decide here whether it really is horizontal or vertical.
      % For now, just take what it seems.
      if seems_horiz && ~isempty (horiz_closest)
        closest = horiz_closest(1);
        edge = horiz_edge;
      elseif ~isempty (vert_closest)
        closest = vert_closest(1);
        edge = vert_edge;
      else
        return
      end

      if ismember ('control', modifiers)
        if (~isfield (rectangles(1), 'alt_days'))
              % TODO: reduce redundancy w.r.t. code below
          day = mod1 (round ((coordinates(1) + coordinates_now(1))/2), 7);
          rectangles(closest).alt_days = [day; mod1(day+1,7)];
          rectangles(closest).alt_hrs = rectangles(closest).on_off([2,4]);
        end
        if ~isempty (horiz_closest)
          closest = horiz_closest(1);
        end
        A = 0;
        B = 0;
                % if seems moving a vertical edge
        if ~seems_horiz && ~isempty (rectangles(closest).alt_hrs)
          A = between (rectangles(closest).alt_hrs(1), coordinates(1),...
                       rectangles(closest).on_off(2), 24);
          B = between (rectangles(closest).alt_hrs(2), coordinates(2),...
                       rectangles(closest).on_off(4), 24);
        end
        if A + B > 0
          edge = -3;
        else
          edge = horiz_edge;
          if horiz_edge > 0
            if isempty (rectangles(closest).alt_days)
              % If this rectangle has no alt days, copy from an overlap
              for i = 1:horiz_candidates
                if ~isempty (rectangles(i).alt_days)
                  rectangles(closest).alt_days = rectangles(i).alt_days;
                  rectangles(closest).alt_hrs=rectangles(closest).on_off([2,4]);
                  break
                end
              end
              % If no overlap has alt days either, take from mouse
              % TODO: reduce redundancy w.r.t. code above
              if isempty (rectangles(closest).alt_days)
                day = mod1 (round ((coordinates(1) + coordinates_now(1))/2), 7);
                rectangles(closest).alt_days = [day; mod1(day+1,7)];
                rectangles(closest).alt_hrs = rectangles(closest).on_off([2,4]);
              end
            end
            dist = abs (rectangles(closest).alt_hrs - coordinates(2));
            if length (dist) < 2 || dist(1) < dist(2)
              edge = -1;
            else
              edge = -2;
            end
          end
        end
      end

      % Move the edge by the amount of the mouse move
      add_undo ({'replace', closest, rectangles(closest)});

      if edge > 0       % a non-alt setting
        old_value = rectangles(closest).on_off(edge);
        rectangles(closest).on_off(edge) = ...
          rectangles(closest).on_off(edge) + move(1 + seems_horiz);
      else
        if edge > -3    % changing alt_hrs
          new_time = rectangles(closest).alt_hrs(-edge) + move(1 + seems_horiz);
          [val, pos] = min (abs (new_time - rectangles(closest).on_off([2,4])));
          if val < 1
            new_time = rectangles(closest).on_off(2*pos);
          end
          rectangles(closest).alt_hrs(-edge) = new_time;
                    % remove alt settings if they're not longer needed
          if all (rectangles(closest).alt_hrs(:) ...
                  == rectangles(closest).on_off([2,4]))
            rectangles(closest).alt_hrs = [];
            rectangles(closest).alt_days = [];
            if (isempty ([rectangles(:).alt_days]))
              rectangles = rmfield (rectangles, {'alt_days', 'alt_hrs'});
            end
          end
        else        % changing alt_days
          rectangles(closest).alt_days ...
            = mod1 (rectangles(closest).alt_days + round(move(1)), 7);
        end
      end

      % After a horizontal move (changing start/end date), see if other
      % rectangles abutted this originally, and now overlap.
      % If so, adjust their edges too, unless the user is holding 'SHIFT'.
      if ~seems_horiz && edge > 0 % moving vertical edge -> horizontal move
        if ~ismember('shift',modifiers)
          other_edge = 4 - edge;  % 1 -> 3 and 3 -> 1
          matches = find (abs (on_off (other_edge,:) - old_value) < 2);
          disp (matches);
          for i = matches(:)'
            if i == closest       % avoid issues with very narrow rectangles
              continue
            end
            if ~isempty (intersect (rectangles(closest).burst, ...
                                    rectangles(i).burst))
              oo = rectangles(i).on_off;
              oo(other_edge) = rectangles(closest).on_off(edge);
              if oo(3) > oo(1)    % if we haven't swapped left<->right
                add_undo ({'replace', i, rectangles(i)});
                rectangles(i).on_off(other_edge) ...
                  = rectangles(closest).on_off(edge);
              end
            end
          end
        end
      end

      draw_rectangles;
    end
  end

  function add (~, ~)
    set (f, 'WindowButtonUpFcn', @add_rectangle);
  end

  function [st, en] = align_with_existing (down_pos, up_pos)
    % Find matching start and end times from existing rectangles
    if isempty (rectangles)
      match = [];
    else
      on_off = [rectangles(:).on_off];
      match = find (abs (on_off(3,:) - down_pos(1)) < 2);
    end
    if ~isempty (match)
      st = on_off(3, match(1));
    else
      st = round (down_pos(1));
    end
    if ~isempty (rectangles)
      match = find (abs (on_off(1,:) - up_pos(1)) < 2);
    end
    if ~isempty (match)
      en = on_off(1, match(1));
    else
      en = round (up_pos (1));
    end
  end


  function add_rectangle (~,~)
    coordinates_now = get (ha, 'CurrentPoint');
    coordinates_now = coordinates_now (1, 1:2);

    [st, en] = align_with_existing (coordinates, coordinates_now);

    % even if drawn right-to-left, we need st < en
    % If drawn bottom to top, it is a rectange than spans midnight.
    if st > en
      [st, en] = deal (en, st);
    end

    % 0.5 because rectangle specifies midpoint of pixel.
    on_off = [round(st); coordinates(2)+0.5; round(en); coordinates_now(2)+0.5];

    new_rectangle.on_off = [];           %initialize
    for fn = fieldnames (rectangles(1))'
      new_rectangle.(fn{1}) = [];
    end

    new_rectangle.on_off = on_off;
    new_rectangle.power = mean ([rectangles(:).power]);
    new_rectangle.missed = [];
    if on_off(2) < on_off(4)
      new_rectangle.burst = ceil(on_off(2)):floor(on_off(4));
    else
      new_rectangle.burst = [ceil(on_off(2)):size(cv,1), 1:floor(on_off(4))];
    end
    new_rectangle.burst = round (new_rectangle.burst);
    new_rectangle.trust = mean ([rectangles(:).trust], 2);

    rectangles = [rectangles; new_rectangle];

    add_undo ({'delete', length(rectangles), []});

    % Next time, go back to nudging rectangles
    set (f, 'WindowButtonUpFcn', @final_select_or_move);
    draw_rectangles;
  end

  function toggle_alt_days (object_handle,~)
    if ~isempty (selected_rect)
      if isempty (rectangles (selected_rect).alt_days)
        days_to_hide = [];
      else
        days_to_hide = setdiff (1:7, rectangles (selected_rect).alt_days);
      end
    else
      modifiers = get_modifiers (object_handle);
      if ismember ('shift', modifiers)
        days_to_hide = [];
      else
        days_to_hide = setdiff (1:7, days_to_hide);
      end
    end
    selected_rect = [];
    draw_rectangles;
  end

  function add_undo (action)
    undo_stack{end + 1} = action;
  end

  function undo (~, ~)
    if isempty (undo_stack)
      return
    end
    undo_action = undo_stack{end};
    undo_stack(end) = [];
    r = undo_action{2};
    switch undo_action{1}
      case 'replace'
        rectangles (r) = undo_action{3};
      case 'delete'
        rectangles (r) = [];
      case 'insert'
        rectangles(r+1:end+1) = rectangles(r:end);
        rectangles(r) = undo_action{3};
    end
    draw_rectangles;
  end

  function done (~, ~)
    delete (f);
    f = 0;
  end

  function cancel (~, ~)
    ok = false;
    done;
  end
end
