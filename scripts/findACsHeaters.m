function [state] = findACsHeaters (data, cor_solar, cor_pool, mild_v, ind, state, meta, candidates)
%FINDACSHEATERS Identify corrections for some air conditioners and heaters
%   Find which customers have peaks in usage on hot days, and
%   identify the corresponding air conditioner power
%   Find customers that have identifiable large heaters.


  % Find steps
  % use loop to avoid creating huge array
  % Ignore 11pm spikes from hot water
  %meta.peakhours = (4*meta.SamPerDay/24-1):(23*meta.SamPerDay/24);
  ph = length(meta.peakhours);
  % hot_season = meta.summer;   % Alternative
  % hot_season = find (meta.daylight_saving);
  warm_days = meta.max_temp > 24;
  warm_days = rolling_min (-rolling_min (-warm_days, 3));
  hot_season = find (warm_days' | meta.daylight_saving);
  non_hot = setdiff (hot_season, meta.hotDays);
  thresh_ac = zeros (1, size (data, 1));   % can on/off be found by thresholding?
  conf_ac = thresh_ac;                     % certainty device exists

  % Do big heaters in the same loop
  thresh_am = zeros (1, size (data, 1));
  conf_am = thresh_am;
  am_range = floor (size (data,3)/6):ceil (size (data,3)/3);

  thresh_pm = zeros (1, size (data, 1));
  conf_pm = thresh_pm;
  pm_range = floor (size (data,3) / 3):size (data,3);

  pv = cell(size (data, 1), 1);
  cntv = pv;
  sepv = pv;
  score = zeros (size (data, 1), 1);

  corrected = data;
  corrected(:,:) = corrected(:,:) + full (cor_solar) ... % + full (cor_pool) ...
                                  + full (state.cor_hw(ind,:));

  all_cor = zeros (size(data, 3), size(data, 1));
  for i = candidates
    vamp = full(cumsum(state.vampires(ind(i),:)));

    d = squeeze (corrected(i,:,:));

    userData = 2*(squeeze (d)'+ones (size (d,2),1)*state.vacation(i,:));
    userData(1:end-1) = min (userData(1:end-1), userData(2:end));


%      v = userData(1:size(userData, 1)/4.8, meta.winter);
%      [av, bv] = hist (rolling_min (v(:), 5), 20);
%      [pv{i}, cntv{i}, sepv{i}, rectv] = secondPeak (av + 0.1*([av(2:end),0] + [0,av(1:end-1)]),5);
%      score(i) = max(bv(pv{i})) * max (rectv(:));
% %       disp ([bv(pv); cntv]);
% %       if rectv > evalin('base', 'rv')
% %         keyboard
% %       end
% %       df = bv(pv(2)) - bv(pv(1));
% %       if th > df && 0.7 * th < df
% %         th = df;
% %       end

    % TODO 13 Feb 2017:
    % Record conf_ac, shape of cor_ac for each user in batch

    if any (any (isfinite (userData(meta.peakhours, meta.hotDays))))
      [cor_ac, conf_ac(i), thresh_ac(i)] = bigDevices (userData, ...
                                              vamp, ...
                                              meta.hotDays, ...
                                              meta.peakhours, ...
                                              hot_season, ...
                                              1:size (userData,1));
      all_cor(:, i) = max (-cor_ac, [], 2);
    end

%{
    if any (any (isfinite (userData(am_range, meta.winter))))
      [cor_am, thresh_am(i), conf_am(i)] = bigDevices (userData, ...
                                              vamp, ...
                                              meta.winter, ...
                                              am_range, ...
                                              meta.winter, ...
                                              am_range);
    end

    if any (any (isfinite (userData(am_range, meta.winter))))
      [cor_pm_full, thresh_pm(i), conf_pm(i)] = bigDevices (userData, ...
                                              vamp, ...
                                              meta.winter, ...
                                              pm_range, ...
                                              100:meta.winter(end), ...
                                              1:size(data,3));
      heat_days(i) = sum (any (cor_pm_full));
      heat_slots(i) = sum (cor_pm_full(:) ~= 0);
      cor_pm_full = cor_pm_full';
      cor_pm(i,:) = cor_pm_full(:);
    end
%}

    % TODO: What do we do with bins?
    %binsj = bins{1};
    %binsl = bins{2};

    %fprintf('Estimated AC power of %state: %f L %f(%f) J %f(%f)\n', state.NMIs{i}, airconPower_(i,1), airconPower_(i,2), airconPower_2(i,2), airconPower_(i,3), airconPower_2(i,3));
%keyboard

  end         % for

%    [aa, idx] = sort (score);

  % Find average dependency of air con power on time of day
  idx = false (1, size (data, 1));
  for z = 1:size (data, 1)
    %if conf_ac(i) > 0.1
      idx(z) = (any (diff (all_cor(all_cor(:,z) > 0, z))));
    %end
  end
  scaled = bsxfun (@rdivide, all_cor, max (all_cor));
  means = sum (scaled(:, idx), 2) ./ sum (scaled(:, idx) > 0, 2);

  order = find (any (all_cor));
  [~, idx1] = sort (mean (max (data(order,:,:), [], 3), 2) ./ max (0, thresh_ac(order)' - 0.1));
  order = order(idx1);

  next_time = false (1, size (data, 1));
  ac_count = zeros (size (data, 3), size (data, 2));
  ac_sum = ac_count;

  state.mild_summer = intersect (state.mildDays, hot_season);
  state.mild_winter = setdiff (state.mildDays, hot_season);
  global_excess = mean (state.gridLoadAll(meta.hotDays,:)) - mean (state.gridLoadAll(state.mildDays,:));
  gla = state.gridLoadAll';
  gex = bsxfun (@minus, gla, mean (gla(:, state.mildDays), 2));

  remainder = 0 * gla;
  sz = size (data);
  cor_ac = zeros (sz(end:-1:1));
  pl = zeros (1, size (data, 1));
  pla = pl;
  res_on_raw = pl;
  levels = pl;

  % Mask out any cooling reported over winter
  no_cooling = true (size (ac_count));
  no_cooling(:, hot_season) = false;

  hot_times = (gex > 0.5 * max (gex(:)));
  reduced_thresh = false (1, size (data, 1));
  clear_cut = reduced_thresh;

  % Hints:
  % - high power at start then trailing off suggests bringing to a
  %   comfortable temperature and then maintaining it
  % - Extend to at least day 100 in autumn and day 310 in spring
  % - Use temperature more explicitly nearer winter.  (Avoid cold peaks.)
  % - Check for spikes after a cooling run
  % - Check for similar time of day to other runs / non-runs

  % Regress against time-of-day means and temperature
  input_data = permute (corrected, [3 2 1]);
  for pass = 1:2
    for i = candidates
      vamp = full(cumsum(state.vampires(ind(i),:)));
    % for i = [2294, 1250, 1261, 225, 390, 204, 1331, 1998, 962, 613, 1331, 1998, 962, 613, 1803, 1769,  2442, 86, 2214]
      % power more than average for mild day, oriented as (hours, days)
      me_year = input_data(:,:,i);

      me = me_year(:, meta.hotDays);
      all_excess = me_year - 0.5 * repmat (mild_v(i,:)', 1, size (me_year,2));
      excess = all_excess (:, meta.hotDays);
      % Assume only daily max is known.  Regress against it and means
      daily_max = repmat (meta.max_temp(meta.hotDays)', size (me,1), 1);
      all_means = repmat (means(:), 1, size (me, 2));
      if pass == 1 || retry(i)
        ac_times = (me > thresh_ac(i) / 2);
      else
        ac_times = cor_ac(:, meta.hotDays, i) > 0;
      end
      ac_times = ac_times(:);
      if sum (ac_times) <= 1
        cor_ac(:, :, i) = 0;
        figure(10); imagesc (max (0, me_year));
        continue;
      end

      % Don't fit to spurious estimates of airconditioning.
      false_ac = (excess(ac_times) < 0.5 * mean (excess(ac_times)));
      if sum (false_ac) > 0.2 * sum (ac_times)
        aa = find (ac_times);
        ac_times(aa(false_ac)) = false;
      end

      % TODO: regress against
      %  (a) grid power - grid mild-day power
      %  (b) current temperature, half hour ago, 2 hours ago
      %  (In paper, compare the "performance" of different regressions,
      %  possibly measured by simlarity of distribution after correction
      %  to distribution on mild days,
      %  possibly by data sets with ground truth, if big enough,
      %  possibly amount of "correction > observed".)
      % Find new threshold?
      %    Correlate with reconstruction?
      % Find "jumps"
      %    How can we use these to calibrate the regression / thresholds?
      % Clip above with aircon power.
      % Clip below with power of "idling" aircon.

      % TODO:
      % Fix over-estimated confidence: (< means underestimated corr)
      % 2961 (pl right)
      % 113, 562, 1929, 987(est0), 347<, 1203<, 1552<, 2487(miss overnight),
      % 457<, 1145<
      % Mean removal makes very much worse
      % 1221, 1564, (2283), 602
      % Under-confidence
      % 685, 1114, 1625, 2657

      % Fix: 562 (timed device), 2444+1803 (weekends only),
      %      159 (not household), 1164 (solar error?)
      %      2487+2882+1119 (missing overnight),
      %      525+2864+2375+2939+2167+111(!)+355+509(am) (threshold too low),
      %      493+2893+825(!)+1509+1579(!)+2928(!)+1462 (threshold too low),
      %      1216(?)+182+588(?)+1650(?)+730(!) (threshold too low),

      %      2677 (aircon-like on days 100-150)
      %      1125 (periodic signal caught)
      % Reverse cycle:
      % 1878, 1227, 1201, 2487, 2513?, 1780, 1623?, 935, 1152, 2003, 40,
      % 1452, 1264, 2616(heat_only?), 600
      dmm = repmat (meta.max_temp(:)', size (data, 3), 1);
%       spike_clusters;

      recon = fit_power (excess, ac_times, daily_max, all_means, gla, hot_season, meta.hotDays, dmm, thresh_ac(i));

      % TODO: re-fit if any recon < 0 (or even close)
      if any (recon(:) < 0.1 * max (recon(:)))
%         fprintf ('recon < 0.1 * max(recon)\n');
        recon = max (0.1 * max (recon(:)), recon);
%         next_time(i) = true;
%         continue;
      end

      if i == 47
        recon1 = recon;
      end

      myc = me_year(:);
      [~, jump_times, signed_jumps] = join_jumps(myc);
  %    jump_times(signed_jumps > 0) = jump_times(signed_jumps > 0) + 1;    % HACK!!  find_jumps seems to return the wrong value
      jj = zeros (size (myc));
      jj(round (jump_times)) = signed_jumps;
      jt = zeros (size (myc));
      jt(round (jump_times)) = jump_times;

      est = zeros (size (recon));
      est(:, hot_season) = me_year(:, hot_season) > recon(:, hot_season);
      est(est > 0) = recon(est > 0);

%      est = min (est, bsxfun (@minus, me_year, vamp));
%      [pl(i), pla(i)] = mismatch_ratio (est, me_year, state.mild_summer, state.mild_winter);
%      figure(110); plot (sum ((me_year(:, meta.hotDays)-est(:,meta.hotDays)))', meta.max_temp (meta.hotDays), 'o');

      % TODO: Find better criterion and better response to fitting too much
      % TODO: Fit exponential(?) to mild days, and look at probability
      %       of observing data without aircon.
      % TODO: allow for NaNs?  Allow for vacations
      len = sum (~isnan (data(i,hot_season,1)));
      day_ratio = sum (any (est > 0)) / len;
%{
      if day_ratio > 0.5
        % To decide:
        %   Is an aircon present at all?
        %   If so, what power would it be at what time?
        %      and when is it on?
        %
        my_excess = mean (excess, 2, 'omitnan');
        corr_excess = corr(my_excess, global_excess');
        has_aircon = corr_excess > 0.5;

        % If some correlation with temperature,
        % see if few days have long runs or possible aircon use.
%         figure (103); imagesc (me_year);
%         figure (106); imagesc (-est);
%         figure (202); imagesc (max (0, me_year - est))
        if ~has_aircon && corr_excess > 0
          est(:) = -rolling_min(-rolling_min(est(:), 3), 3); % Remove spikes
          figure (203); imagesc (-est);
          day_ratio = sum (any (est > 0)) / len;
          if day_ratio < 0.5 && (corr (sum (est, 2), global_excess') > 0.5 ...
                                 || corr (sum (est(:, meta.hotDays))', ...
                                          state.gridLoad(meta.hotDays)) > 0.3)
            has_aircon = true;
          end
          if corr_excess > 0.3    % Delayed to here to get rolling_min
            has_aircon = true;
          end
        end

        % TODO: Check for peak times in mild days, like morning peak
        % TODO: Check correlation of times with peak use times
        %       e.g., mean(consumption(est ~= 0)) - mean(consumption(est==0))
        %       (No, that misses falsely detected morning peaks)
        %       corr between sum(est        (:, hot_season))
        %                and sum(consumption(:, hot_season))
        if ~has_aircon
          order(order == i) = [];
          figure(10); imagesc (max (0, me_year));
          continue;
        elseif day_ratio > 0.6
          [est, level] = harder_ACs (me_year, thresh_ac(i), all_excess, est, ...
                            day_ratio, corr_excess, gex, ind, state, meta);
          if isempty (est)
            [est, est_on, est_off] = check_overnight (me_year, ...
                                                     zeros (size (me_year)), ...
                                                     [], [], jj, jt, ...
                                                     hot_season, ...
                                                     state, meta);
            if sum (est) == 0
              figure(10); imagesc (max (0, me_year));
              continue;
            end
          else
            levels(i) = level;
          end
        end
      end

%      figure (203); imagesc (max (0, me_year - est))
%}
      new_est = est;
      pr = max (gex / max (gex(:)), 0);
      % Find pre-wakeup rise and post-wakeup dip
      m = mean (me_year(:, hot_season), 2, 'omitnan');
      m(2:end) = 0.5 * (m(1:end-1) + m(2:end));
      d = diff (m);
      r = @(x)(round(x * meta.SamPerDay/48));
      [first_max, rise] = max (d(r(3):r(16)));
      rise = rise + 2;
      peak = m(rise);

      % Check that we didn't miss an earlier rise time
      if rise > r(15) && m(rise) > 1.2 * min(m)
        [pre_max, pre_rise] = max (d(1:rise-1));
        if pre_max > 0.8 * first_max && pre_rise >= r(10)
          rise = pre_rise;
        end
      end

      while rise > r(3) && m(rise) > 0.5 * peak && d(rise - 1) > 0.5 * d(rise)
        rise = rise - 1;
      end
      th = 1.1 * max (m(ceil (rise / 2) : rise - 1));
      if m(rise) < th && m(rise + 1) > m(rise)
        rise = rise + 1;
      end

      [~, fall] = min (d(rise + 1:rise + 8));
      fall = fall + rise;
      while fall < length (d) && d(fall + 1) < 0.3 * d(fall)
        fall = fall + 1;
      end
      fall = min (fall, r(24));
%      figure(11); plot (1:48, m, [rise, rise], [0, max(m)], [fall, fall], [0, max(m)]);

      % daylight saving compensation
      dls = me_year;
      if ~meta.timestamp_uses_daylightsaving
        sam_per_hr = meta.SamPerDay/24;
        %range = 100 * size (me_year, 1) + 2 : 275 * size (me_year, 1) + 2;
        range = find (meta.daylight_saving, 1) * size (me_year, 1) + sam_per_hr ...
                : find (meta.daylight_saving, 1, 'last') * size (me_year, 1) + sam_per_hr;
        dls(range) = dls(range + round (sam_per_hr));
      end
%      figure(10); imagesc (max (0, dls));
      smoothed = -rolling_min (-rolling_min (dls, 3), 3);

      [slopes2, cutoff2, coeff2] = fit_piecewise_exp (smoothed(1:rise-1, state.mildDays));
      [slopes3, cutoff3, coeff3] = fit_piecewise_exp (smoothed(rise:fall, state.mildDays));
      [slopes4, cutoff4, coeff4] = fit_piecewise_exp (smoothed(fall+1:r(34)-1, state.mildDays));
      [slopes5, cutoff5, coeff5] = fit_piecewise_exp (smoothed(r(34):end, state.mildDays));

      A1 = 1:rise-1;
      A2 = rise:fall;
      A3 = fall+1:r(34)-1;
      A4 = r(34):size (new_est, 1);
      new_est(A1, hot_season) = max_likelihood (smoothed(A1, hot_season), recon(A1, hot_season), pr(A1, hot_season), slopes2, cutoff2, coeff2);
      new_est(A2, hot_season) = max_likelihood (smoothed(A2, hot_season), recon(A2, hot_season), pr(A2, hot_season), slopes3, cutoff3, coeff3);
      new_est(A3, hot_season) = max_likelihood (smoothed(A3, hot_season), recon(A3, hot_season), pr(A3, hot_season), slopes4, cutoff4, coeff4);
      new_est(A4, hot_season) = max_likelihood (smoothed(A4, hot_season), recon(A4, hot_season), pr(A4, hot_season), slopes5, cutoff5, coeff5);

      if sum (new_est(:) ~= 0) < 20
        tmp = all_excess;
        tmp(~hot_times) = NaN;
        a = mean (tmp, 2, 'omitnan');
        a(a == max (a)) = NaN;
        a(a == min (a)) = NaN;
        x = mean (a, 'omitnan');
        sd = sqrt (var (a, 'omitnan'));
        if x > 2 * sd && x > 0.2
          recon2 = x * ones (size (me_year, 1), length (hot_season));
          new_est2 = new_est;
          A1 = 1:r(13)-1;
          A2 = r(13):r(20)-1;
          A3 = r(20):r(34)-1;
          A4 = r(34):size (new_est2, 1);
          new_est2(A1, hot_season) = max_likelihood (smoothed(A1, hot_season), recon2(A1, :), pr(A1, hot_season), slopes2, cutoff2, coeff2);
          new_est2(A2, hot_season) = max_likelihood (smoothed(A2, hot_season), recon2(A2, :), pr(A2, hot_season), slopes3, cutoff3, coeff3);
          new_est2(A3, hot_season) = max_likelihood (smoothed(A3, hot_season), recon2(A3, :), pr(A3, hot_season), slopes4, cutoff4, coeff4);
          new_est2(A4, hot_season) = max_likelihood (smoothed(A4, hot_season), recon2(A4, :), pr(A4, hot_season), slopes5, cutoff5, coeff5);
          new_est = max (new_est, min (new_est2, medfilt1 (new_est2, 7)));
          reduced_thresh(i) = true;
        end
      end

      [new_ee, est_on, est_off] = find_edges (me_year, new_est, recon, jj, jt);
      [extra, ~, ~] = check_overnight (me_year, new_ee, ...
                                               est_on, est_off, jj, jt, ...
                                               hot_season, ...
                                               state, meta);
      % Maximum likelihood of check_overnight results
      f = find (extra);
      [rows, cols] = ind2sub (size (extra), f);
      idx = f(rows < rise);
      extra(idx) = max_likelihood (smoothed(idx), extra(idx), pr(idx), slopes2, cutoff2, coeff2);
      idx = f(rows >= rise & rows <= fall);
      extra(idx) = max_likelihood (smoothed(idx), extra(idx), pr(idx), slopes2, cutoff2, coeff2);
      idx = f(rows > fall & rows < r(34));
      extra(idx) = max_likelihood (smoothed(idx), extra(idx), pr(idx), slopes2, cutoff2, coeff2);
      idx = f(rows > r(34));
      extra(idx) = max_likelihood (smoothed(idx), extra(idx), pr(idx), slopes2, cutoff2, coeff2);
      new_est = new_ee + extra;

%      [ee, est_on, est_off] = find_edges (me_year, est, recon, jj, jt);
%      %[est, est_on, est_off] = check_overnight (me_year, ee,
%      [est, ~, ~] = check_overnight (me_year, ee, ...
%                                               est_on, est_off, jj, jt, ...
%                                               hot_season, ...
%                                               state, meta);
      est = new_est;
      est(no_cooling) = 0;
%      new_est(no_cooling) = 0;
%{
      figure(210); image (21 * (1 + (est ~= 0) * 2 - (new_est ~= 0)));
      figure(10); imagesc (max (me_year, 0));
      figure(11); imagesc (max (me_year - new_est, 0));
      figure(12); imagesc (max (me_year - est, 0));
      figure(13); imagesc (max (me_year - min (est, new_est), 0));
      figure(14); imagesc (min (est, new_est));
%}

      % TODO:
      % - Check for false positives
      %   = Neighbours are not less

      smoothed = me_year;
      smoothed(2:end-1) = min (smoothed(2:end-1), max (smoothed(1:end-2), smoothed(3:end)));
      smoothed(no_cooling) = 0;
      ee = est;
      my = smoothed;
      ee(no_cooling) = 0;
      ee = ee(ceil (end/4):end, :);
      my = my(ceil (end/4):end, :);
      idx1 = (ee(:) > 0);

      pwr = mean (ee(idx1));
      if ~isfinite (pwr)
        pwr = 0;
      end


      est = max (0, min (est, bsxfun (@minus, me_year, vamp)));
      [pl(i), pla(i), res_on_raw(i)] = mismatch_ratio (est, me_year, state.mild_summer, state.mild_winter);

      % Too much mismatch between corrected and other times of year
%       if pla(i) + 4 * pl(i) > 7 || (pla(i) + 4 * pl(i) > 6 && pla(i) > 0.1)
       figure(1); imagesc (me_year);
       figure(2); imagesc (max (0, me_year - est));
       if 5 * pla(i) + pl(i) > 6
         next_time(i) = true;
         cor_ac(:, :, i) = NaN;
         continue;
       end

       % Too small to be plausibly an air conditioner, or reliably detected
       if max (max (est)) < 0.1
         cor_ac(:, :, i) = 0;
         continue;
       end

      % Mark as "clear cut"
      % those with a big difference between aircon and background
      s = sum (my(~idx1) > pwr);
      s1 = sum (~idx1 & my(:) > 0.5 * pwr);

      if s < 50 && s1 < 192
        if sum (my(~idx1) > 0.8 * pwr) < 20
          pwr = 0.8 * pwr;
        elseif sum (my(~idx1) > 0.9 * pwr) < 20
          pwr = 0.9 * pwr;
        end
        extra_est = pwr * ((smoothed > pwr) & ~est);
        est = est + extra_est;
        clear_cut(i) = true;
      end

      ac_count(est > 0) = ac_count(est > 0) + 1;
      ac_sum(est > 0) = ac_sum(est > 0) + all_excess(est > 0);

%       figure(111); plot (sum ((me_year(:, meta.hotDays)-est(:,meta.hotDays)))', meta.max_temp (meta.hotDays), 'o');
%       figure (103); imagesc (me_year);
%       figure (104); imagesc (me_year - reshape (est, size (me_year)));
%       figure (105); plot (me_year - reshape (est, size (me_year)));
%       figure (106); imagesc (-est);
%       figure (202); imagesc (max (0, me_year - est))
%       figur e (107); imagesc (ac_count);
%       figure (108); imagesc (ac_sum);
      remainder = remainder + (max (0, me_year - est));
figure(200); imagesc (remainder);
      cor_ac(:, :, i) = est;
fprintf ('i=%d ', i);
    end
    if pass == 1
      cor_ac_first_pass = cor_ac;
      % number of airconditioners on a each time.
      aircon_use = sum (cor_ac(:, :, :) ~= 0, 3);

      % If correlation of a user's aircon with aircon_use is much less than
      % the correlation of the residual with the aircon_use,
      % then "discard" user: subtract the aircon estimate and look for aircon again.
      aircon_cor = corr (reshape (cor_ac, (size (cor_ac, 1) * size (cor_ac, 2)), size (cor_ac, 3)), aircon_use(:));
      residual = input_data;
      residual(isnan (residual)) = 0;
      residual = residual - cor_ac;
      resid_cor =  corr (reshape (residual, (size (cor_ac, 1) * size (cor_ac, 2)), size (cor_ac, 3)), aircon_use(:));

      % Probably discovered a non-aircon device
      discard = aircon_cor < resid_cor / 2 & resid_cor > 0.5 * mean (aircon_cor, 'omitnan');
      input_data(:, :, discard) = input_data(:, :, discard) - cor_ac (:, :, discard);

      % May have discovered aircon, but missed much of it
      retry = aircon_cor < resid_cor & resid_cor > 0.25 * mean (aircon_cor, 'omitnan');
      id = input_data(:, :, retry);
      smoothed = medfilt1 (id, 3, [], 2);
      rs = reshape (smoothed, (size (cor_ac, 1) * size (cor_ac, 2)), sum (retry));
      id(isnan (id)) = 0;
      rs(isnan (rs)) = 0;
      if isempty (rs)
        continue;
      end
      smoothed_cor =  corr (rs, aircon_use(:));
      input_cor = corr (reshape (id, (size (cor_ac, 1) * size (cor_ac, 2)), sum (retry)), aircon_use(:));
      idx = input_cor > smoothed_cor;
      fr = find (retry);
      thresh_ac(fr(idx)) = 0.75 * thresh_ac(fr(idx));
      thresh_ac(fr(~idx)) = mean (rs(:,~idx), 1, 'omitnan') + sqrt (var (rs(:, ~idx), 0, 1, 'omitnan'));

      % if we're discarding, and smoothed seems more like most aircon use,
      % replace input by smoothed input
      fr = fr(~idx);
      fr = fr(discard(fr));
      input_data(:, :, fr) = medfilt1 (input_data(:, :, fr), 3, [], 2);

      % Recaulculate aircon use with discarded users removed
      aircon_use = sum (cor_ac(:, :, ~discard) ~= 0, 3);

      % Replace global excess by mean (global excess, aircon use)
      ratio = max (gex(:)) / max (aircon_use(:));
      gex = 0.5 * (gex + ratio * aircon_use);

      % Recalculate mean dependence of aircon power on time of day
    else
      % Choose results from the pass which seems most like typical aircon
      aircon_cor_2 = corr (reshape (cor_ac, (size (cor_ac, 1) * size (cor_ac, 2)), size (cor_ac, 3)), aircon_use(:));
      idx = (aircon_cor > aircon_cor_2) | all (isnan (aircon_cor_2), 1);
      cor_ac(:, :, idx) = cor_ac_first_pass (:, :, idx);
      cor_ac = min (cor_ac, permute (corrected, [3,2,1]));
    end
  end

  % High use on cold-ish days is more likely to be heating than cooling
  max_temp = meta.max_temp;
  prev_max = [max_temp(1); max_temp(1:end-1)];
  cold_days = (max_temp < 15 ...
            | (max_temp < 20 & prev_max < 25));
  cool_days = (max_temp < 25 & max_temp + prev_max < 45);
  cold_hrs = meta.temperatures' < 18;
  cold_hrs = cold_hrs & bsxfun (@and, cool_days', true (size (cor_ac, 1), 1));
  cold_hrs = cold_hrs | bsxfun (@and, cold_days', true (size (cor_ac, 1), 1));
  old_size = size (cor_ac);
  cor_ac = reshape(cor_ac, numel(cor_ac(:, :, 1)), size (cor_ac, 3));
  cor_ac(cold_hrs,:) = 0;
  cor_ac = reshape(cor_ac, old_size);

  % Remove outliers with "too many" days using aircon for a given number of
  % hours
  curves = zeros (length (hot_season), size (cor_ac, 3));
  for i = 1:size (cor_ac, 3)
    c = cor_ac(:, hot_season, i);
    s = sum (c ~= 0);
    curves(:, i) = sort(s)';
  end

  [d, u] = find (diff (curves) ~= 0);
  vals = curves (sub2ind (size (curves), d+1, u));
  times = zeros (meta.SamPerDay, size (curves, 2));
  deleted = [];
  for i = 1:meta.SamPerDay
    day = (vals == i);
    times (i, u(day)) = d(day);
    times (i, times(i,:) == 0) = NaN;
    [smallest, idx] = min (times(i,:));
    [next, idx1] = min (times(i, 1: size (times, 2) ~= idx));
    while next - smallest > 4
      fprintf ('%d\t', i);
      deleted = [deleted, idx];
      idx = (u == idx);
      d = d(~idx);
      u = u(~idx);
      vals = vals(~idx);
      times(i, idx) = NaN;
      smallest = next;
      idx = idx1;
      [next, idx1] = min (times(i, 1: size (times, 2) ~= idx));
    end
  end
  cor_ac(:, :, deleted) = NaN;

  all_count = squeeze (sum (~isnan (data), 1))';
  state.pl(ind) = pl;

  % Identify times with much airconditioner use
  as = ac_sum(:);
  mm = max (min (ac_sum(:, hot_season)));
  as = max (0, ac_sum - mm);
  as = -rolling_min (-rolling_min (as(:), 9), 9);
  common_times = ac_sum;
  common_times(:) = as;

  % statistics
  % Mean power at each time slot
  capacity = 2 * squeeze (max (max (cor_ac)));
  tot_power  = 2 * sum (cor_ac, 3);
  count = sum (cor_ac ~= 0, 3);
  mean_power = tot_power ./ (count);

  save ([meta.metaDataPath 'aircon_checkpoint1'], 'cor_ac', 'corrected', 'clear_cut', 'pl', 'pla', 'meta');
  %generate_aircon_training (cor_ac, max (corrected, 0), sum (cor_ac, 3), false, false, meta, clear_cut);

  new_remainder = zeros (size (remainder));
  for i = candidates
    me_year = squeeze (corrected(i,:,:))';
    vamp = full(cumsum(state.vampires(ind(i),:)));

    tmp = cor_ac(:, meta.hotDays, i);
    tmp(tmp == 0) = NaN;
    my_means = mean (tmp, 2, 'omitnan');
    idx = isfinite (my_means);
    my_means_mean = mean (my_means(idx));
    all_means_mean = mean (all_means(idx,1), 'omitnan');
    ratio = my_means_mean / all_means_mean;
    if isnan(ratio)
      ratio = 1;
    end;
    my_means(~idx) = medfilt1 (all_means(~idx, 1) * ratio);
    my_means = rolling_min (-rolling_min (-my_means, 5), 5);
    my_means(2:end-1) = (my_means(1:end-2) + my_means(2:end-1) + my_means(3:end))/3;
    my_means = repmat (my_means, 1, size (me, 2));

    ac_times = (cor_ac (:, meta.hotDays, i) ~= 0);

    all_excess = me_year - 0.5 * repmat (mild_v(i,:)', 1, size (me_year,2));
    excess = all_excess (:, meta.hotDays);

    recon = fit_power (excess, ac_times, daily_max, my_means, mean_power, hot_season, meta.hotDays, dmm, thresh_ac(i));

    pr = count / size (data, 3);
    % pr = max (gex / max (gex(:)), 0);
    smoothed = -rolling_min (-rolling_min (me_year, 3), 3);
    r = @(x)(x * meta.SamPerDay / 48);
    [slopes2, cutoff2, coeff2] = fit_piecewise_exp (smoothed(1:r(13)-1, state.mildDays));
    [slopes3, cutoff3, coeff3] = fit_piecewise_exp (smoothed(r(13):r(20)-1, state.mildDays));
    [slopes4, cutoff4, coeff4] = fit_piecewise_exp (smoothed(r(20):r(34)-1, state.mildDays));
    [slopes5, cutoff5, coeff5] = fit_piecewise_exp (smoothed(r(34):end, state.mildDays));

    new_est(1:r(13)-1,    hot_season) = max_likelihood (smoothed(1:r(13)-1,    hot_season), recon(1:r(13)-1,    hot_season), pr(1:r(13)-1,    hot_season), slopes2, cutoff2, coeff2);
    new_est(r(13):r(20)-1,hot_season) = max_likelihood (smoothed(r(13):r(20)-1,hot_season), recon(r(13):r(20)-1,hot_season), pr(r(13):r(20)-1,hot_season), slopes3, cutoff3, coeff3);
    new_est(r(20):r(34)-1,hot_season) = max_likelihood (smoothed(r(20):r(34)-1,hot_season), recon(r(20):r(34)-1,hot_season), pr(r(20):r(34)-1,hot_season), slopes4, cutoff4, coeff4);
    new_est(r(34):end,    hot_season) = max_likelihood (smoothed(r(34):end,    hot_season), recon(r(34):end,    hot_season), pr(r(34):end,    hot_season), slopes5, cutoff5, coeff5);
    [ee, est_on, est_off] = find_edges (me_year, new_est, recon, jj, jt);
    [new_est, ~, ~] = check_overnight (me_year, ee, ...
                                             est_on, est_off, jj, jt, ...
                                             hot_season, ...
                                             state, meta);
    new_est(no_cooling) = 0;
    est = new_est;

    est = max (0, min (est, bsxfun (@minus, me_year, vamp)));
    [pl(i), pla(i)] = mismatch_ratio (est, me_year, state.mild_summer, state.mild_winter);

    ac_count(est > 0) = ac_count(est > 0) + 1;
    ac_sum(est > 0) = ac_sum(est > 0) + all_excess(est > 0);

    new_remainder = new_remainder + (max (0, me_year - est));
    new_cor_ac(:, :, i) = est;
  end

  generate_aircon_training (cor_ac, max (corrected, 0), sum (cor_ac, 3), true, false, meta);

  all_count = squeeze (sum (~isnan (data), 1))';
  state.pl(ind) = pl;

  % Identify times with much airconditioner use
  as = ac_sum(:);
  mm = max (min (ac_sum(:, hot_season)));
  as = max (0, ac_sum - mm);
  as = -rolling_min (-rolling_min (as(:), 9), 9);
  common_times = ac_sum;
  common_times(:) = as;

  % statistics
  % Mean power at each time slot
  capacity = 2 * squeeze (max (max (cor_ac)));
  tot_power  = 2 * sum (cor_ac, 3);
  count = sum (cor_ac ~= 0, 3);
  mean_power = tot_power ./ (count);

  used_capacity = cor_ac;
  for i = 1:size (cor_ac, 3)
    tmp = cor_ac(:, :, i);
    tmp(tmp ~= 0) = capacity(i);
    used_capacity(:, :, i) = tmp;
  end
  mean_capacity = sum (used_capacity, 3) ./ (count);

  temperatures = meta.temperatures(meta.hotDays, :)';
  powers = mean_power(:, meta.hotDays);
  tot_pwr = tot_power(:, meta.hotDays);
  capacities = mean_capacity(:, meta.hotDays);
  cnt = count(:, meta.hotDays);

  figure (1); plot (sort (capacity (capacity > 0)));
  figure (2); imagesc (mean_power);
  figure (3); imagesc (mean_capacity);
  figure (4); plot (temperatures(:), powers(:), '.');
              xlabel ('temperature (C)');
              ylabel ('current power (kW)')
  figure (5); plot (temperatures(:), capacities(:), '.');
              xlabel ('temperature (C)');
              ylabel ('mean capacity (kW)');
  figure (6); plot (temperatures(:), tot_pwr(:), '.');
              xlabel ('temperature (C)');
              ylabel ('total power (kW)');
  figure (7); plot (temperatures(:), cnt(:), '.');
              xlabel ('temperature (C)');
              ylabel ('number of customers');
  figure (8); plot (temperatures(:), cnt(:));
              xlabel ('temperature (C)');
              ylabel ('number of customers');

  figure (11); plot (sort (pl));
               ylabel ('KS ratio');
  gde = 2 * bsxfun (@minus, state.gridLoadAll(meta.hotDays,:), mean (state.gridLoadAll(state.mildDays,:)))';
  figure (14); plot (gde(:), powers(:), '.');
               xlabel ('mean excess power (kW)');
               ylabel ('current power (kW)')
  figure (15); plot (gde(:), capacities(:), '.');
               xlabel ('mean excess power (kW)');
               ylabel ('mean capacity (kW)');
  figure (16); plot (gde(:), tot_pwr(:), '.');
               xlabel ('mean excess power (kW)');
               ylabel ('total power (kW)');
  figure (17); plot (gde(:), cnt(:), '.');
               xlabel ('mean excess power (kW)');
               ylabel ('number of customers');
  figure (18); plot (gde(:), cnt(:));
               xlabel ('mean excess power (kW)');
               ylabel ('number of customers');

  level0 = pl(levels == 0);
  level1 = pl(levels == 1);
  level2 = pl(levels == 2);

  figure (20); plot (sort (level0), (1:length (level0)) / length (level0), ...
                     sort (level1), (1:length (level1)) / length (level1), ...
                     sort (level2), (1:length (level2)) / length (level2));

  figure (21); plot (pl(:), capacity(:), '.');

  remainder_totals = sum (remainder);
  gex_totals = sum (gex);
  corr (remainder_totals(hot_season)', gex_totals(hot_season)')

  remainder_totals = max (remainder);
  gex_totals = max (gex);
  corr (remainder_totals(hot_season)', gex_totals(hot_season)')

  remainder_totals = sum (remainder(1:r(13)-1,:));
  gex_totals = sum (gex(1:r(13)-1,:));
  corr (remainder_totals(hot_season)', gex_totals(hot_season)')

  remainder_totals = sum (remainder(r(13):r(25)-1,:));
  gex_totals = sum (gex(r(13):r(25)-1,:));
  corr (remainder_totals(hot_season)', gex_totals(hot_season)')

  remainder_totals = sum (remainder(r(25):r(37)-1,:));
  gex_totals = sum (gex(r(25):r(37)-1,:));
  corr (remainder_totals(hot_season)', gex_totals(hot_season)')

  remainder_totals = sum (remainder(r(37):end,:));
  gex_totals = sum (gex(r(37):end,:));
  corr (remainder_totals(hot_season)', gex_totals(hot_season)')

  cc = zeros (1, 3000);
  for ii = 1:3000
    rem = squeeze (data(ii,:,:))';
    cc(ii) = corr (rem(:), gex(:));
  end
  [~, id] = sort(-cc);
  display (id(1:10));

  [~, peak] = max (gex(:));
  cor = reshape (cor_ac, [size(cor_ac,1)*size(cor_ac,2), size(cor_ac,3)]);
  [~, peak_day] = max (max (gex));
  on_at_peak = (cor(peak,:) > 0);
  on_in_peak_day = any (cor_ac(:, peak_day,:) > 0);
  on_in_peak_day = on_in_peak_day(:);
  oipd = sum (on_in_peak_day);
  count_on_peak_day = sum (cor_ac(:, peak_day, :) > 0, 3);
  figure (30);
  plot ([1, length(count_on_peak_day)], oipd([1,1]), ...
         1:length(count_on_peak_day), count_on_peak_day);

   max_count_each_day = max (sum (cor_ac > 0, 3));
   tot_count_each_day = sum (any (cor_ac > 0), 3);
   [ss, idx] = sort (tot_count_each_day);
   idx = idx(ss > 0);
   figure (33); plot (1:length (idx), max_count_each_day(idx), ...
                      1:length (idx), tot_count_each_day(idx))

   no_ac = squeeze ((all (all (cor_ac == 0))));
   figure (34); imagesc (squeeze (sum (corrected(no_ac,:, :), 'omitnan'))');

  % Regress against mean aircon power, time-of-day means and temperature
  for i = order(:)'
    % power more than average for mild day, oriented as (hours, days)
    me_year = squeeze (corrected(i,:,:))';

    me = me_year(:, meta.hotDays);
    excess = me - repmat (mild_v(i,:)', 1, size (me,2));
    % Assume only daily max is known.  Regress against it and means
    daily_max = repmat (meta.max_temp(meta.hotDays)', size (me,1), 1);
    all_means = repmat (means(:), 1, size (me, 2));
    ac_times = (me > thresh_ac(i) / 2);
    ac_times = ac_times(:);
    if sum (ac_times) <= 1
      disp (sum (ac_times));
      cor_ac(:, :, i) = 0;
      continue;
    end
    % TODO: regress against
    %  (a) grid power - grid mild-day power
    %  (b) current temperature, half hour ago, 2 hours ago
    %  (In paper, compare the "performance" of different regressions,
    %  possibly measured by simlarity of distribution after correction
    %  to distribution on mild days,
    %  possibly by data sets with ground truth, if big enough,
    %  possibly amount of "correction > observed".)
    % Find new threshold?
    %    Correlate with reconstruction?
    % Find "jumps"
    %    How can we use these to calibrate the regression / thresholds?
    % Clip above with aircon power.
    % Clip below with power of "idling" aircon.
    if sum (ac_times) > 4
      dm = daily_max(ac_times);
      am = all_means(ac_times);
      g = gla(ac_times);
      B = robustfit ([dm, am, g], excess(ac_times));
      dmm = repmat (meta.max_temp(:)', size (data, 3), 1);
      amm = repmat (means(:), 1, size (data, 2));
      recon = B(1) + B(2) * dmm + B(3) * amm + B(4) * gla;
      r = recon(:, hot_season);
      too_small = any (r(:) < 0.1 * max (r(:))) ...
                  || any (max (r) < 0.25 * max (r(:)));
      while B(3) < 0 || B(4) < 0 || too_small
        changed = false;
        if B(3) < 0
          am = zeros (size (am));
          B(3) = 0;
          changed = true;
        end
        if B(4) < 0
          g = zeros (size (g));
          B(4) = 0;
          changed = true;
        end
        if too_small && ~changed  % if fitting to variations makes some too small
          if B(3) + B(4) > 0      % ... then don't fit to the variations.
            am = zeros (size (am));
            g  = zeros (size (g));
          else
            dm = zeros (size (dm));
          end
        end
        B = robustfit ([dm, am, g], excess(ac_times));
        recon = B(1) + B(2) * dmm + B(3) * amm + B(4) * gla;
        r = recon(:, hot_season);
        too_small = any (r(:) < 0.1 * max (r(:))) ...
                    || any (max (r) < 0.25 * max (r(:)));
        if too_small && ~any (B(2:4))
          B(1) = 0.1;
          too_small = false;
        end
      end
    else
      recon = mean (excess(ac_times)) * ones (size (gla));
    end

    est = zeros (size (recon));
    est(:, hot_season) = me_year(:, hot_season) > recon(:, hot_season);

    myc = me_year(:);
    [~, jump_times, signed_jumps] = join_jumps(myc);
%    jump_times(signed_jumps > 0) = jump_times(signed_jumps > 0) + 1;    % HACK!!  find_jumps seems to return the wrong value
    jj = zeros (size (myc));
    jj(round (jump_times)) = signed_jumps;
    jt = zeros (size (myc));
    jt(round (jump_times)) = jump_times;

    [ee, est_on, est_off] = find_edges (me_year, est, recon, jj, jt, common_times);
    est = ee;

    ac_count(est > 0) = ac_count(est > 0) + 1;
    all_excess = me_year - repmat (mild_v(i,:)', 1, size (me_year,2));
    ac_sum(est > 0) = ac_sum(est > 0) + all_excess(est > 0);

%     figure (103); imagesc (me_year);
%     figure (104); imagesc (me_year - reshape (est, size (me_year)));
%     figure (105); plot (me_year - reshape (est, size (me_year)));
%     figure (106); imagesc (-est);
%     figure (202); imagesc (max (0, me_year - est))
%     figure (107); imagesc (ac_count);
%     figure (108); imagesc (ac_sum);
  end

  state.airconPower(ind) = -sum (mean (cor_ac(cor_ac(:) < 0)));
  cor_ac = permute (cor_ac, [3, 2, 1]);
  state.cor_ac(ind,:) = sparse (double (cor_ac(:, :)));
end

function recon = fit_power (power, ac_times, daily_max, all_means, gla, hot_season, hot_days, dmm, thresh_ac)
  g = gla(ac_times);
  if sum (isfinite (power(ac_times) + g)) > 5
    dm = daily_max(ac_times);
    am = all_means(ac_times);
    admm = all_means(:,1) * dmm(1,:);
    adm = admm(:, hot_days);
    adm = adm(ac_times);
    B = robustfit ([dm, am, adm, g], power(ac_times));
    amm = repmat (all_means(:, 1), 1, size (gla, 2));
    recon = B(1) + B(2) * dmm + B(3) * amm + B(4) * admm + B(5) * gla;

    if sum (B(3:5)) >= 0 && max (abs (B)) < 8
      rr = recon(:, hot_days);
      recon = max (recon, min (rr(:)) / 3);
      figure(1); plot (recon);
      figure(2); imagesc (recon);
      return;
    end

    r = recon(:, hot_season);
    too_small = any (r(:) < 0.1 * max (r(:))) ...
                || any (max (r) < 0.25 * max (r(:)));
    while B(3) < 0 || B(4) < 0 || B(5) < 0 || too_small
      changed = false;
      if B(4) < 0
        adm = zeros (size (adm));
        B(4) = 0;
        changed = true;
      end
      if B(3) < 0 && ~changed
        am = zeros (size (am));
        B(3) = 0;
        changed = true;
      end
      if B(5) < 0 && ~changed
        g = zeros (size (g));
        B(5) = 0;
        changed = true;
      end
      if too_small && ~changed  % if fitting to variations makes some too small
        rr = recon(:, hot_days);
        recon = max (recon, min (rr(:)) / 3);
        return;

        if B(3) + B(5) > 0      % ... then don't fit to the variations.
          if B(3) < B(5) && B(3) > 0
            am = zeros (size (am));
          else
            g = zeros (size (g));
          end
        else
          if B(2) < B(4)
            dm  = zeros (size (dm));
          else
            adm = zeros (size (adm));
          end
        end
      end
      B = robustfit ([dm, am, adm, g], power(ac_times));
      % Prevent vastly implausible estimate of constant threshold
      if all (B(2:end) == 0)
        if B(1) < thresh_ac/4
          B(1) = thresh_ac / 2;
        elseif B(1) > thresh_ac
          % TODO: Find a better way to choose between B(1) and thresh_ac(i).
          % TODO: We should distinguish between
          %       the threshold and reconstruction.
          B(1) = (B(1) + thresh_ac) / 2;
        end
      end

      recon = B(1) + B(2) * dmm + B(3) * amm + B(4) * admm + B(5) * gla;
      r = recon(:, hot_season);
      too_small = (any (r(:) < 0.1 * max (r(:))) ...
                  || any (max (r) < 0.25 * max (r(:)))) && max(r(:)) > 0;
    end
  else
    recon = mean (power(ac_times)) * ones (size (gla));
  end
end

function [slopes, cutoffs, coeffs] = fit_piecewise_exp (vec)
  % [SLOPES, CUTOFFS] = FIT_PIECEWISE_EXP (vec)
  op = sort (vec(isfinite (vec)));
  cutoff = ceil (length (op)/5);
  op1 = op(cutoff:end) - op(cutoff);
  slopes(1)  = mean (op1);
  cutoffs(1) = op(cutoff) + log (length (op1) / length (op)) * slopes(1);
  coeffs(1) = exp(cutoffs(1) / slopes(1));

  % Fit second component
  if  log (length (op1)) < op1(end) / slopes(1)
    gap = -(log ((length (op1):-1:1) / length (op1)) + op1(:)' / slopes(1));
    cutoff = find (gap > 0, 1, 'last');
    op2 = op1(cutoff:end) - op1(cutoff);
    slopes(2)  = mean (op2);
    coeffs(2)  = (length (op2) / length (op));
    cutoffs(2) = op1(cutoff) + op(ceil(end/5));
    coeffs(2)  = coeffs(2) * exp (cutoffs(2) / slopes(2));
    % TODO: Recalculate slope1

    final = log (length (op));
    if abs (final +([op2(1), op2(end)]) / slopes(2) + log (length (op2) / length (op))) ...
        > abs (final + ([op1(1), op1(end)]) / slopes(1) + log (length (op1) / length (op)))
      slopes(2) =  slopes(1);
    end
%{
    % Check fit
    figure(fig_num);
    plot (op, log ((length (op):-1:1) / length (op)), ...
          op1([1,end]) + op(ceil(end/5)), -([op1(1), op1(end)]) / slopes(1) + log (length (op1) / length (op)), ...
          op2([1,end]) + op1(cutoff) + op(ceil(end/5)), -([op2(1), op2(end)]) / slopes(2) + log (length (op2) / length (op)), ...
          op1 + cutoffs(1), log(coeffs(1)) - (op1 + cutoffs(1)) / slopes(1), ...
          op2 + cutoffs(2), log(coeffs(2)) - (op2 + cutoffs(2)) / slopes(2));
  else
  figure(fig_num); plot (op, log ((length (op):-1:1)/length (op)), ...
                     op1([1,end])+op(ceil (end/5)), -([op1(1), op1(end)]) / slopes + log (length (op1) / length (op)));
%}
  end
  if all (slopes == slopes(1))
    slopes = slopes(1);
    cutoffs = cutoffs(1);
  end
end

function new_est = max_likelihood (me_year, recon, pr, slopes, cutoff, coeff)
  new_est = zeros (size (me_year));
  log_likelihood_off = -inf (size (me_year));
  log_likelihood_on  = log_likelihood_off;
  me_year_on = me_year - 0.8 * recon;
  for i = 1:length (slopes)
    mask_off = me_year    > cutoff(i); % max (cutoff(i), 0.8 * recon);
    mask_on  = me_year_on > cutoff(i); % max (cutoff(i), 0.8 * recon);
    log_likelihood_on (mask_on)  = log (coeff(i) / slopes(i) * (pr(mask_on)))      - me_year_on(mask_on) / slopes(1);
    log_likelihood_off(mask_off) = log (coeff(i) / slopes(i) * (1 - pr(mask_off))) - me_year(mask_off)   / slopes(1);
  end
  mask = log_likelihood_on > log_likelihood_off;
  new_est(mask) = recon(mask);
end

function [pl_means, pl, resid_on_raw] = mismatch_ratio (est, me_year, similar_baseline, other_baseline)
  % MISMATCH_RATIO - Evaluate statistical fit of cancelling a device

  % Omit any missed days.
  % (Also omit vacations?)
  if any (isnan (me_year(1,:)))
    days = size (me_year, 2);
    valid_days = ~isnan (me_year(1,:));
    est = est(:, valid_days);
    me_year = me_year(:, valid_days);

    tmp = false (1, days);
    tmp(similar_baseline) = true;
    tmp = tmp(valid_days);
    similar_baseline = find (tmp);

    tmp = false (1, days);
    tmp(other_baseline) = true;
    tmp = tmp(valid_days);
    other_baseline = find (tmp);
  end

  % For each time-of-day, find number that have aircon ON (A).
  on_per_slot = min (sum (est ~= 0, 2), length (other_baseline));
  ops = find (on_per_slot);
  A = me_year(est ~= 0) - est(est ~= 0);
  D = max (0, me_year(est ~= 0));

  % Select (B) same number of slots on mild summer days with aircon OFF
  % Select (C) same number of slots on mild non-summer days with airdon OFF
  B = zeros (size (A));
  C = zeros (size (A));
  sources = (est ~= 0)+0;
  j = 0;
  sim = me_year(:, similar_baseline);
  sim(est(:, similar_baseline) ~= 0) = NaN;   % Don't treat "est" as baseline
  for i = ops'
    tmp = find (~isnan (sim(i, :)));
    range = 1:on_per_slot(i);
    if length (tmp) < length (range)
      if ~isempty (tmp)
        tmp = repmat (tmp, 1, ceil (length (range) / length (tmp)));
      else
        continue;     % skip this row.  TODO: get B from another row?
      end
    end
    B(j + range) = sim(i, tmp(range));
    C(j + range) = me_year(i, other_baseline(range));
    sources(i, similar_baseline(tmp(range))) = 2;
    sources(i, other_baseline(range)) = 3;
    j = j + on_per_slot(i);
  end

  if j == 0
    pl = NaN;
    pl_means = NaN;
    resid_on_raw = NaN;
    return;
  end

  figure (200); imagesc (sources);

  % Find (KS?) measure of difference between
  % distribution of residual and days without big device
  [~, ~, residual] = kstest2 (A, B);

  % Find (KS?) measure of difference between
  % distributions of different days within days without big device
  [~, ~, without] = kstest2 (B, C);
  old_without = without;

  pl = residual / (without + 0.05);

  fprintf ('with %8g  without %8g  pl = %8g\n', residual, without, pl);

  if mean (A) > mean (B)
    if mean (A) ~= 0 && mean (B) ~= 0 && mean (C) ~= 0
      [~, ~, residual]    = kstest2 (A / mean (A, 'omitnan'), B / mean (B, 'omitnan'));
      [~, ~, new_without] = kstest2 (B / mean (B, 'omitnan'), C / mean (C, 'omitnan'));
      without = min (without, new_without);
      pl_means = residual / (without + 0.05);
      fprintf ('with* %8g  without* %8g  pl* = %8g\n', residual, without, pl_means);
    else
      pl_means = pl;
    end
  else
    pl_means = pl;
  end

  if pl_means > 1.8
    figure (109); plot (sort (A), (1:length (A)) / length (A), ...
                        sort (B), (1:length (B)) / length (B));
  end

  [~, ~, total] = kstest2 (B, D);
  pl = old_without / total;
  resid_on_raw = residual / total;

  if pl > 1
    figure (109); plot (sort (A), (1:length (A)) / length (A), ...
                        sort (B), (1:length (B)) / length (B), ...
                        sort (C), (1:length (C)) / length (C), ...
                        sort (D), (1:length (D)) / length (D));
    legend ('A', 'B', 'C', 'D');
    figure (203); imagesc (est)
    figure(23); imagesc (me_year - est)
    figure(100); imagesc (me_year)
  end

  % TODO? Check slots just before/after turn on?
  % TODO? Check distributions of hot / mild summer / mild non-summer days

  % TODO? Check time-of-fortnight matches?

  % TODO? Check distribution of durations above quantiles
  % (deciles? quartiles)?
end

function [est, on, off] = find_edges (me_year, est, recon, jj, jt, common_times)
  % FIND_EDGES -- find A/C on/off jumps (not based on threshold).

  est(est ~= 0) = recon(est ~= 0);
  recon = recon(:);
  myc = me_year(:);
  % Find first and last estimated-on times in each burst
  est_on = (est > 0);
  est_off = find (est_on(1:end-1) > est_on(2:end));       % last ON
  est_on  = find (est_on(2:end) > est_on(1:end-1)) + 1;   % first ON
  if isempty (est_off)
    if isempty (est_on)
      on = [];
      off = [];
      return;
    else
      on = est_on;
      off = [];
      est(est_on:end) = recon(est_on:end);
      return;
    end
  elseif isempty (est_on)
    on = [];
    off = est_off;
    est(1:est_off) = recon(1:est_off);
    return;
  end

  % Match them so that each "on" is before the corresponding "off".
  if est_off(1) < est_on(1)
    est_on = [1, est_on];
  end
  if length (est_on) > length (est_off)
    est_off(end+1) = length (myc);
  end
  on = NaN (size (est_on));
  off = on;

  vampires = 0;
  jumps = zeros (length (est_on), 4);   % seen-on, seen-off, recon-on, recon-off

  for i = 1:length (est_on)
    if ~isnan (est_on(i))
      if i == 1
        first = 1;
      else
        first = est_off(i-1)+1;
      end

      % Only look at extension to points with aircon plausibly on.
      range = est_on(i):-1:first;
      small =  find (myc(range) < 0.6 * recon(range), 1);
      if ~isempty (small)
        first = est_on(i) + 1 - small;
      end

      range = first:min (est_on(i) + 1, length(me_year(:)));
      r1 = max (1, range - 1);
      mismatch = abs (jj(r1) - recon(range));
      mismatch(jj(r1) <= 0) = NaN;
      mismatch = mismatch .* sqrt (3 + (length (mismatch) : -1 : 1))';
      [m, pos] = min (mismatch);
      if isnan(m)
        if i > 1 && first == est_off(i-1)+1
          range = est_off(i-1):est_on(i);
          if all (isfinite (range)) && all (myc(range) > 0.5 * recon(range))
            est(range) = min (recon(range), myc(range) - vampires);
            %est(range) = recon(range);
          end
          est_off(i-1) = NaN;
          est_on(i) = NaN;
          off(i-1) = NaN;
          on(i) = NaN;
        elseif i == 1
          est_on(i) = 1;
          on(i) = 1;
        else
          on(i) = est_on(i);
        end
      else
        range = first + pos - 1 : est_on(i);
        %est(range) = min (recon(range), myc(range) - vampires);
        est(range) = recon(range);
        if ~isempty (range)
          on(i) = jt(max (1, range(1) - 1));
          f = floor (on(i));
          c = ceil  (on(i));
          if f ~= c
            est(f + 1) = (f + 1 - on(i)) * recon(c);
          end
        else
          on(i) = est_on(i);
        end
      end
      r = min (max (round (on(i)), 1), length (est(:)));
      jumps(i, [1,3]) = [max(jj(max(r-1, 1)), jj(r)), recon(r)];  % HACK! - find out correct index, r or r-1
    end

%  for i = 1:length (est_off)
    if isnan (est_off(i))
      continue;
    end
    if i == length (est_off)
      last = NaN;       % force "one day" below.
    else
      last = est_on(i+1)-1;
    end

    % Only look at extension to points with aircon plausibly on.
    if isnan (last)       % if i = length(est_off) or est_on(i+1) is NaN
      last = est_off(i) + size (me_year, 1);
    end
    last = min (last, length (myc));

    range = est_off(i):last;
    small =  find (myc(range) < 0.6 * recon(range), 1);
    if ~isempty (small)
      last = est_off(i) - 1 + small;
    end

    range = est_off(i):last;
    mismatch = abs (jj(range) + recon(range));
    %mismatch(jj(max (range - 1, 1)) >= 0) = NaN;
    mismatch(jj(range) >= 0) = NaN;
    mismatch = mismatch .* sqrt (3 + (1:length (mismatch)))';
    [m, pos] = min (mismatch);
    rrange = est_off(i):last;
    if isnan(m) || (abs (jj(range(pos))) < 0.3 * recon(range(pos)) ...
                    && all (myc(rrange) > 0.5 * recon(rrange)))
      if i < length (est_off) && last == est_on(i+1)-1
        rrange = [rrange, range(end)+1];
        if all (isfinite (rrange)) && all (myc(rrange) > 0.5 * recon(rrange))
          est(rrange) = min (recon(rrange), myc(rrange) - vampires);
          %est(rrange) = recon(rrange);
        end
        est_on(i+1) = NaN;
        est_off(i) = NaN;
        on(i+1) = NaN;
        off(i) = NaN;
      elseif i == length (est_off)
        est_off(i) = size (me_year, 1);
        off(i) = size (me_year, 1);
      else
        off(i) = est_off(i);
      end
    else
      if isempty (pos)
        fprintf ('Error: pos is empty\n');
        pos = 1;
      end
      range = est_off(i) : est_off(i) + pos - 1;
      %est(range) = min (recon(range), myc(range) - vampires);
      est(range) = recon(range);
      if ~isempty (range)
        off(i) = jt(range(end));
        f = floor (off(i));
        c = ceil  (off(i));
        if f ~= c
          est(c) = (off(i) - f) * recon(c);
        end
      else
        off(i) = est_off(i);
      end
    end
    r = min (max (round (off(i)), 1), length (est(:)));
    jumps(i, [2,4]) = [min(jj(max(r-1, 1)), jj(r)), recon(r)];

  end

  on_idx  = isfinite (on);
  off_idx = isfinite (off);
  jumps = [jumps(on_idx,1), jumps(off_idx,2), jumps(on_idx,3), jumps(off_idx,4)];
  est_on  = est_on(on_idx);
  est_off = est_off(off_idx);
  on  = on(on_idx);
  off = off(off_idx);
  % TODO: Allow thermostatic on/off of over 0.5*recon after a big on,
  %       if recon is much larger than mean+2sigma of mild days.

  % Adjust based on actual observed jumps for on/off
  a = abs(jumps(:,1:2)) - jumps(:, 3:4);
  [b, c, d] = robust_mean_sd (a);
  if abs (b) > 0.1
  end


  % TODO: if this ON/OFF is short, and not in a common aircon time,
  %       omit it.
  both = isfinite  (est_on) & isfinite (est_off);
  on  = est_on(both);
  off = est_off(both);
  if exist ('common_times', 'var')
    for i = 1:length (on)
      if ~any (common_times(on(i):off(i))) ...
          && off(i) - on(i) < 5
        est(on(i):off(i)) = 0;
        on(i)  = NaN;
        off(i) = NaN;
      end
    end
    on  = on(isfinite  (on));
    off = off(isfinite (off));
  end

end



function [extra, est_on, est_off] = check_overnight (me_year, est, ...
                                                   est_on, est_off, ...
                                                   jj, ~, hot_season, ...
                                                   state, meta)
  extra = zeros (size (me_year));
  rest = me_year - est;    % Only look for overnights that aren't yet counted
  off = mean (rest(   end, state.mild_summer), 'omitnan');
  on  = mean (me_year(end, meta.hotDays), 'omitnan');

  after_hot = meta.hotDays + 1;
  if after_hot(end) > size (est, 2)
    after_hot = after_hot(1:end-1);
  end
  % mild = mean (me_year (:, meta.mildDays), 2);
  % [mn, idx] = min (mild);
  low_time = (1 * meta.SamPerDay/24):(5 * meta.SamPerDay / 24); % TODO: base this on  mild
  hot_morn = max (0, rest (low_time, after_hot));
  [am_hist, bins] = hist (hot_morn(:), 20);
%   peak = secondPeak (pm_hist(:)', 2);
  [peak, ~, ~, rectl] = secondPeak (am_hist(:)', 2);
  peak = sort (peak);
  bp = bins(peak);
  no_hist =  rectl(end, end) < 0.15 * (peak(2) - peak(1)) * am_hist(peak(2));

  if (on < 2*off) && (off < 0.1) && no_hist
    return;
  end
  if no_hist
    thresh = (off + on) / 2;
    thresh2 = (thresh + off) / 2;
  else
    thresh = (bp(2) + bp(1)) / 2;
    thresh2 = (thresh + bp(1)) / 2;
  end

  % Prevent excessively low threshold.
  if thresh < 0.1 * max (est(:))
    re = rest(end, hot_season);
    thresh = mean (re(re > max(est(:)) * 0.25)) * 0.5;
  end

  overnight = (rest(end, hot_season) > thresh);
  morning = any (rest(1:low_time(end), hot_season) > thresh);
  if sum (overnight) > 0.5 * length (hot_season) ...
     || (sum (overnight) > 0.25 * length (hot_season) && on < off + 0.1)
    overnight = false (size (overnight));
  end
  if sum (morning) > 0.5 * length (hot_season) ...
     || (sum (morning) > 0.25 * length (hot_season) && bp(2) < bp(1) + 0.1)
    return;
  end
  baseline = mean (rest (:, [false, ~overnight(1:end-1)]), 'omitnan');

  slots = size (rest, 1);
  for i = find (overnight(:)' | [morning(2:end), false])
    d = hot_season(i);   % actual day
    h = slots * d;
    % Is it already counted?  If so, continue
    if overnight(i) && any (est_on <= h & est_off > h)
      continue;
    end

    % Check for start
    %  - Is est non-zero for this day?
    %  - If so, was the step down much smaller than est?
    if overnight(i)
      prev = find (ceil (est_off / slots) == d, 1, 'last');
    else
      prev = find (ceil (est_off / slots) == d+1, 1);
      hh = find (rest (:, d+1) > thresh, 1);
      if ~isempty (hh)
        h = h + hh;
      end
    end
    if ~isempty (prev)
      step_down = min (est_off(prev) + 1, length (est(:)));
      if me_year(step_down) - me_year(step_down-1) > 0.7 * est(step_down-1) ...
          || (~overnight(i) && h < est_off(prev))
        prev = [];
      end
    end
    if isempty (prev)
      range = h : -1 : h - slots + 1;
    else
      range = h : -1 : est_off(prev);
    end
    if all (rest (range) > thresh)
      start = est_off(prev);
    else
      mismatch = abs (jj(range) + (thresh + me_year(h) - off) / 2);
      %mismatch(jj(max (range - 1, 1)) >= 0) = NaN;
      mismatch(jj(range) >= 0) = NaN;
      mismatch = mismatch .* sqrt (sqrt (4 + (1:length (mismatch))))';
      mismatch = mismatch(1:find(rest(range) < thresh2), 1);
      [m, pos] = min (mismatch);
      if isfinite (m)
        start = range(pos);
      else
        start = h;
      end
    end
    % Check for end
    %  - Is there a step down roughly equal to (thresh - off)?
    %  - Does it stay above thresh until the next ON period?
    last = est_on(min (prev+1, length (est_on))) - 1;
    if isempty (last)
      last = min (h + size (me_year, 1) / 2, length (jj));
    end
    range = h : last;

    mismatch = abs (jj(range) + (thresh + me_year(h) - off) / 2);

    rr = 1:find(rest(range) < thresh2, 1);
    if ~isempty (rr)
      mismatch = mismatch(rr);
      range = range(rr);
    end

    if ~isempty (prev) && start == est_off(prev) ...
        && ~isempty (range) && range(end) == last ...
        && -min (jj(range)) < 0.8 * thresh
      mismatch(:) = NaN;
      mismatch(find(jj(range) <= 0, 1, 'last')) = 0;
    end

    % Extra weight to times after which the aircon is off
    thresholds = (rest(min (range+1, length (rest(:)))) < thresh2) ...
                 & (rest(range) > thresh);
    mismatch(thresholds) = mismatch(thresholds) / 2;

    %mismatch(jj(max (range - 1, 1)) >= 0) = NaN;
    mismatch(jj(range) > 0) = NaN;

%     if me_year(end,d) - off > 3 * thresh ...
%         && all (rest(range) > 0.75 * me_year(end, d)) ...
%         && max (mismatch) < 1.5 * min (mismatch)
%       mismatch(:) = NaN;
%       mismatch(find(jj(range) <= 0, 1, 'last')) = 0;
%     end

    mismatch = mismatch .* sqrt (sqrt (4 + (1:length (mismatch))))';

    [m, pos] = min (mismatch);
    if isnan (m)
      pos = max (length (mismatch) - 1, 1);
    end
    if ~isempty (pos) && pos > slots / 2 ...
                      && (pos > 0.7 * slots || range(1) + pos ~= last)
      continue;
    end

    % How do we determine the power?
    % - take average of
    range = start:range(pos);
    thresh_a = max(est(:)) / 10;
    if ~isempty (range)
      if overnight(i)
        new = min (me_year (range), me_year(end, d)) - off;
      else
            % median not mean to ignore half-on slots at start/end
        new = min (me_year (range), median (me_year (range)) - off);
      end
      if new < thresh_a     % Don't add tiny powers if daytime power is large
        new = 0;
      end
      extra(range) = max (new, est(range));
      if ~isempty (prev)
        est_off(prev) = range(end);
      else
        est_on(end+1) = range(1);
        est_off(end+1) = range(end);
      end
    end
  end
  extra (est ~= 0) = 0;       % Avoid overwriting anything

  est_on  = sort (est_on);
  est_off = sort (est_off);
end
