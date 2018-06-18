function [est, type] = harder_ACs (me_year, thresh_ac, all_excess, est, ...
                              day_ratio, corr_excess, gex, ind, state, meta)
  % Look for right threshold?
  % Omit common times of day
  figure(10); imagesc (me_year);

  % len = sum (~isnan (data(i,meta.summer,1)));
  len = sum (~isnan (me_year(1,meta.summer)));

  ee = est;
  ee(:) = -rolling_min(-rolling_min(ee(:), 3), 3); % Remove spikes
  % TODO: base on number of non-hot days with ee>0
  day_ratio = sum (any (ee > 0)) / len;

  % TODO: reinstate spikes during high solar times

  % TODO: Check for actual jumps near on/off times.

  if day_ratio < 0.6
    figure(20); imagesc (est);
    est = ee;
    figure(21); imagesc (ee);
    figure(22); imagesc (me_year - est);
    type = 1;
    return;
  end

  gex_col = gex(:, meta.summer);
  ee_col = me_year(:, meta.summer);
  glob_corr = corr (ee_col(:), gex_col(:));
  if glob_corr < 0.2
    est = mean_day (me_year, all_excess, gex, meta);
    if isempty (est)
      type = 0;
    else
      type = 8;
    end
    return;
  end

  ee = est;
  ee(me_year < thresh_ac / 2) = 0;
  ee(:) = -rolling_min(-rolling_min(ee(:), 3), 3); % Remove spikes
  day_ratio = sum (any (ee > 0)) / len;

  ee_col = ee(:, meta.summer);
  glob_corr = corr (ee_col(:), gex_col(:));
  if (glob_corr > 0.5 && day_ratio < 0.8) || day_ratio < 0.6
    % TODO: rolling_min, and reinstate spikes during high solar times.
    figure(20); imagesc (est);
    est = ee;
    figure(21); imagesc (ee);
    figure(22); imagesc (me_year - est);
    type = 2;
    return;
  end

  % TODO:
  % Assign belief to each run, based on:
  %  - global excess
  %  - height
  %  - clean on/off times jumps
  %  - Tail of PDF on mild days
  %  - length
  % Delete low-belief runs

  type = 2;
  for ratio = [0.5, 0.7, 0.8, 0.9, 1]
    type = type + 1;
    ee(:, meta.summer) = thresh_ac/2 ...
                       * (all_excess(:, meta.summer) > thresh_ac * ratio);
    ee(:) = -rolling_min(-rolling_min(ee(:), 3), 3); % Remove spikes
    day_ratio = sum (any (ee > 0)) / len;
    ee_col = ee(:, meta.summer);
    glob_corr = corr (ee_col(:), gex_col(:));

    if glob_corr > 0.2 && glob_corr - day_ratio > -0.1
      figure(20); imagesc (est);
      est = ee;
      figure(21); imagesc (ee);
      figure(22); imagesc (me_year - ee);
      return;
    end
  end

  figure(20); imagesc (est);
  figure(21); imagesc (ee);
  figure(22); imagesc (me_year - ee);
  est = mean_day (me_year, all_excess, gex, meta);
  if isempty (est)
    type = 0;
  else
    type = 8;
  end
end

function est = mean_day (me_year, all_excess, gex, meta)
  est = [];
  range = round (0.5 * size (me_year, 1)) : 0.75 * size (me_year, 1);
  me  = median (all_excess(range, meta.hotDays));
  all = median (gex(range, meta.hotDays));
  if corr (me, all) < 0.4
    return;
  end

  afternoons = median (all_excess(range, meta.summer));
  [pm_hist, bins] = hist (afternoons, 20);
%   peak = secondPeak (pm_hist(:)', 2);
  [peak, ~, ~, rectl] = secondPeak (pm_hist(:)', 2);
  peak = sort (peak);
  bp = bins(peak);
  if rectl(end, end) < 0.15 * (peak(2) - peak(1)) * pm_hist(peak(2))
    return;
  end
  thresh = (bp(1) + bp(2)) / 2;
  days = meta.summer(afternoons > thresh);   % Days with aircon
  vals = max (0, all_excess(:, days));
  smoothed = -rolling_min (-rolling_min (vals, 7), 7);

  [all_hist, bins] = hist (smoothed(:)', 20);
  [peak, ~, ~, rectl] = secondPeak (all_hist(:)', 2);
  peak = sort (peak);
  bp = bins(peak);
  if rectl(end, end) < 0.15 * (peak(2) - peak(1)) * all_hist(peak(2))
    return;
  end
  thresh = (bp(1) + bp(2)) / 2;
  found = false;
  for i = 1:5
    smoothed(smoothed <= thresh) = 0;

    if max (smoothed(:)) == 0
      break;
    end

    mm = mean (smoothed(smoothed > 0));
    sd = sqrt (var (smoothed(smoothed > 0)));
    smoothed (smoothed > mm + 2*sd) = mm + 2*sd;

    est = zeros (size (me_year));
    est(:, days) = smoothed;

    a = sum (any (est > 0)) / length (meta.summer);
    b = corr (sum (est(:, meta.summer))', ...
              sum (gex(:, meta.summer), 'omitnan')');
    if i == 1 && b < 0.7 && ...
        (max (est(:)) < min (0.1, 0.5 * mean (me_year (me_year > 0.05))))
      break;
    end

    if (a < 0.5 && b > 0.46) || (b / a > 4.5) || b > 0.9
      found = true;
      if max (est(:)) < 4 * min (est(est > 0))
        break;
      end
    end

    thresh = mean (est (est > 0));
  end

  figure(10); imagesc (max (0, me_year));
  figure(20); imagesc (est);
  figure(21); imagesc (max (0, me_year - est));
  figure(22); hist (afternoons(:)', 20);
  figure(23); hist (smoothed(:)', 20);

  if ~found
    est = [];
  end

end
