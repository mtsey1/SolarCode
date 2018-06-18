function [m, s, outliers] = robust_mean_sd (data, margin)
  d = data(:);
  if nargin < 2
    margin = 2;
  end
  done = false;
  while ~done
    m = mean (d, 'omitnan');
    s = sqrt (var (d, 'omitnan'));
    omit = (abs (d - m) > margin * s);
    done = all (isnan (d(omit)));
    d(omit) = NaN;
  end
  outliers = reshape (isnan (d), size (data));
end