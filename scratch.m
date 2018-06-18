load ([evalin('base', 'bulkDataPath'), 'aircon'])

f = zeros (1, size (cor_ac, 3));
recall = f;
precision = f;
pure_precision = f;
pure_recall = f;
pure_f = f;

total_pure_true = 0;
total_pure_fake = 0;
total_pure_miss = 0;
cumsum_mean = 0;
cumsum_mine = 0;

for i = 1:size (cor_ac,3)
  j = cellfun (@(x)(isequal(x, state.NMIs{i})), ground_NMIs);
  j = find (j, 1);
  if isempty (j)
    pure_precision(i) = NaN;
    pure_recall(i) = NaN;
    pure_f(i) = NaN;
    precision(i) = NaN;
    recall(i) = NaN;
    f(i) = NaN;
    continue;
  end
  ground = squeeze (aircon(j, :, :))';
  ground(:, [1:30, 340:end]) = 0;             % omit winter reverse-cycle
  estimate = cor_ac(:, :, i);
  me_year = squeeze (corrected(i, :, :))';

  [est_mean, MAE_mean] = Wytock_mean (me_year, ground);
  MAE_mine = mean (abs (estimate(:) - me_year(:)), 'omitnan');
  
  if isfinite (MAE_mean)
    cumsum_mean = cumsum_mean + MAE_mean;
  end
  if isfinite (MAE_mine)
    cumsum_mine = cumsum_mine + MAE_mine;
  end
  fprintf ('i = %d:  mean %g (%g) mine %g (%g)\n', ...
    i, MAE_mean, cumsum_mean / i, MAE_mine, cumsum_mine / i);

  pure_ground = ground;
  pure_ground(isnan (pure_ground)) = 0;
  pure_est = estimate;
  pure_est(isnan (pure_est)) = 0;
  
  estimate(isnan (ground)) = NaN;
  ground(isnan(estimate)) = NaN;
  if all (isnan (ground(:))) || max (me_year(:)) == 0 || ~(max (max (ground(:, 125:280))) > 0) || sum (any (isnan (me_year))) > 250 || any ((isnan (me_year(:)) & ground(:) > 0) | (ground(:) > 2 * me_year(:) + 0.1))
    pure_precision(i) = NaN;
    pure_recall(i) = NaN;
    pure_f(i) = NaN;
    precision(i) = NaN;
    recall(i) = NaN;
    f(i) = NaN;
    continue;
  end

 %{
  figure(1); plot (ground, estimate, 'o')
  figure(2); hist (ground(ground > 0.1 & estimate < 0.1))
  figure(3); hist (estimate(estimate > 0.1 & ground < 0.1))
  figure(4); imagesc (ground - estimate)
  figure(5); imagesc (squeeze (corrected(i, :, :))')
  figure(6); imagesc (ground);
  figure(7); imagesc (estimate);
 %}

  tp = min (ground, estimate);
  fp = ground - tp;
  fn = estimate - tp;
  precision(i) = sum (tp(:), 'omitnan') / sum (tp(:) + fp(:), 'omitnan');
  recall(i) =  sum (tp(:), 'omitnan') / sum (tp(:) + fn(:), 'omitnan');
  f(i) = 2 / (1 / precision(i) + 1 / recall(i));

  miss = pure_ground .* (pure_ground & ~pure_est);
  fake = pure_est .* (pure_est & (pure_ground < 0.01));
  true = min (pure_ground, pure_est);
  true = sum (true(:), 'omitnan');
  fake = sum (fake(:), 'omitnan');
  miss = sum (miss(:), 'omitnan');
  pure_precision(i) = true / (true + fake);
  pure_recall(i) = true / (true + miss);
  pure_f(i) = 2 / (1 / pure_precision(i) + 1 / pure_recall(i));
  
  if isfinite (true)
    total_pure_true = total_pure_true + true;
  end
  if isfinite (fake)
    total_pure_fake = total_pure_fake + fake;
  end
  if isfinite (miss)
    total_pure_miss = total_pure_miss + miss;
  end
  
  disp([i, precision(i), recall(i), f(i), pure_precision(i), pure_recall(i), pure_f(i)])
  %keyboard
end

fprintf ('mean pred %g mean rec %g ovall prec %g overall rec %g\n', ...
  mean (pure_precision, 'omitnan'), mean (pure_recall, 'omitnan'), ...
  total_pure_true / (total_pure_true + total_pure_fake), ...
  total_pure_true / (total_pure_true + total_pure_miss));

figure(10); plot (precision, recall, 'o');
figure(11); plot (pure_precision, pure_recall, 'o');

%for i = find(pure_precision > 0.1 & pure_precision < 0.5)
[~, idx] = sort (precision + recall);
for i = idx
  j = cellfun (@(x)(isequal(x, state.NMIs{i})), ground_NMIs);
  j = find (j, 1);
  if isempty (j)
    continue;
  end
  ground = squeeze (aircon(j, :, :))';
  estimate = cor_ac(:, :, i);
  estimate(isnan (ground)) = NaN;
  ground(isnan(estimate)) = NaN;
  me_year = squeeze (corrected(i, :, :))';
  
  if all (isnan (ground(:))) || max (me_year(:)) == 0 || ~(max (max (ground(:, 125:280))) > 0) || sum (any (isnan (me_year))) > 250 || any ((isnan (me_year(:)) & ground(:) > 0) | (ground(:) > 2 * me_year(:) + 0.1))
    pure_precision(i) = NaN;
    pure_recall(i) = NaN;
    pure_f(i) = NaN;
    precision(i) = NaN;
    recall(i) = NaN;
    f(i) = NaN;
    continue;
  end
  
  figure(1); plot (ground);
  figure(2); plot (estimate);
  figure(3); plot (ground, estimate, 'o');
  figure(4); imagesc (ground - estimate)
  figure(5); imagesc (me_year)
  figure(6); imagesc (ground);
  figure(7); imagesc (estimate);
  figure(8); imagesc (max (0, me_year - ground));
  figure(9); imagesc (max (0, me_year - estimate));
  
  disp([ground_NMIs(j), precision(i), recall(i), f(i), pure_precision(i), pure_recall(i), pure_f(i)])
  keyboard
end

figure(5); plot (precision, recall, 'o');


function [aircon, MAE] = Wytock_mean (me_year, truth)
  fraction = mean (truth(:), 'omitnan') / mean (me_year(:), 'omitnan');
  aircon = fraction * me_year;
  
  MAE = mean (abs (aircon(:) - truth(:)), 'omitnan');
end