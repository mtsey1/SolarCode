[~, biggest] = max (max (sum (dat, 3), [], 1));
smth = 0.5 * (dat + medfilt1 (dat, 'omitnan'));
smth = medfilt1 (smth, 'omitnan');

d = diff (squeeze (smth(:, biggest, :)));
[~, mx] = max (d);
indices = mx + (biggest - 1) * size (dat, 1);
indices = indices(:)' + numel (dat(:, :, 1)) * (0:size (dat, 3)-1);
extent = 5 * size (dat, 1)/24;
range = -extent:1+extent;
indices = bsxfun (@plus, indices(:), range);
traces = smth(indices);
traces = bsxfun (@minus, traces, traces(:, extent+1));
traces = bsxfun (@rdivide, traces, traces(:, extent+2));

figure(1); plot (traces');
figure(2); plot (mean (traces, 'omitnan'));
figure(3); plot (range * 24 / size (dat,1), quantile (traces, [0.25, 0.5, 0.75])');
xlabel ('offset (hours)');
ylabel ('normalised power')
axis ([-6, 6, -1, 2]);
print ('-depsc', sprintf ('jump_neighbourhood_%s.eps', meta.dataset_name));

sq = (squeeze (smth(:, biggest, 1:10)));
figure(11); plot (bsxfun (@rdivide, sq, max (sq)))

for k = 1:size (dat, 3)
  figure(1); plot (traces(k,:));
  figure(2); plot ([dat(:, biggest, k), smth(:, biggest, k)]);
end