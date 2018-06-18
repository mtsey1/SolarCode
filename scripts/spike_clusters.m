A = me_year;
layer = zeros ([size(A), 7]);
for len = 3:2:size (layer,3)
  B = -rolling_min (-rolling_min (A, len), len);
  C = A - B;
  layer(:, :, len) = C;
  figure(11); imagesc (B);
  figure(12); imagesc (C);
  A = B;
end

C = layer (:, :, 5);
D = C;
D(D == 0) = NaN;
V = var (D, 0, 2, 'omitnan');
D = mean (D, 2, 'omitnan') + sqrt (V);
D = repmat (D, [1, size(C,2)]);
has = (C > D);
has = -rolling_min (-rolling_min (has, 3), 3);
figure(12);
imagesc (C .* has);

r = ranges(find (has(:)));
len = diff (r) + 1;
d = reshape (diff([0; has(:)]), size (has));
[hr, day] = find (d > 0);
power = C(round ((r(1,:) + r(2,:)) / 2));

all_powers = me_year .* (has);
cs = cumsum (all_powers(:), 'omitnan');
sums = len;
sums(2:end) = cs(r(2,2:end)) - cs(r(1,2:end)-1);
sums(1) = cs(r(2,1));

cs2 = cumsum (all_powers(:) .* all_powers(:), 'omitnan');
sums2 = len;
sums2(2:end) = cs2(r(2,2:end)) - cs2(r(1,2:end)-1);
sums2(1) = cs2(r(2,1));

sd = max (0, sums2 ./ len - (sums ./ len) .^ 2);

X = [hr, day, len(:), power(:), sd(:)];
sp = @spike_dist;
p = pdist (X, @(x, y)(spike_dist (x, y, size (me_year, 1))));
link = linkage (X, 'average', {@(x, y)(spike_dist (x, y, size (me_year, 1)))});
figure(13); dendrogram (link, size (X,1));
T = cluster (link, 'maxclust', 3);
mean_pwr = zeros (3,1);
for i = 1:3
  mean_pwr(i) = mean (X(T == i, 4), 'omitnan');
end
[~, idx] = max (mean_pwr);
mask = zeros (size (me_year));
for i = 1:size (r, 2)
  mask (r(1,i):r(2,i)) = 1;
end
my = me_year;
my(isnan(my)) = 0;
corr (mask(:), my(:))
a = clusterdata (X, 'cutoff', 80, 'distance', @(x, y)(spike_dist (x, y, size (me_year, 1))));

function dist = spike_dist (X, Y, sam_per_day)
  % Columns of X are
  % 1. hour
  % 2. day
  % 3. length
  % 4. power
  % 5. standard deviation
  %
  % spike_dist^2 is
  % 1 * min (abs (hour(X) - hour(Y)), 24 - abs (hour(X) - hour(Y)))
  % + 1/100 * min (abs (day(X) - day(Y)), 365 - abs (day(X) - day(Y)))
  % + 1 * (len (X) - len (Y)) .^ 2
  % + 8 * (power (X) - power(Y)) .^ 2
  % + 8 * (standard_deviation(X) - standard_deviation(Y)) .^ 2
  gap = X - Y;
  tmp = abs (gap (:, 1:2));
  tmp(:, 2) = mod1 (tmp(:, 2), 365);
  gap(:, 1:2) = min (tmp, bsxfun (@minus, [sam_per_day, 365], tmp));
  gap(:,1) = gap(:,1) .^ 1.5;
  gap(:,5) = sqrt (gap(:,5));
  gap = gap .* gap;
  
  dist = sqrt (gap * [1; 0.01; 1; 8; 8]);
end