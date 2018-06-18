%%{
load '../../../../../GoogleDrive/FIT1041/heaters/data';
data(~isfinite (data)) = 0;
data = data(:, :, all (all (data >= 0)));   % Get rid of solar users
data_vec = reshape (data, [size(data,1)*size(data,2), size(data,3)]);
hrs = size (data_vec, 1);

load temperatures2013.txt;
temperatures = temperatures2013';
temperatures = temperatures(:);

if 0
  cooling = max (0, temperatures - 21);
  heating = max (0, 15 - temperatures);
else
  [H, W] = nnmf (data_vec, 2);
  cooling = double (H(:,2));        % may turn out swapped with heating
  heating = double (H(:,1));        % may turn out swapped with cooling
  idx1 = cooling < heating * 1;
  idx2 = heating < cooling * 1;
  cooling(idx1) = 0;
  heating(idx2) = 0;
%   cooling = cooling * mean (cooling1) / mean (cooling);
%   heating = heating * mean (heating1) / mean (heating);
%   cooling2 = cooling2 * mean (cooling1(cooling1 ~= 0)) / mean (cooling2(cooling2 ~= 0));
%   heating2 = heating2 * mean (heating1(heating1 ~= 0)) / mean (heating2(heating2 ~= 0));
end

% Find Base
total = mean (data, 3, 'omitnan');
[~, mild_days] = sort (mean (total, 1, 'omitnan'));
base = double (mean (total(:, mild_days(1:ceil(end/2))), 2));
base = kron (ones(size (data,2), 1), base);
%}

%for i = [8, 13, 14, 20, 22, 26, 28, 34, 36, 40, 46, 47, 48, 50, 1:size(data, 3)]
iters = min (size (data, 3), 2000);
all_cooling_theta = zeros (hrs, iters);
all_cooling = all_cooling_theta;
all_heating = all_cooling;
all_heating_theta = all_heating;
figure(1); imagesc (reshape (cooling, size (data (:, :, 1))));
for i = 1:iters
  d = double (data_vec(:, i));
  d(~isfinite (d)) = 0;
  % See http://web.cvxr.com/cvx/doc/
  cvx_solver sedumi
  cvx_begin
      variable Y(hrs, 4)
      variable theta(hrs, 3)
      variable tmp_c(hrs)
      variable tmp_h(hrs)

      minimize( norm (tmp_c + [tmp_c(2:end); 0], 1) ...
              + norm (tmp_h + [tmp_h(2:end); 0], 1) ...
              + norm (Y(:,3) - base .* theta(:,3), 1) ...
              + norm (Y(:,4), 1) ...
              + norm (theta(2:end,1) - theta(1:end-1, 1), 1) ...
              + norm (theta(2:end,2) - theta(1:end-1, 2), 1) ...
              + pow_pos (norm (theta(2:end,3) - theta(1:end-1, 3), 2), 2) ...
              + 0 * norm (Y(2:end,4) - Y(1:end-1, 4), 1) );
       subject to
          sum (Y, 2) == d
          tmp_c == Y(:,1) - cooling .* theta(:,1);
          tmp_h == Y(:,2) - heating .* theta(:,2);
          theta >= 0;
          Y >= 0;
      disp (i)
  cvx_end
  
  all_cooling(:, i) = Y(:,1);
  all_cooling_theta(:, i) = theta(:,1);
  all_heating(:, i) = Y(:,2);
  all_heating_theta(:,i) = theta(:,2);
  
  [U, X, V] = svd (all_cooling, 0);
  figure(3);
  imagesc (reshape (U(:,1), size (data(:, :, 1))));
  disp (i);
  drawnow;

%   for k = 1:4
%     figure(k); imagesc (reshape (Y(:,k), size (data(:, :, 1))));
%   end
%   figure(5); imagesc (data(:, :, i));
end