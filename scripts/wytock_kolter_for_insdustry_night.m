%%{
load ../data;
data(~isfinite (data)) = 0;
data = data(:, :, all (all (data >= 0)));   % Get rid of solar users
data_vec = reshape (data, [size(data,1)*size(data,2), size(data,3)]);
hrs = size (data_vec, 1);

load ../temperatures2013.txt;
temperatures = temperatures2013';
temperatures = temperatures(:);

%%%%%%%%%%%%%%% Lachlan's code for optimisation. %%%%%%%%%%%%%%%%%%%%%%%%%%
% sum all of the days for each house, then sum all houses for each day.
% this command automatically reshapes if there are 2 dimensions of length 1.
daily_totals = sum (sum (data, 1), 3);
% sort days into increasing order and store the indexes in sorted_days.
[sorted_vals, sorted_days] = sort (daily_totals);
% mild day values are the bottom quarter.
mild_days = sorted_days(1:floor(end/4));
% get the average consumption over the mild days for each 48 hour slot.
% this means you have 48 averages for each half hour for each individual
% house. So this calculates each houses individual average consumption for
% each time slot for those days we think they are unlikely to use heating
% and cooling.
typical_use = mean (data(:, mild_days, :), 2);
% data = data minus each houses typical day use for each
% timeslot.
data = max (0, bsxfun (@minus, data, typical_use));

%%%% My initial method put afterwards to try and further reduce noise.%%%%
% sum all 48 hour values.
data_temp = sum(data, 1); 
% combine them to represent one value per day, value being total consumption for that day.
data_daily = reshape (data_temp, [size(data_temp,1)*size(data_temp,2), size(data_temp,3)]); 
% sort daily consumption into increasing order.
sorted_data_daily = sort(data_daily);
% only consider the lower quarter consumption values.
sorted_data_daily([92:365],:) = [];
% create vector data_avg, which is the average half hour consumption on
% "mild days" for each house.
data_avg = sum(sorted_data_daily, 1);
data_avg = data_avg(:);
a = 91 * 48;
data_avg = data_avg / a;
b = data_avg(2, 1); % so this is how you access house 2's mild half hour consumption.
% minus each houses respective mild average from every value of it's data
% vector.
%{
for x = [1: size(data_vec, 2)]
    for y = [1: size(data_vec, 1)]
        data_vec(y, x) = data_vec(y, x) - data_avg(x, 1);
        if data_vec(y, x) < 0
            data_vec(y, x) = 0;
        end
    end
end
%}
data_vec = max (0, bsxfun (@minus, data_vec, data_avg'));

if 0
  cooling = max (0, temperatures - 21);
  heating = max (0, 15 - temperatures);
else
  [H, W] = nnmf (data_vec, 2);
  cooling = double (H(:,2));        % may turn out swapped with heating
  heating = double (H(:,1));        % may turn out swapped with cooling
  heating = heating * norm (cooling, 1) / norm (heating, 1);
  idx1 = cooling < heating * 1;
  idx2 = heating < cooling * 1;
  cooling(idx1) = 0;
  heating(idx2) = 0;
  cooling(cooling < 0.1 * max (cooling)) = 0;
  heating(heating < 0.5 * max (heating)) = 0;
end

figure(11); imagesc(reshape (heating, size (data(:, :, 1))));
figure(12); imagesc(reshape (cooling, size (data(:, :, 1))));
drawnow

% Find Base
total = mean (data, 3, 'omitnan');
[~, mild_days] = sort (mean (total, 1, 'omitnan'));
base = double (mean (total(:, mild_days(1:ceil(end/2))), 2));
base = kron (ones(size (data,2), 1), base);
%}

for i = 2166 % 500 %[1:size(data, 3)]
  d = double (data_vec(:, i));
  d(~isfinite (d)) = 0;
  global to_classify
  d = to_classify(:);
  % See http://web.cvxr.com/cvx/doc/
  cvx_solver sedumi
  cvx_begin
      variable Y(hrs, 4)
      variable theta(hrs, 3)
      variable tmp_c(hrs)
      variable tmp_h(hrs)
%       variable by_week(48, 365)
      
      minimize( norm (tmp_c + [tmp_c(2:end); 0], 1) ...
              + norm (tmp_h + [tmp_h(2:end); 0], 1) ...
              + norm (Y(:,3) - base .* theta(:,3), 1) ...
              + norm (Y(:,4), 1) ...
              + norm (theta(2:end,1) - theta(1:end-1, 1), 1) ...
              + 30 * norm (theta(2:end,2) - theta(1:end-1, 2), 1) ...
              + 0.1 * pow_pos (norm (theta(2:end,3) - theta(1:end-1, 3), 2), 2) ...
              + 0 * norm (Y(2:end,4) - Y(1:end-1, 4), 1) );
       subject to
          sum (Y, 2) == d
          tmp_c == Y(:,1) - cooling .* theta(:,1);
          tmp_h == Y(:,2) - heating .* theta(:,2);
          theta >= 0;
          Y >= 0;
  cvx_end

  for k = 1:4
    figure(k); imagesc (reshape (Y(:,k), size (data(:, :, 1))));
  end
  figure(5); imagesc (to_classify);
  figure(6); imagesc (to_classify - reshape (Y(:, 1), size (data(:, :, 1))));
  figure(7); imagesc (to_classify - reshape (Y(:, 2), size (data(:, :, 1))));
end
