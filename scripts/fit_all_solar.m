function [cap, az, ze, seen, capFactor, daily_min, vampires, mismatch] = fit_all_solar (data, meta, s)
  global figure1 figure2 figure3 figure10 figure11
  figure1 = figure(1);   plot(1); figure1  = get (figure1,  'CurrentAxes');
  figure2 = figure(2);   plot(1); figure2  = get (figure2,  'CurrentAxes');
  figure3 = figure(3);   plot(1); figure3  = get (figure3,  'CurrentAxes');
  figure10 = figure(10); plot(1); figure10 = get (figure10, 'CurrentAxes');
  figure11 = figure(11); plot(1); figure11 = get (figure11, 'CurrentAxes');

  % (don't use squeeze, as it fails if  size (data,1)==1)
  %max_solar = reshape (-2*min (data,[],2), [size (data,1), meta.SamPerDay]);

  dark = [1:s.dark_end-1, s.dark_start+1:meta.SamPerDay];
  vampires = find_vampires (reshape (permute (data(:,:,dark), [2, 1, 3]), size (data,1)*meta.Days, length(dark)));
  vampires = reshape (vampires, [meta.Days, size(data,1)])';
  %solar_capacity = max (max_solar(:,s.solar_start:s.solar_end),[],2) + vampires;

  smr = [meta.January, meta.February, meta.November, meta.December];
  wtr = [meta.June,meta.July];
  broad = [meta.May, meta.June, meta.July, meta.August, meta.September];

%   % Record maximum in summer, winter and a broadly-defined winter,
%   % along with the day of each.
%   sz = [size(data,1), size(s.SunPos, 2)];
%   sz_year = [size(data,2), size(s.SunPos, 2)];
%   hrs = repmat (1:size (s.SunPos, 2), size (data, 1), 1);
%   solar_range = hrs(1,:) + s.solar_start - 1;

  % Remove vampires to allow seeing generation but no consumption
  data_no_vamp = bsxfun (@minus, data, vampires);

%   % Sun positions on peak generation days in summer and winter
%   [~, spos] = min (data_no_vamp(:, smr, :), [], 2);
%   spos = s.SunPos (sub2ind (sz_year, smr(squeeze (spos(:, :, solar_range))), hrs));
%   spos = reshape (spos, sz);
%
%   [~, wpos] = min (data_no_vamp(:,wtr,:), [], 2);
%   wpos = s.SunPos (sub2ind (sz_year, wtr(squeeze (wpos(:, :, solar_range))), hrs));
%   wpos = reshape (wpos, sz);
%
%
%   %broad_winter = reshape (max (-meta.SamPerDay/24*bsxfun (@minus,data(:,broad,:), vampires(:,broad)),[],2), [size(data,1), meta.SamPerDay]);
%   [~, bpos] = min (data_no_vamp(:,broad,:), [], 2);
%   bpos = s.SunPos (sub2ind (sz_year, broad(squeeze (bpos(:, :, solar_range))), hrs));
%   bpos = reshape (bpos, sz);

  [spos, wpos, bpos] = compute_peak_sunpos(data, data_no_vamp, meta, s);

  % Estimate orientation of solar panel
  % Must estimate three parameters: zenith angle, azimuth and capacity
  fprintf('Fitting solar in summer / winter\n')
  cap = zeros(1,size (data,1));

  daily_min = squeeze (min (data,[],3));
  [~,first_solar_day] = max (daily_min < 0,[],2);
  %[~,last_solar_day]  = max (daily_min (:,end:-1:1) < 0,[],2);
  disconnected = ~(daily_min < 0)';       % ~< instead of >= -> count NaN
  % users for which not all days generate.
  % Precalculate to avoid calling rolling_min on fully-connected users
  not_all = (any (disconnected));
  gap = 11;         % If no net generation for >gap days, assume no solar
          % gap must be odd.
  disconnected(:,not_all) = -rolling_min (-rolling_min (disconnected(:,not_all), gap), gap);
  not_all = (any (disconnected));     % recalc after "noise" removed

  az = cap; ze = cap;

  % 62: New capacity.
  % 106, 1683: ditto

  % 229: Discontinued capacity.

  % 465: winter shade
  % 1920, 1946, 2237, 2245, 2299, 2325, 2348, 2419, 2476: ditto
  % 2625, 2672, 2788, 2874: ditto

  % 2811: summer shade

  % 731: bad fit because of because winter overnight heaters.
  % 1033: ditto

  % 235, 255, 1154: clipped (by inverter capacity?)
  % 1651, 1700, 1755, 1902, 1966, 1982, 2012, 2049, 2056, 2069, 2261?: ditto
  % 2306, 2314, 2564, 2603, 2717: ditto

  % 1464: Meter totally wrong on days 184,185
  %features = zeros (size(data, 1), 11);
%   load features
%   start = find (features(:,1), 1, 'last');
%   if isempty (start)
%     start = 1;
%   end
  start = 1;
  %for i = [1154, 1651, 1755, 1902, 1966, 1982, ...
  %         2012, 2049, 2056, 2069, 2306, 2314, 2564, 2603, 2717]
  stats = zeros (size (data, 1), 8);
  mismatch = zeros (size (data, 1), 1);

  save_file = '../data/saved_orientation.mat';
  if(exist(save_file,'file'))
      orientation_data_cell = load(save_file);
  else
      orientation_data_cell.orientation_data = [];
  end

  jump_to_unsaved = false;
  i = start;
  while(i <= size(data,1))
%   for i = start:size (data,1)
    r1 = smr;
    r2 = wtr;
    % initial guess at: capacity: max generation (assumes sun overhead)
    %               azimuth:  time of max generation

    % Find days near solstices for which solar is connected
    if ~not_all(i)
      wpos = bpos;

      sps = spos(i,:);
      spw = wpos(i,:);
    else
      % Find suitable summer time
          % if seen enough, or summer is all there is, use it
      if sum (disconnected(r1, i)) <= length(smr)/2 ...
          || first_solar_day(i)+length(smr)/2 > meta.Days
        r1 = r1(~disconnected(r1, i));      % only connected days
        sps = spos(i,:);
      else    % Not enough in summer, but maybe some earlier
        range = length(r1);
        d = find(rolling_min (~disconnected(:,i), range));
        if any (d)
          [~, pos] = max (abs(d-wtr(floor (end/2))));
          d = max (d(pos) - floor (range/2), 1);
        else
          filtered = cumsum (~disconnected(:,i));
          filtered = filtered(range+1:end)-filtered(1:end-range);
          [~,d] = max (filtered);
        end
        r1 = d:d+range-1;
        sps = s.SunPos(r1(floor (end/2)),:);
      end

      % Find suitable winter time
      if ~any (disconnected(r2, i))
        spw = wpos(i,:);
      else
        range = length(r2);
        if first_solar_day(i)+range > meta.Days
          % only a few days of solar, so just use summer.
          r2 = r1;
          spw = sps;
        else
          for dummy = 1   % allow "break" in the code below
            if first_solar_day(i) > r2(1)
              r2 = first_solar_day(i):first_solar_day(i)+range;
              if ~any (disconnected(r2, i))
                break;
              end
            end
            start = find(~disconnected(r2,i), 1, 'first');
            if ~isempty(start)
              r2 = r2(start):r2(start)+range-1;
              if ~any (disconnected(r2, i))
                break;
              end
            end
            % Find most wintry run of range connected days
            d = find(rolling_min (~disconnected(:,i), range));
            if any (d)
              [~, pos] = min (abs(d-wtr(floor (end/2))));
              d = max (d(pos) - floor (range/2), 1);
              r2 = d:d+range-1;
            else
              r2 = r1;     % (spw will be set below)
            end
          end
          r2 = r2(r2<= meta.Days);
          if length(r2) < range/2
            r2 = r1;
            spw = sps;
          else
            spw = s.SunPos(r2(floor (end/2)),:);
          end
        end
      end
    end

    % Compute r1, r2, sps and spw
    % [r1, r2, sps, spw] = compute_r1_r2_sps_spw(spos, wpos, bpos, disconnected);

%     % Within r1 and r2, find the three biggest samples,
%     % and the biggest without excluding vampires.
%     %
%     % Find top three in r1 and r2
%     d = sort (data_no_vamp(i, r1, :));
%     d = d(1, 1:min(3, size (d,2)), :);
%     d = d * -meta.SamPerDay/24;
%     nv = min (data(i, r1, :), [], 2) * -meta.SamPerDay/24;
%     a = [d, nv];
%     seen_summer = reshape(a, [size(a,2), 1, size(d,3)]);
%
%     % For visualising
% %     figure;
% %     plot(linspace(0,23,48),[squeeze(d)' squeeze(nv)]);
% %     xlabel('Time of day');
% %     ylabel('Generation - Consumption');
% %     title(sprintf('House %d - Summer',i));
% %     legend({'Largest gen (no vampire)','2nd largest gen (no vampire)','3rd largest gen (no vampire)','Largest gen (with vampire)'},'Location','best');
%
%     d = sort (data_no_vamp(i, r2, :));
%     d = d (1, 1:min (3, size (d,2)), :);
%     d = d * -meta.SamPerDay/24;
%     nv = min (data(i, r2, :), [], 2) * -meta.SamPerDay/24;
%     a = [d, nv];
%     seen_winter = reshape(a, [size(a,2), 1, size(d,3)]);

    [seen_summer, seen_winter] = compute_seen(data(i,:,:), data_no_vamp(i,:,:), r1, r2, meta);

    % Assess smoothness of summer and winter curves
    % and allocate weight to each.
    % (What about inaccurate vampire estimates.)
    seen = [seen_summer, seen_winter];
    %seen2 = -2 * [min(squeeze (data(i, meta.summer, :)))
    %              min(squeeze (data(i, meta.winter, :)))];
    clip_seen = squeeze (max (seen(1,:,:), 0));
    d1seen = diff (clip_seen, 1, 2);
    d1filt = 0.5 * d1seen(:, 2:end-1) ...
             + 0.25 * (d1seen(:, 1:end-2) + d1seen(:, 3:end));
    sumsqd1seen = diag (d1seen * d1seen') ./ (max (squeeze (seen(1,:,:)), [], 2) .^ 2);
    sumsqd1filt = diag (d1filt * d1filt') ./ (max (squeeze (seen(1,:,:)), [], 2) .^ 2);
    [~, pos] = max (squeeze (seen(1,:,:)), [], 2);
    means = sum (bsxfun (@times, clip_seen, 1:size(clip_seen, 2)), 2) ...
            ./ sum (clip_seen, 2);
    stats(i,:) = [sumsqd1seen', sumsqd1filt', pos', means'];

    if(jump_to_unsaved)
%         if(size(orientation_data_cell.orientation_data,1) >= i && ~isequal(orientation_data_cell.orientation_data(i,:),[0 0 0 0]))
        if(length(orientation_data_cell.orientation_data) >= i && ~isempty(orientation_data_cell.orientation_data{i}))
            i = i + 1;
            continue;
        else
            jump_to_unsaved = false;
        end
    end

    f = find(~disconnected(:, i));
    centrality = -182:182;    % 365 days, 0 at middle
    [~, use] = sort (abs (centrality(~disconnected(:, i))));

    max_len = 5;             % number of days to take max over

    len =  max_len * floor (      length (use)/max_len);
    len2 = max_len * floor (0.5 * length (use)/max_len);
    % should check hemisphere
    use_hot_days  = sort (use(floor (end/2)+1:end));
    use_cold_days = sort (use(1:floor (end/2)));
    hot_days  = floor (use_hot_days( max_len:max_len:end) / 5);
    cold_days = floor (use_cold_days(max_len:max_len:end) / 5);
    % New "seen".  with or without vampires, summer or winter, half hour      clear seen
    seen = zeros (2, len / max_len, size (data_no_vamp, 3));
    tmp = squeeze (data_no_vamp(i, f(1:len), :));
    tmp = reshape (tmp(:), max_len, numel (tmp) / max_len);
    seen(1,:) = min (tmp, [], 1);
    tmp = squeeze (data(i, f(1:len), :));
    tmp = reshape (tmp(:), max_len, numel (tmp) / max_len);
    seen(2,:) = min (tmp, [], 1);
    seen = -seen;
    sp = s.SunPos(f(max_len:max_len:end), :);
%    solar_range = bsxfun (@plus, (s.solar_start:s.solar_end)', size (data, 3) * (0:len/max_len-1));
    solar_range = s.solar_start:s.solar_end;
    solar_range = solar_range(:);

    autofit = false;
    if isempty (cold_days)
      i = i + 1;
      continue;
    end
    my_data = squeeze (data(i,:,:));
    my_data = my_data(:, s.dark_end:s.dark_start);

    tic
    [cap(i), az(i), ze(i), ~, mismatch(i), funcs] = findSolarOrientation (seen, sp, meta.SamPerDay, solar_range, my_data, s, hot_days, cold_days);

    if ~autofit
      if ~exist ('solar_by_pc', 'var')
        load solar_by_pc;
      end
      irradiation = ones (size (squeeze (data(i,:,:))'));
      irradiation(s.dark_end:s.dark_start,:) = 2 * ...
      squeeze (solar_by_pc (meta.pclist == meta.postcode(i),:,:))';
      [az(i), ze(i), cap(i), abort, orientation_data, jump_to_unsaved, jump_to_house] = manual_solar_orientation (az(i), ze(i), cap(i), seen, sp, squeeze (data(i,:,:))', squeeze (data_no_vamp(i,:,:))', s, hot_days, cold_days, s.solar_start, s.solar_end, solar_range, meta, meta.UserOffset + i, orientation_data_cell, disconnected(:,i), irradiation, funcs);
      if ~abort
        save(save_file,'orientation_data');
        orientation_data_cell.orientation_data = orientation_data;
      end
    else
      jump_to_house = NaN;
    end
    if(~isnan(jump_to_house))
      i = jump_to_house;
      continue;
    else
      i = i + 1;
    end


%       % Check if output is clipped by the inverter
%       max_seen = max (seen, [], 2);
%       [~, first_80] = max (seen > 0.8 * max_seen, [], 2);
%       [~, last_80]  = max (seen(:, end:-1:1) > 0.8 * max_seen, [], 2);
%       last_80 = size (seen, 2) + 1 - last_80;
%       range_80 = last_80 - first_80;
%       if max (range_80) > 8
%         [~, first_90] = max (seen > 0.9 * max_seen, [], 2);
%         [~, last_90]  = max (seen(:, end:-1:1) > 0.9 * max_seen, [], 2);
%         last_90 = size (seen, 2) + 1 - last_90;
%         range_90 = last_90 - first_90;
%
%         sum_80 = sum (seen > 0.8 * max_seen, 2) ./ (range_80 + 1);
%         sum_90 = sum (seen > 0.9 * max_seen, 2) ./ (range_90 + 1);
%         mx_jump = max (abs (d1seen (:, first_80:min (last_80, size (d1seen, 2)))), [], 2);
%
%         features(i, 2:11) = [range_80', range_90', sum_80', sum_90', mx_jump'];
%         figure (20); imagesc (squeeze (data(i, :, :))');
%         figure (21); plot    (squeeze (data(i, meta.winter, :))');
%         figure (22); plot    (squeeze (data(i, meta.summer, :))');
%         guess = all (features(i,:) >= [0, 12, 7, 11, 0, 0.9, 0.6, 0.4, 0.6, 0.15, 0.2]);
%         features(i, 2:end)
%         fprintf ('%d  my guess: %d.  ', i, guess);
%         g = input ('clipped? (y, n, u) ', 's');
%         if g == 'y'
%           features(i,1) = 1;
%         elseif g == 'n'
%           features(i,1) = -1;
%         elseif g == 'u'
%           features(i,1) = 0;
%         elseif g == 'b'
%           keyboard
%         elseif isempty (g)
%           features(i,1) = 2 * guess - 1;
%         end
%         if mod (i, 100) == 0
%           save features
%         end
%       end

%     if(size (data,1) == 3000) % Hardcoded for now. Only for the solar data
%         save('../data/auto_fit_orientation_data.mat','az','cap','ze');
%     end

%{
    [cap(i), az(i), ze(i), ~, mismatch] = findSolarOrientation ([seen_summer; seen_winter],...
      [sps; spw], meta.SamPerDay, s.solar_start, meta.UserOffset + i, squeeze (data(i, :, :)), r1, r2);

    %if first_solar_day(i) > wtr(1)
    %    fprintf('solar: first %d cap %g  azimuth %g   zenith %g\n', ...
    %        first_solar_day(i), cap(i), az(i), ze(i));
    %end
    if mismatch > 0.25
      % fit using data throughout the year
      % At points spaced one week apart,
      % take minimum over a month either side
      half_range = 30;
      j = 1;
      samples=first_solar_day(i)+half_range:28:meta.Days-half_range;
      if ~isempty (samples)
        seen_u = zeros (length(samples), meta.SamPerDay);
        sp = repmat (s.SunPos(1,1), length (samples), size (s.SunPos, 2));
        for u = samples
          r1 = (u - half_range) : (u + half_range);
          max_solar_part = reshape (max (-meta.SamPerDay/24*data(i,r1,:), [],2), [1, meta.SamPerDay]);
          seen_u(j,:) = max_solar_part - max (max_solar_part([1:s.dark_end-1,s.dark_start+1:end]), [],2);
          sp(j,:) = s.SunPos(u,:);
          j = j + 1;
        end
        [cap(i), az(i), ze(i)] = findSolarOrientation (seen_u, sp, meta.SamPerDay, s.solar_start, meta.UserOffset+i, cap(i), az(i), ze(i));
      end
    end
%}
%    stats (i,:) = [sumsqd1seen', sumsqd1filt', pos', means'];
  end

  az(az>180) = az(az>180) - 360;  % make westerly numbers easier to see

  seen = zeros (size (data,1), meta.Days, (s.dark_start - s.dark_end)+1);
  capFactor = seen;
  for i = 1:meta.Days
    for j = 0:(s.dark_start-s.dark_end)
      % Estimate generation at this time, assuming no cloud
      p1 = cosd(ze);
      p2 = sind(ze).*cosd(s.pp(i,j+1) - az);
      capFactor(:,i,j+1) ...
        = max (0,   bsxfun (@times, s.s1(i,j+1), p1) ...
        + bsxfun (@times, s.s2(i,j+1), p2))';

      % (Actual generation) - (estimate of vampires)
      % Don't correct by meta.SamPerDay/24, as this is to be added to data
      seen(:,i,j+1) = max ((vampires(:,i)-data(:,i,s.dark_end+j)) ./ cap', 0);
    end
  end
  seen(cap < 0.1,:,:) = 0;      % ignore noisy measurements from tiny solar installations

  vampires = [vampires(:,1), diff(vampires,1,2)];
end
