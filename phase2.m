function state = phase2 (s, meta, data)
% phase 2.
% Estimate use at particular times of day, to get a better solar estimate,
% from which to get better estimates for things like pool pumps.

% Provides:
%    (See variables in "isempty(s)")
% Needs:
%    meta.Days
%    meta.peakhours

%pragma warning on Octave:syntax-extension
%global iif;
global fout_s;

PRINT_THIS_PHASE = 1;

state = s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2                  % Set up data structures
    %%%% data to pass to later phases
    state.airconPower = [];
    state.cooling = [];

    state.heatPower = [];
    state.heating = [];

    state.poolPump = sparse(0);
    state.pools = [];
    state.has_am = [];        % List of NMIs that have a morning load

    state.solarUsed = sparse(1,1);

if PRINT_THIS_PHASE
    fout_s = fopen(meta.outCSV, 'w');
    if fout_s == -1
        error ('could not open %s for writing', meta.outCSV);
    end
    fprintf(fout_s,'NMI,Category (1=Solar 2=HW),Days seen,Max demand (kW), Max demand 3pm-9pm (kW),Size Category,Average Daily Energy Demand (kWh),Hot Day Average Daily Energy Demand (kWh),Cold Day Average Daily Energy Demand (kWh),AC - peak (kW),AC expected hot day cooling energy(kWh),AC -  Hot day Cooling Shape Cluster Number,Heating - Expected Cold Day Heating Energy (kWh),Heating - Cold Day Heating Shape Cluster Number,Days of Vacation, Stay at home,Pool pump size (kW),Average Daily Pool Pump Energy (kWh),Hot Water system Size (kW), Average Daily Hot Water Heating Energy (kWh),Solar Size (kW), Solar Energy Produced (Gross) (kWh/year),Solar Energy Consumed by Customer (kWh/year), Solar orientation (degrees)');
    fprintf(fout_s,',ave%d (kWh)', 1:meta.SamPerDay);
    fprintf(fout_s,',hot%d (kWh)', 1:meta.SamPerDay);
    fprintf(fout_s,',cold%d (kWh)', 1:meta.SamPerDay);
    fprintf(fout_s,'\n');

    state.fout_solar = fopen(meta.out_solarCSV, 'w');
    if state.fout_solar == -1
        error ('could not open %s for writing', meta.out_solarCSV);
    end

    state.fout_pool = fopen(meta.out_poolCSV, 'w');
    if state.fout_pool == -1
        error ('could not open %s for writing', meta.out_poolCSV);
    end

    state.fout_hw = fopen(meta.out_hwCSV, 'w');
    if state.fout_hw == -1
        error ('could not open %s for writing', meta.out_hwCSV);
    end
end     % if PRINT_THIS_PHASE


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif isempty(data)            % Tidy up at end

    state = find_solar_by_pc (s, meta);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else                            % Main process, called once per pass
    ind = meta.UserOffset + (1:size(data,1));
    candidates = 1:size (data, 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate solar correction based on phase1's tentative estimates
    [su, cor_solar] = calc_cor_solar (data, s, meta, ind);
    state.cor_solar(ind, :) = cor_solar;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate pool pumps
    fprintf('Detecting pool pumps...\n');
    pump = zeros(1,size(data,1));
    pump_pattern_class = zeros (1, size(data,1));
    %records = zeros(length(ind),4);
%big_spikes = [ 554 6454 5680 4378 4789 4537 5173 7791 7379 2775 7293 7501 6932 6696 276 486 2488 985 8850 8899 5701 173 2988 762 4505 2069 7510 9137 5058 5515 1658 3412 2798 1603 8384 2952 361 3935 7063 1072 1667 3452 8989 9103 3948 7629 508 7303 3875 2331 2172 3465 6443 1679 3640 4824 1000];
%num_NMIs = cellfun (@str2double, s.NMIs(ind));
%for j = 1:length (big_spikes)
%  i = find (num_NMIs == big_spikes(j));
%  if ~isempty (i)
%    corrected = squeeze(data(i,:,:));
%    corrected(:) = corrected(:)' + s.cor_hw(ind(i),:) + cor_solar(i,:);
%    figure(1); imagesc (corrected');
%    figure(2); plot (corrected'(:,200:210));
%    figure(3); plot (corrected'(:,200:210)(:));
%    keyboard
%  end
%end
    has_am = false (1, size (data,1));
    cor_pool = zeros (size(data(:,:), 2), ceil (0.04 * size (data,1)));
    pool_count = 0;
    has_pool = zeros (ceil (0.04 * size (data,1)), 1);

    regress_test = false || false;
    %%% Regression testing framework.
    %%% To compare before/after a patch, set a breakpoint here and
    %%% set regress_test to true.
    if regress_test
      base      = load ('rectangles_base.mat');
      base_pump = base.pump;
      base      = base.rectangles;
      patched   = load ('rectangles_patched.mat');
      patched_pump = patched.pump;
      patched   = patched.rectangles;
      if length (base_pump) ~= length (patched_pump)
        fprintf ('Mismatched lengths.  Check manually\n');
      end

      fprintf ('%d had different pump status\n', ...
               sum (~base_pump ~= ~patched_pump));
      ranking = -2 * (base_pump | patched_pump) + (base_pump & patched_pump);
      [~, idx] = sort (ranking);
      for i = idx(:)'
        if isempty (base{i}) && isempty (patched{i})
          continue;
        end

        ad = [isfield(base{i}, 'alt_days'), isfield(patched{i}, 'alt_days')];
        if length (base{i}) ~= length (patched{i}) ...
           || any (any ([base{i}.on_off] ~= [patched{i}.on_off])) ...
           || ad(1) ~= ad(2) ...
           || (all(ad) && any(~isequal ([base{i}.alt_days], [patched{i}.alt_days])))

          fprintf ('%d: ranking %d\n', i, ranking(i));
          corrected = double (squeeze(data(i,:,:)));
          corrected(:) = corrected(:)' + full(s.cor_hw(ind(i),:)+cor_solar(i,:));

          figure(100); show_rectangles (base{i},    corrected', []);
          figure(101); show_rectangles (patched{i}, corrected', []);
          if length (base{i}) == length (patched{i})
            bb = base{i};
            pp = patched{i};

            for k = 1:min (length (base{i}), length (patched{i}))
              printed_rectangle = false;
              if any (bb(k).on_off ~= pp(k).on_off)
                fprintf ('%2d: %f %f %f %f\n -> %f %f %f %f\n', ...
                         k, bb(k).on_off, pp(k).on_off);
                printed_rectangle = true;
              end
              if all(ad)
                if ~isequal (bb(k).alt_days, pp(k).alt_days)
                  if ~printed_rectangle
                    fprintf ('%2d: %f %f %f %f\n', k, bb(k).on_off);
                  end
                  bbad = [bb(k).alt_days; 0; 0];
                  ppad = [pp(k).alt_days; 0; 0];
                  fprintf ('%2d: Alt Days %d %d -> %d %d\n', ...
                           k, bbad(1:2), ppad(1:2));
                end
              elseif ad(1)
                if ~isempty (bb(k).alt_days)
                  if ~printed_rectangle
                    fprintf ('%2d: %f %f %f %f\n', k, bb(k).on_off);
                  end
                  fprintf ('%2d: Alt Days %d %d ->\n', ...
                           k, bb(k).alt_days);
                end
              elseif ad(2)
                if ~isempty (pp(k).alt_days)
                  if ~printed_rectangle
                    fprintf ('%2d: %f %f %f %f\n', k, bb(k).on_off);
                  end
                  fprintf ('%2d: Alt Days -> %d %d\n', ...
                           k, pp(k).alt_days);
                end
              end
            end
          end
          keyboard;
        end
      end
    end

    rectangles = cell(1, size (data, 1));
    % Examples to fix:
    % 1  74 109 156(extended too far) 232(DoW) *238*(DoW,twice per day)
    %for i = [1, 74, 109, 111, 156, 232, 238]
%    for i = [74, 84, 238]

%if exist([meta.metaDataPath 'cache_cor_pool.mat'])
%  load ([meta.metaDataPath 'cache_cor_pool.mat'])
%else
if true
  cor_pool = zeros (size (data));
  state.poolPump = false (size (data, 1));
  else
    for i = candidates
%{
%    % Solar entities that seem to have pools
    for i = [6, 1, 51, 74, 84, 109, 111, 123, 148, 154, 156, 168, 232, 238, 254, 263, ...
             313, 314, 349, 352, 353, 359, 361, 373, 375, 384, 388, 402, 457, 463, 464, ...
             530, 533, 548, 557, 559, 568, 569, 572, 581, 584, 595, 599, 628, 645, ...
             709, 716, 728, 735, 752, 754, 767, 774, 794, 820, 822, 833, 836, 838, 839, 847, 848, 855, 857, 861, 865, 868, 872, 876, 883, 890, 892, ...
             904, 905, 913, 949, 973, 975, 977, 1007, 1027, 1044, 1046, 1047, 1055, 1058, 1062, 1064, 1066, 1068, 1075, 1079, 1086, ...
             1107, 1124, 1133, 1150, 1159, 1161, 1167, 1224, 1226, 1256, 1271, 1280, 1285, ...
             1335, 1375, 1407, 1413, 1432, ...
             1525, 1577, 1597, 1669, 1708, 1716, 1718, 1719, 1751, 1790, ...
             1816, 1848, 1849, 1851, 1856, 1864, 1870, 1892, 1898, 1899, ...
             1945, 1964, 2025, 2080, 2087, 2090, ...
             2137, 2155, 2156, 2160, 2167, 2170, 2175, 2179, 2195, 2196, ...
             2200, 2204, 2228, 2235, 2249, 2251, 2259, 2260, 2261, 2263, 2271, 2273, 2276, ...
             2307, 2325, 2330, 2358, 2360, 2366, 2371, 2379, 2397, ...
             2400, 2407, 2412, 2427, 2437, 2452, 2487, 2495, ...
             2511, 2515, 2517, 2528, 2578, 2598, ...
             2609, 2612, 2619, 2646, 2655, 2662, 2676, 2694, ...
             2707, 2710, 2713, 2731, 2734, 2736, 2741, 2770, 2773, 2777, ...
             2808, 2825, 2835, 2857, 2861, 2887, 2891, ...
             2904, 2987, 2990, 2999]
%}
%{
    for i = [774, 794, 820, 822, 833, 836, 838, 839, 847, 848, 855, 857, 861, 865, 868, 872, 876, 883, 890, 892, ...
             904, 905, 913, 949, 973, 975, 977, 1007, 1027, 1044, 1046, 1047, 1055, 1058, 1062, 1064, 1066, 1068, 1075, 1079, 1086, ...
             1107, 1124, 1133, 1150, 1159, 1161, 1167, 1224, 1226, 1256, 1271, 1280, 1285];
%}
    % To do:
    %        Find out why overlapping rectangles are not cleared
    %        Better test for on/off times for "new" rectangles
    %        Use fractional on/off times in set_cor_pool_precise
    %        Recalculate power
    %        Look for new rectangles in "clean" areas
    %        Tweak start/end and on/off times
    %        For long bursts: look for brief changes
    %                         look again for DoW changes
%{
    for i = [1027, 1044, 1046, 1047, 1055, 1058, 1062, 1064, 1066, 1068, 1075, 1079, 1086, ...
             1107, 1124, 1133, 1150, 1159, 1161, 1167, 1224, 1226, 1256, 1271, 1280, 1285];
%}
%{
    % Have day-of-week, or are mistaken for day-of-week
    for i = [6, 238, 628, 767, 861, 949, 1079, 1086, 1271, 822, 892, 1044 1159, 1280]
%}

%{
    % Constant settings interfere with solar estimation
    for i = [1055, 1124, 1159, 352, 569, 794, 905]
%}

% HW: 342, 574, 650, 891, 1190, 1258, 1264?, 1286
%{
    % missing rectangles / days
    for i = [402, 464, 568, 123, 238, 263, 833, 890, 1075, 1079, 1133, 1161 ...
             373, 533, 584, 595, 628, 820, 822, 892, 1047, 1086]
%}
%{
     % morning taken as edge
     for i = [111, 232, 709, 752, 872, 905, 913, 1044]
%}
%{
     % other rectangle taken as edge
     for i = [6]
%}
%{
     % other duration too long / short
     for i = [109, 313, 457, 735, 863, 1161, ...
              794, 1256]
%}
%{
     % outlier in fictitious gap
     for i = [74, 599, 973, 1027]
%}
%{
     % over-correct for under-estimated solar
     for i = [349]
%}

% TODO:
% - display heat map of scalar_quality
% - Estimate power fron "best" edge?
%   (i.e., one with best Kolmogorov-Smirnoff match to current edges)
% - Estimate "morning" using daylight saving,
%   and penalise it in edge quality
% - Avoid outliers in fictitious gaps.  (Avoid creating, then prune.)
        corrected = double (squeeze(data(i,:,:)));
        corrected(:) = corrected(:)' + full(s.cor_hw(ind(i),:)+cor_solar(i,:));

        vamp = full(cumsum(s.vampires(ind(i),:)));
%        try
 %fprintf ('%d ', i);
          [p, cor, power, ~, ~, ~, ~, rect_am, rect] = has_pool_pump(corrected, vamp, meta, ...
                                          s.NMIs{ind(i)}, cor_solar(i,:));
drawnow;
          if ~isempty (rect_am)
            has_am(i) = true;
          end
          % Regression testing framework
          rectangles{i} = rect;

%if any (big_spikes == str2double (s.NMIs{ind(i)}))
% figure(1); imagesc (corrected');
% keyboard
%end

          if (~p)
            [p, cor, power] = has_untimed_pump (corrected + cor, vamp, ...
                                          meta, ind(i), any(su == i), state);
            cor = cor(1,:);
            if (p)
              % figure(1); imagesc ((corrected)');
              %figure(2); plot ((corrected)');
%              figure(3); imagesc (reshape (cor, [365, 48])');
              kb_later = true;
            end
          else
              kb_later = false;
          end

          pump(i) = power(1);
          pump_pattern_class(i) = p;
          if p>0
              pool_count = pool_count + 1;
              cor_pool(1:length (cor(:)), pool_count) = cor(:);
              has_pool(pool_count) = i;
              fprintf('%d:%d:%g\n', ind(i), p, power);
              if kb_later
%                keyboard
              end
          end

%        catch
%           vv = cumsum(s.vampires(ind(i),:));
%           save (['bug_',state.NMIs{ind(i)}], 'corrected', 'vv', 'meta', 'ind', 'i', 'su');
%           fprintf ('Error detecting pool pump for user %s\n', state.NMIs{ind(i)});
%           s = lasterror;
%           fprintf ('%s in %s line %s\n', s.identifier, s.stack(1).name, s.stack(1).line);
%           keyboard
%        end
    end

    %Regression testing framework
    %save rectangles_base.mat rectangles pump;
    save ([meta.metaDataPath 'rectangles_patched.mat'], 'rectangles', 'pump');

%    [~, idx] = sort (-pump);
    state.has_am = [state.has_am, s.NMIs(ind(has_am))];
    state.poolPump(ind) = pump;
    state.poolPumpPattern(ind) = pump_pattern_class;
    days = size(cor_pool, 1);
    has_pool = has_pool(1:pool_count);
    cor_pool = cor_pool(:,1:pool_count);
    cor_pool = sparse (kron(has_pool, ones(days, 1)), repmat (1:days, 1, pool_count), cor_pool(:));
    if size (cor_pool, 1) < size (data, 1)
      cor_pool(size (data, 1), size (cor_solar, 2)) = 0;
    end
    fprintf('\n');
    save ([meta.metaDataPath, 'cache_cor_pool.mat'], 'cor_pool');
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Estimate airconditioner peak power, and power used on cooling
    if meta.do_cooling
      fprintf('cooling... ');

      % Power consumption on hot days ("cooling") and hot non-vacations ("cooling_v")
      [i2, i1] = meshgrid(1:meta.SamPerDay, meta.hotDays);
      s2i = sub2ind ([meta.Days meta.SamPerDay], i1(:), i2(:));
      hot_no_hw_solar = reshape(data(:,s2i) + full(s.cor_hw(ind,s2i) + cor_solar(:,s2i)), [size(data,1), length(s2i)/meta.SamPerDay, meta.SamPerDay]);
      cooling    = squeezeNaN(@mean, hot_no_hw_solar, 2);
      cooling_v  = squeezeNaN(@mean, hot_no_hw_solar + repmat(full(s.vacation(ind,meta.hotDays)),[1,1,meta.SamPerDay]), 2);
      clear hot_no_hw_solar;

      % Subtract power consumption on mild days ("cooling") and mild non-vacations ("cooling_v"), to estimate the actual cooling power
      [i2m, i1m] = meshgrid(1:meta.SamPerDay, s.mildDays);
      s2im = sub2ind ([meta.Days meta.SamPerDay], i1m(:), i2m(:));
      mild_no_hw_solar = reshape(data(:,s2im) + full(s.cor_hw(ind,s2im) + cor_solar(:,s2im)), [size(data,1), length(s2im)/meta.SamPerDay, meta.SamPerDay]);
      mild   = squeezeNaN(@mean, mild_no_hw_solar, 2);
      mild_v = squeezeNaN(@mean, mild_no_hw_solar + repmat(full(s.vacation(ind,s.mildDays)),[1,1,meta.SamPerDay]), 2);
      clear mild_no_hw_solar;

      %state.hot(ind,:) = cooling;
      cooling   = cooling   - mild;
      cooling_v = cooling_v - mild_v;

      airconPower_(:,1) = meta.SamPerDay/24*squeezeNaN(@max,cooling_v,2);
      % Avoid undersized array if for loop skips users
      airconPower_2 = zeros(size(airconPower_));
      state.mildDay(ind,:) = mild;

      fprintf ('Found difference profile.  Finding on times... ');
      state = findACsHeaters (data, cor_solar, cor_pool, mild_v, ind, state, meta, candidates);
      fprintf('Done\n');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Should estimate heating for each hour of the day, to correct winter
    % solar estimates

if PRINT_THIS_PHASE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Recalculate solar, correcting for periodic loads

    solar_user_range = find (s.solar_users >  meta.UserOffset & ...
                             s.solar_users <= meta.UserOffset + size (data, 1));
    if ~isempty(solar_user_range)
      solar_users = s.solar_users (solar_user_range);
      corrected = squeeze(data(su,:,:));
      %corrected(:) = corrected(:,:) + full(s.cor_hw(ind(su),:)+cor_pool(su,:));
      [cap, az, ze, seen, capFactor, daily_min, vampires] ...
                          = fit_all_solar (corrected, meta, s);

      state.solar_cap(solar_user_range) = cap;
      state.solar_az (solar_user_range) = az;
      state.solar_ze (solar_user_range) = ze;
      state.seen     (solar_user_range,:,:) = seen;
      state.capFactor(solar_user_range,:,:) = capFactor;
      state.daily_min(solar_user_range,:) = daily_min;
      state.vampires (solar_users+meta.UserOffset,:) = vampires;
    end
else
%%% In phase 3 if doing re-estimation of solar etc
    %summary.ACpower (ind)= state.airconPower(ind);

    %    %%% Partial cluster these users, to 10 clusters
    %    normHourly = bsxfun(@rdivide, cooling, squeezeNaN(@mean, cooling,2));
    %    n = 10;
    %    [classes, centroids, sumd] = kmeans (normHourly, n, 'EmptyAction','singleton', 'Replicates', 100);
    %    counts = zeros(1,n);
    %    for i = 1:n
    %        counts(i) = sum(classes == i);
    %    end
    %    state.coolCentroids(ind,:) = centroids;
    %    state.coolCounts(ind,:) = counts;
    if meta.do_cooling
      state.cooling(ind) = max(sum(cooling,2),0);
      cooling = max(cooling, 0);
      normHourly = bsxfun(@rdivide, cooling, squeezeNaN(@max, cooling,2)+1e-6);
      n = 6;
      if ~isfield(state, 'coolCentroids')  % if first pass, generate clusters
	  [classes, centroids] = kmeans (normHourly, n, 'EmptyAction','singleton', 'Replicates', 100);
	  counts = zeros(1, n);
	  for i = 1:n
	      counts(i) = sum(classes == i);
	  end
	  [~, sorted] = sort (-counts);
	  state.coolCentroids = centroids(sorted,:);
	  fout = fopen (meta.clusterCSV, 'w');
	  if ~fout
	      error ('Error opening %s for writing\n', meta.clusterCSV);
	  end
	  fprintf(fout, 'Cooling\n');
	  for i = 1:n
	      fprintf(fout, '%g,', centroids(i,1:end-1));
	      fprintf(fout, '%g\n', centroids(i,end));
	  end
	  fclose(fout);
      end
      dists = zeros(length(ind), n);
      for i = 1:n
	  dists(:,i) = sum (bsxfun(@minus, normHourly, state.coolCentroids(i,:)).^2, 2);
      end
      [~, state.coolPattern(ind)] = min (dists, [], 2);
    else
      state.airconPower(ind) = 0;
      state.cooling(ind) = 0;
      state.coolPattern(ind) = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Estimate power used on heating
    if meta.do_heating
      fprintf('heating... ');

      % Power consumption on cold days ("heating") and hot non-vacations ("heating_v")
      [i2, i1] = meshgrid(1:meta.SamPerDay, meta.coldDays);
      s2i = sub2ind ([meta.Days meta.SamPerDay], i1(:), i2(:));
      cold_no_hw_solar = reshape(data(:,s2i) + full(s.cor_hw(ind,s2i)...
						    + cor_solar(:,s2i)),...
				 [size(data,1), length(s2i)/meta.SamPerDay,...
				  meta.SamPerDay]);
      heating    = squeezeNaN(@mean, cold_no_hw_solar, 2);
  %    heating_v  = squeezeNaN(@mean, cold_no_hw_solar ...
  %                 + repmat(full(s.vacation(ind,meta.coldDays)),...
  %                          [1,1,meta.SamPerDay]), 2);
      clear cold_no_hw_solar;

      % Subtract power consumption on mild days ("heating") and mild non-vacations ("heating_v"), to estimate the actual heating power
      %state.cold (ind,:) = reshape(heating, [size(heating,1), meta.SamPerDay]);
      heating   = heating - mild;
      %heating_v = heating_v - mild_v;

      fprintf('Done\n');
      %heatPower_(:,1) = 2*squeezeNaN(@max,heating_v,2);   % '2*' to convert to kW

      %    % Find steps
      %    % use loop to avoid creating huge array
      %   % Ignore 11pm spikes from hot water
      %    %meta.peakhours = (4*meta.SamPerDay/24-1):(23*meta.SamPerDay/24);
      %    ph = length(meta.peakhours);
      %    for i = 1:size(data,1)
      %   userData = 2*(squeeze(data(i,hotDays,meta.peakhours))'+ones(ph,1)*vacation(i,hotDays));
      %   userData = min(userData(1:end-1), userData(2:end));
      %   if all(isnan(userData))
      %       continue
      %   end
      %   userData = [userData userData(end)];    % make up to full size
      %
      %   [dev use alt conf jumps bins] = bigDevices(userData);
      %   heatPower_ (i,3) = alt(1,1);
      %   heatPower_2(i,3) = alt(1,2);
      %   heatPower_ (i,2) = alt(2,1);
      %   heatPower_2(i,2) = alt(2,2);
      %   binsj = bins{1};
      %   binsl = bins{2};
      %
      %   fprintf('Estimated Heating power of %d: %f %f(%f) %f(%f)\n', i, heatPower_(i,1), heatPower_(i,2), heatPower_2(i,2), heatPower_(i,3), heatPower_2(i,3));
      %
      %    end            % for

      %    state.heatPower(ind) = squeezeNaN(@mmedian,(heatPower_+heatPower_2),2);

      %    summary.heatpower (ind)= state.heatPower(ind);
      %summary.heatenergy(ind)= sum(state.heating(ind,:),2);

      %    %%% Partial cluster these users, to 10 clusters
      %    normHourly = bsxfun(@rdivide, heating, squeezeNaN(@mean, heating,2));
      %    n = 10;
      %    [classes, centroids, sumd] = kmeans (normHourly, n, 'EmptyAction','singleton', 'Replicates', 100);
      %    counts = zeros(1,n);
      %    for i = 1:n
      %        counts(i) = sum(classes == i);
      %    end
      %    state.heatCentroids(ind,:) = centroids;
      %    state.heatCounts(ind,:) = counts;
      state.heating(ind) = max(sum(heating,2),0);
      heating = max(0, heating);
      normHourly = bsxfun(@rdivide, heating, squeezeNaN(@max, heating,2)+1e-6);
      n = 6;
      if ~isfield(state, 'heatCentroids') % if first pass, generate clusters
	  [classes, centroids] = kmeans (normHourly, n, 'EmptyAction','singleton', 'Replicates', 100);
	  for i = 1:n
	      counts(i) = sum(classes == i);
	  end
	  [~, sorted] = sort (-counts);
	  state.heatCentroids = centroids(sorted,:);

	  fout = fopen (meta.clusterCSV, 'a');
	  if ~fout
	      error ('Error opening %s for appending\n', meta.clusterCSV);
	  end
	  fprintf(fout, 'Heating\n');
	  for i = 1:n
	      fprintf(fout, '%g,', centroids(i,1:end-1));
	      fprintf(fout, '%g\n', centroids(i,end));
	  end
	  fclose(fout);
      end
      for i = 1:n
	  dists(:,i) = sum (bsxfun(@minus, normHourly, state.heatCentroids(i,:)).^2, 2);
      end
      [~, state.heatPattern(ind)] = min (dists, [], 2);
    else
      state.heating(ind) = 0;
      state.heatPattern(ind) = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % estimate stay-at-home
    days = intersect(meta.weekends,s.mildDays);
    [i2, i1] = meshgrid(1:meta.SamPerDay, days);
    s2i = sub2ind ([meta.Days meta.SamPerDay], i1(:), i2(:));
    WEmildPatterns = reshape(data(:,s2i) + full(s.cor_hw(ind,s2i) + cor_solar(:,s2i)), [size(data,1), length(s2i)/meta.SamPerDay, meta.SamPerDay]);
    WEmildPatterns = squeezeNaN(@mean, WEmildPatterns,2);

    days = intersect(meta.weekdays,s.mildDays);
    [i2, i1] = meshgrid(1:meta.SamPerDay, days);
    s2i = sub2ind ([meta.Days meta.SamPerDay], i1(:), i2(:));
    WDmildPatterns = reshape(data(:,s2i) + full (s.cor_hw(ind,s2i) + cor_solar(:,s2i)), [size(data,1), length(s2i)/meta.SamPerDay, meta.SamPerDay]);
    WDmildPatterns = squeezeNaN(@mean, WDmildPatterns,2);

    %    WEmildPatterns = squeezeNaN(@mean, data(meta.other_users,intersect(meta.weekends,s.mildDays),:), 2);
    %    WDmildPatterns = squeezeNaN(@mean, data(meta.other_users,intersect(meta.weekdays,s.mildDays),:), 2);
    % probably stay-at-home if
    % (a) time of morning peak is similar on weekend and weekdays, and
    % (b) the business-hours use is not much below morning peak, and
    % (c) there is a lunch-time peak.
    % (d) the difference between weekday and weekend profile is "small".

    % (a) Times of peaks
    [~, WEamPeakTime] = max (WEmildPatterns(:,5*meta.SamPerDay/24:11*meta.SamPerDay/24), [], 2);
    [~, WDamPeakTime] = max (WDmildPatterns(:,5*meta.SamPerDay/24:11*meta.SamPerDay/24), [], 2);
    aa = max (WEamPeakTime - WDamPeakTime, 0);

    % (b) business hours vs morning peak
    business = mean (WDmildPatterns(:,9*meta.SamPerDay/24:17*meta.SamPerDay/24), 2);
    %bb = max(WDamPeak ./ business, 1);

    % (c) lunch-time peak
    lunchPeak = max (WDmildPatterns(:,11*meta.SamPerDay/24:14*meta.SamPerDay/24), [], 2);
    cc = max(lunchPeak ./ business, 1);

    % (d) difference between weekend and weekday
    dd = abs(squeezeNaN(@sum, (WEmildPatterns - WDmildPatterns) ./ (WEmildPatterns + WDmildPatterns), 2));
    dd(isnan(dd)) = 0;

    %state.stayAtHomePeakDiff (ind) = aa/2;
    %state.stayAtHomeLunchPeak(ind) = 1 ./ cc;
    %state.stayAtHomeWeekend  (ind) = iif(isnan(dd), 0, dd / 2);
    %state.stayAtHome(ind) = (1 ./ cc + aa/2 + dd / 2) < 3;
    stayAtHome = (1 ./ cc + aa/2 + dd / 2) < 3;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Record average power consumption at each time of day
    ave  = squeezeNaN(@mean, data, 2);
    hot  = squeezeNaN(@mean, data(:,meta.hotDays, :), 2);
    cold = squeezeNaN(@mean, data(:,meta.coldDays,:), 2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Total solar power that was consumed locally instead of exported
    slr = ismember(ind, state.solar_users);
    state.solarUsed(ind(slr)) = sum(cor_solar(slr,:) ...
                                    + double(min(0,data(slr,:))), 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output stats
        % hot water users are listed in order.  Track them to avoid search
        % Add sentinels
peak   = meta.SamPerDay/24*max(data(:,:),[],2);
peak39 = meta.SamPerDay/24*max(max(data(:,:,15*meta.SamPerDay/24:21*meta.SamPerDay/24),[],2), [], 3);
fprintf ('Outputing data for this pass (%d rows)... ', length(ind));
    if ~isempty(meta.hw_users)
        hw_ptr = find(s.hw_users == meta.hw_users(1)+meta.UserOffset);
    else
        %hw_ptr    = 1;
        hw_ptr = find (s.has_hws_est > meta.UserOffset, 1);
        if isempty (hw_ptr)
            hw_ptr = 1;
        end
    end
    if isempty (state.hw_users) || state.hw_users(end) >= 0
        state.hw_users   (length(state.hw_users)+1) = -1;
    end
    if isempty (state.has_hws_est) || state.has_hws_est(end) >= 0
        state.has_hws_est(length(state.has_hws_est)+1) = -1;
    end

    su = meta.solar_users;
    if ~isempty(su)
        if isnan (su)           % if we don't trust ESD, solar = producer
            su = find (min (data(:,:), [], 2) < 0);
        end
        first_solar_data = find(s.solar_users == su(1)+meta.UserOffset);
    else
        first_solar_data = 1;
    end
    solar_ptr = first_solar_data;
    if isempty (state.solar_users) || state.solar_users(end) >= 0
        state.solar_users(length(state.solar_users)+1) = -1;
    end

    for i = 1:length(ind)
        j = i-1 + ind(1);
        if state.solar_users(solar_ptr) == j
            sol = solar_ptr;
            solar_ptr = solar_ptr + 1;
        else
            sol = 0;
        end
        %if state.hw_users(hw_ptr) == j
        if state.has_hws_est(hw_ptr) == j
            hw = hw_ptr;
            hw_ptr = hw_ptr + 1;
        else
            hw = 0;
        end
        fprintf(fout_s, '%s,', state.NMIs{j});
        fprintf(fout_s, '%g,', ...
            (sol ~= 0) + 2*(hw ~= 0), ...
            state.daysValid(j), ...
            peak(i), ...
            peak39(i), ...
            (peak(i) > 4.5) + (peak(i) > 6), ...
            sum(ave(i,:)), ...
            sum(hot(i,:)), ...
            sum(cold(i,:)), ...
            state.airconPower(j), ...
            state.cooling(j), ...
            state.coolPattern(j), ...
            state.heating(j), ...
            state.heatPattern(j), ...
            state.emptyDays(j), ...
            stayAtHome(i) ...
            );
        if state.poolPump(j) > 0
            fprintf(fout_s, '%g,%g,', full(state.poolPump(j)), ...
                -full(sum(cor_pool(i,:)))/sum(~isnan(data(i,:,1))));
            fprintf(state.fout_pool, '%s,', state.NMIs{j});
            fprintf(state.fout_pool, '%g,', reshape(-full(cor_pool(i,:)), [meta.Days,meta.SamPerDay])');
            fprintf(state.fout_pool, '\n');
        else
            fprintf(fout_s, ',,');
        end
        if hw
            fprintf(fout_s, '%g,%g,', -2*full(min(state.cor_hw(j,:))), -full(sum(state.cor_hw(j,:)))/sum(~isnan(data(i,:,1))));
            fprintf(state.fout_hw, '%s,', state.NMIs{j});
            fprintf(state.fout_hw, '%g,', reshape(-full(state.cor_hw(j,:)), [meta.Days,meta.SamPerDay])');
            fprintf(state.fout_hw, '\n');
        else
            fprintf(fout_s, ',,');
        end
        if sol
            fprintf(fout_s, '%g,%g,%g,%g', 2*full(max(cor_solar(i,:))), full(sum(cor_solar(i,:))), full(state.solarUsed(j)), full(state.solar_az(sol)));
            fprintf(state.fout_solar, '%s,', state.NMIs{j});
            fprintf(state.fout_solar, '%g,', reshape(-full(cor_solar(i,:)), [meta.Days,meta.SamPerDay])');
            fprintf(state.fout_solar, '\n');
        else
            fprintf(fout_s, ',,,');
        end
        fprintf(fout_s,',%g', ave(i,:));
        fprintf(fout_s,',%g', hot(i,:));
        fprintf(fout_s,',%g', cold(i,:));
        fprintf(fout_s, '\n');
    end

    % Remove sentries
    if (state.has_hws_est(end) == -1)
      state.has_hws_est(end) = [];
    end
    if (state.solar_users(end) == -1)
      state.solar_users(end) = [];
    end

    fprintf ('Done\n');
end

end

% vi: set expandtab
