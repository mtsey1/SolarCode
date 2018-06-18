function [su, cor_solar] = calc_cor_solar (data, s, meta, ind)
  cor_solar = sparse (length (ind),meta.Days*meta.SamPerDay);
  su = meta.solar_users;
  first_solar_data = 1;
  
  % to investigate:
  dubious = [111, 557, 559,  709, 716, 794, 1027, 1058];
  % (Possible features: max-min range of daily min / jump of rolling_min
  %                     flatness of demand during sunlight hours)
  % 74: second "discon" is a pool.
  % 150: solar rarely works.  "surround" contains non-solar days, due to "gap" constraint.
  % 208: Solar efficiency drops massively in winter.
  % 211: Too westerly.  Possibly too high a capacity.
  % 244: Too easterly?  Solar for very short duration.  Shade?  Second gap wrong
  % 269: First FP. Fifth FN (given factor of 10).  Should merge over NaN
  % 315: Too easterly.  Too high capacity?
  % 380: Too westerly.  Too high capacity, or capFactor1 too high.
  % 388: First, 2nd FP -- adjacent to vacation.  Too easterly -- fixed by removing pool?
  % 392: FP.  (Const reading all day: 0.225)
  % 414, 415: FP
  % 527: Is second unmetered?
  % 549: capacity way too low?
  % 604: FN (due to end of winter?)
  % 612: FP, even with factor of 10.  Curved minimum
  % 644: FP. Curved minimum
  % 667: Third is FN
  % 669: both FP
  % 683: FN
  % 735: FP due to pool pump
  % 737: FN
  % 758: FN
  % 817: Generation drops at day 50
  % 835: FP
  % 837: Second FP (airconditioning day), 3rd FP
  % 872: fourth FP (adjacent to vacation)
  % 872: FP (day with background load, and swimming pool during sunlight)
  % 882: FN maybe caused by too easterly
  % 889: FP
  % 917: FN -- almost whole time.
  % 922: Second FN
  % 928: FN -- almost whole time.  Includes large unmetered time
  % 965: FN at end
  % 995: FN at end
  % 1015: FN at end
  % 1018: "unmetered" for too long
  % 1025: FN
  % 1041: FN at end
  % 1050: Data negated up to day 127, right from 128
  % 1055: FP due to pool
  % 1073: FN at end
  % 1098: FN at end, 1e-112
  % 1112: FN.  Maybe unmetered rather than disconnected.
  % 1177: Capacity seems to increase at about day 75
  % 1183: FN
  % 1226: FP  All-day devices
  % 1267: FN at end.  Intermittent fault
  % 1336: FP. All-day devices
  % 1359: FN
  % 1360: FP?
  % 1407: FP?
  % 1429: FP?
  % 1457: FN (reading constant)
  % 1468: FN. Intermittent fault
  % 1472: unmetered starts too early
  % 1474: FN.  All-day device on adjacent days
  % 1475: FN (reading constant)
  % 1525: FP.  All-day device
  % 1538: FN
  % 1541: FN
  % 1548: Intermittent fault.
  % 1564: FN at end.  Start matches big winter load
  % 1618: FP
  % 1745: First FN?
  % 1762: FP
  % 1804: Data negated up to day 194, right from 195
  % 1821: FN (reading constant)
  % 1848: Frequent intermittent fault?
  % 1880: FN (reading constant).  Massive over-estimation of capacity
  % 1886: FN (reading constant)
  % 1931: FN???
  % 1967: FN (reading constant)
  % 2016: FN at end.
  % 2023: FN (reading constant)
  % 2033: FN -- long run
  dubious_2 = [62, 74, 104, 150, 208, 211, 244, 269, 315, 380, 388, 392, ...
               414, 415, 527, 549, 604, 612, 644, 667, 669, 683, ...
               735, 737, 758, 817, 835, 837, 872, 882, 889, ...
               917, 922, 928, 965, 995, ...
               1015, 1018, 1025, 1041, 1050, 1055, 1073, 1098, ...
               1112, 1177, 1183, 1226, 1267, 1270, 1336, 1359, 1360, ...
               1406, 1407, 1429, 1457, 1468, 1472, 1474, 1475, ...
               1525, 1538, 1541, 1548, 1564, 1618, ...
               1745, 1762, 1804, 1821, 1848, 1880, 1886, 1931, 1967, ...
               2016];
  % 234: Missed HWS, since its turn-on time changes
  dubious_3 = [234];
  
  if ~isempty (su)
    if all (isnan (su))           % if we don't trust ESD, solar = producer
      su = find (min (data(:,:), [], 2) < 0);
    end
    s_su = su + meta.UserOffset;
    first_solar_data = find (s.solar_users == su(1)+meta.UserOffset);
    solar_user_range = first_solar_data:first_solar_data + length (su) - 1;
    solar_correction= zeros (length (su), meta.Days, s.dark_start-s.dark_end+1);
    
    capFactor1(length (setdiff (s.postcode(s_su), [])), ...
               meta.Days,s.dark_start-s.dark_end+1) = 0;
    for k = 1:length (meta.pclist)
      u = find (s.postcode(s_su) == meta.pclist(k));
      %if k==69 keyboard end
      if ~isempty (u)
        capFactor1(u,:,:) = s.capFactor(first_solar_data-1+u,:,:) ...
                           .* repmat (s.solar_by_pc(k,:,:), [length(u), 1, 1]);
        solar_correction(u,:,:) = bsxfun (@times,capFactor1(u,:,:), ...
                                          s.solar_cap(solar_user_range(u))');
      end
      %fprintf ('%f\t', solar_correction(3,22,12));
    end
    
    % second return val of "max" is the *first* occurrence of "true"
    [~,first_solar_day] = max (s.daily_min (solar_user_range,:) < 0,[],2);
    [~,last_solar_day]  = max (s.daily_min (solar_user_range,end:-1:1) < 0,[],2);

    % (~< instead of >= to catch NaN)
    disconnected = ~(s.daily_min (solar_user_range,:) < 0)';
    % users for which not all days generate.
    % Precalculate to avoid calling rolling_min on fully-connected users
    not_all = (any (disconnected));
    gap = 11;         % If no net generation for >gap days, assume no solar
            % gap must be odd.
    disconnected(:,not_all) = -rolling_min (-rolling_min (disconnected(:,not_all), gap), gap);
    not_all = (any (disconnected));     % recalc after "noise" removed

    % More restrictive version: only disconnected if load above
    % vampires, not above zero.  Not yet used.
    daily_min2 = squeeze (min (bsxfun (@minus, ...
                                       data(su,:,:), ...
                                       cumsum (single (full (s.vampires(s_su,:))), 2)), ...
                             [], 3));

    % (~< instead of >= to catch NaN)
    disconnected2 = ~(daily_min2 < 0 & s.daily_min (solar_user_range,:) ~= 0)';
    not_all2 = (any (disconnected2));
    gap = 11;         % If no net generation for >gap days, assume no solar
            % gap must be odd.
    disconnected2(:,not_all2) = -rolling_min (-rolling_min (disconnected2(:,not_all2), gap), gap);
    not_all2 = (any (disconnected2));     % recalc after "noise" removed

    % TODO: Check that days seeming to lack solar really do.
    % TODO: What do we do with days with solar but no "generation" meter?
    % - histogram very peaked?
    % - histogram very different from neighbours?
    % - var (daily_min) much smaller than neighbours?
    % - step change in daily_min?   (could be step change in demand)
    % - period at start or at end of measurements?
    % - gap very long

    discon = false (size(data,1), meta.Days);
    light_hrs = 15:36;
    light_h_cor = light_hrs - s.dark_end + 1;
solar_fails = fopen ('solar_fails.txt', 'w+');
    for i = find (not_all)
if i < 1746, continue; end
      % Find days on which it seems the meter reads 0 for nett generation.
      day_min = min (squeeze (data(su(i), :, light_hrs)), [], 2)';
      hits_zero = find (day_min == 0, 1);
      if ~isempty (hits_zero) && hits_zero < first_solar_day(i)
        unmetered = hits_zero:first_solar_day(i) - 1;
        %first_solar_day(i) = unmetered(1);
        
%{
        small = (day_min >= 0 & day_min < 0.1);
        unmetered = (hits_zero | (small & (hits_zero([1, 1:end-1]) ...
                     | hits_zero([2:end, end]))));
        unmetered = find (unmetered);
        if unmetered(1) < first_solar_day
          first_solar_day = unmetered(1);
        end
        if unmetered(end) > last_solar_day
          last_solar_day = unmetered(end);
        end
%}
        %unmetered = rolling_min (unmetered, 3);
        % What do we do with unmetered?
        % pos(i) = (mean (-i:0) * (len-i+1) + mean (len+1:len+1+i) * (i-1)) / len
        % if (i,j) > 0.1, add solar.  If < 0.1, interpolate up and down
        % Update first_solar_day and/or last_solar_day.
        % Merge unmetered runs if gaps small?
      end
      
      % Find days on which we are disconnected.
      disconnected(any  (isnan (data(su(i),:,:)), 3), i) = true;
      disconnected2(any (isnan (data(su(i),:,:)), 3), i) = true;
      d = disconnected(:,i);
      fd = find (d);
          % If "disconnected" all at start, probably new user
          % Else, check if true disconnect, or all generation used
      if max (fd) ~= length (fd)
        breaks = find (diff([-1;fd;-1],2));
        rr = ranges (fd);
        for r = 1:size(rr, 2)
          if all (all (isnan (data(su(i), rr(1,r):rr(2,r), light_hrs)))) ...
              || all (all (isnan (solar_correction(i, rr(1,r):rr(2,r), ...
                                                   light_h_cor))))
            continue;
          end
          surround = false(meta.Days, 1);
          surround(max (1, rr(1,r) - gap):rr(1,r)-1) = true;
          surround(rr(2,r)+1:min (meta.Days, rr(2,r) + gap)) = true;
	        surround = surround & ~disconnected(:,i);
          surround_dat = data(su(i), surround, light_hrs) ...
                         + solar_correction (i, surround, light_h_cor);
          surround_dat = surround_dat(:);
          
          inside = data(su(i), rr(1,r):rr(2,r), light_hrs);
          [~, ks_without] = kstest2 (surround_dat, inside(:));
          
          inside = inside + solar_correction(i, rr(1,r):rr(2,r), light_h_cor);
          [~, ks_with] = kstest2 (surround_dat, inside(:));

miss_range = [min(s.daily_min(i,rr(1,r):rr(2,r))), ...
	      max(s.daily_min(i,rr(1,r):rr(2,r)))];
miss_jump = (miss_range(1) + miss_range(2)) / 2 ...
	    - mean (s.daily_min(i,surround));
aaa(s.dark_end:s.dark_start, :) = squeeze (solar_correction(i, :, :))';
bbb = aaa;
bbb(:, rr(1,r):rr(2,r)) = 0;
cv = (squeeze (data(su(i), :, :))');
figure(100); imagesc (cv);
figure(102); imagesc (max (0, cv + bbb));
figure(103); imagesc (max (0, cv + aaa));
figure(101); plot ([s.daily_min(i,:)']);
figure(104); plot (cv(:, rr(1,r):rr(2,r)));
fprintf ('with %g  without %g  ratio %g  r %d\n', ...
	 ks_with, ks_without, ks_with / ks_without, r);
disp (rr);
ny = 'NY';
my_guess = (miss_range / miss_jump < 0.2);
p = 'p;
while p == 'p' || p == 'P'
  p = input ('Solar missing? (Y/N/U/P) ', 's');
  if isempty(p), p = ny(my_guess+1); end
  if length(p) > 1,  p=p(1); end
  if p == 'p'
    keyboard;
  end
end
switch upper(p)
  case {'Y','YES'}, is_miss = 1;
  case {'N','NO'},  is_miss = 0;
  case {'U'},       is_miss = -1;
  otherwise
      fprintf('Unknown option.  Assuming  Unknown\n');
      is_miss = -1;
end
fprintf (solar_fails, '%d %d %d %d %g %g %g %g %d\n', ...
	 i, rr(1,r), rr(2,r), is_miss, ...
	 miss_range, miss_jump, ks_with, ks_without, is_miss ~= my_guess);
          
          if (ks_with < ks_without )
            % Probably disconnected
            confirmed = (ks_with * 10 < ks_without);
            if ~confirmed
              % TODO: repeat kstest2 on disconnected2?
            end
            % TODO: Check for no overlap with unmetered.
            if confirmed
    	        discon(i,rr(1,r):rr(2,r)) = true;
            end
          end
        end
%{        
        % for each range
        %   Kolmogorov-Smirnoff test whether match is better with or without
        %        solar correction.
        %   If match to surround is better with correction, ignore this range
        %   If match is better without correction, set discon to true.
        surround = logical(-rolling_min (-disconnected(:,i), 2*gap)) ...
	           & ~disconnected(:,i);
        min_dis = squeeze (min (data(su(i),d,       light_hrs),[],3));
        min_sur = squeeze (min (data(su(i),surround,light_hrs),[],3));
        v_d = var (min_dis);
        v_c = var (min_sur);
        if v_c < 20 * v_d
          a = (squeeze (data(su(i), disconnected(:,i),light_hrs))');
          a = bsxfun (@minus, a, ...
                 cumsum (single (full (s.vampires...
                          (i, disconnected(:,i)))), 2));
          h = hist(a(:),20);

        end
	      discon(i,fd) = true;
%}        
      end
    end
fclose (solar_fails);

    last_solar_day = meta.Days - last_solar_day + 1;
    for i = find (first_solar_day ~= 1)
      solar_correction(i,1:first_solar_day(i)-1,:) = 0;
    end
    for i = find (last_solar_day ~= meta.Days)
      solar_correction(i,last_solar_day(i)+1:end,:) = 0;
    end
    for i = find (any (discon, 2))'
      solar_correction(i, discon(i,:), :) = 0;
    end
    
    tmp = zeros (length (su), meta.Days, meta.SamPerDay);
    tmp(:, :, s.dark_end:s.dark_start) = max (solar_correction,0);
    % any demand less than vampires is due to solar
    % (state.vampires is the diff, to make sparse matrix smaller)
    tmp = max (tmp, -(bsxfun (@minus, data(su,:,:), ...
                cumsum (single (full (s.vampires(s_su,:))),2))));
    cor_solar(su,:) = sparse (double (tmp(:,:)));
  end
end
