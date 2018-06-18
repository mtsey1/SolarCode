function [su, cor_solar] = calc_cor_solar (data, s, meta, ind)
  % 0: use pre-trained SVM.
  % 1: confirm guesses of pre-trained SVM
  % 2: Learn new SVM from confirmed guesses
  learn_fail_mode = 0;
  
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

    disconnected = (s.daily_min (solar_user_range,:) >= 0)';
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
    disconnected2 = (daily_min2 >= 0 | s.daily_min (solar_user_range,:) == 0)';
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
    if learn_fail_mode == 1
      solar_fails = fopen ('solar_fails.txt', 'w+');
    elseif learn_fail_mode == 2
      load solar_fails.txt
      count = 0;
    else
      beta = [0.0361   -1.4144   -0.5224   -0.2956   -0.5230   -0.0302   -1.3614    0.6449    1.0665 0.7955   -0.6134    0.1866]';
      thresh = 4.25;
    end
    vampires = cumsum (single (full (s.vampires(s_su,:))),2);
    count = 0;
    for i = find (not_all)
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
          
          inside_plus = inside + solar_correction(i, rr(1,r):rr(2,r), light_h_cor);
          [~, ks_with] = kstest2 (surround_dat, inside_plus(:));
          if learn_fail_mode < 2
            count = count + 1;
          else
            count = find (solar_fails(:,2) == i & solar_fails(:,4) == rr(2,r));
            if ~any (count)
              continue
            end
          end

          % trim end bits that shouldn't have been included
          orig_len = rr(2,r) - rr(1,r);
          med = median (s.daily_min (solar_user_range(i), rr(1,r):rr(2,r)));
          while (rr(2,r) - rr(1,r) > orig_len / 2)
            m = min (s.daily_min (solar_user_range(i), rr(1,r):rr(2,r)));
            if (m == s.daily_min (solar_user_range(i), rr(1,r)) && m < 0.5 * med) || ~isfinite (s.daily_min (solar_user_range(i), rr(1,r)))
              rr(1,r) = rr(1,r) + 1;
            elseif (m == s.daily_min (solar_user_range(i), rr(2,r)) && m < 0.5 * med) || ~isfinite (s.daily_min (solar_user_range(i), rr(1,r)))
              rr(2,r) = rr(2,r) - 1;
            else
              break
            end
          end
          solar_fails(count, [3,4]) = rr(:,r);
          starts(count,1) = rr(1,r);
          ends(count,1) = rr(2,r);
          customers(count,1) = i;

          if sum (surround) < gap
            surround = false(meta.Days, 1);
            surround(max (1, rr(1,r) - 2*gap):rr(1,r)-1) = true;
            surround(rr(2,r)+1:min (meta.Days, rr(2,r) + 2*gap)) = true;
            surround = surround & ~disconnected(:,i);
          end
          excursion(count,1) = atan (mean (s.daily_min(solar_user_range(i), surround)) / mean (s.daily_min(solar_user_range(i), rr(1,r):rr(2,r))));

          minus_vamp = s.daily_min(solar_user_range(i), rr(1,r):rr(2,r)) - min(data(su(i),rr(1,r):rr(2,r),[1:8,43:48]), [], 3);
          if rr(2,r) - rr(1,r) > 8
            miss_range(count,1:2) = quantile (s.daily_min(solar_user_range(i), rr(1,r):rr(2,r)), [0.25, 0.75]);
            miss_range2(count,1:2) = quantile (minus_vamp, [0.25, 0.75]);
          else
            miss_range(count,1:2) = [min(s.daily_min(solar_user_range(i),rr(1,r):rr(2,r))), ...
                                     max(s.daily_min(solar_user_range(i),rr(1,r):rr(2,r)))];
            miss_range2(count,1:2) = [min(minus_vamp), max(minus_vamp)];
          end
          jump_base = mean (s.daily_min(solar_user_range(i),surround));
          if rr(1,r) > 1
            jump_base = [jump_base, s.daily_min(solar_user_range(i), rr(1,r) - 1)];
          end
          if rr(2,r) < size (data, 2)
            jump_base = [jump_base, s.daily_min(solar_user_range(i), rr(2,r) + 1)];
          end
          miss_jump(count) = (miss_range(count,1) + miss_range(count,2)) / 2 - mean (jump_base, 'omitnan');
          ks_ww(count,1:2) = [ks_with, ks_without];
          len(count, 1) = rr(2,r) - rr(1,r);
          % Look for big jumps indicating devices on timers (e.g., pool pumps)
          pump = 1;
          b      = quantile (inside,      [0.1, 0.5], 2);
          b_plus = quantile (inside_plus, [0.1, 0.5], 2);
          base      = b(:,1,:);
          base_plus = b_plus(:,1,:);
          st = [];
          st(2,:) = base(3:end) - base(1:end-2);
          st(1,1:end+1) = diff (base);
          st = st(:);
          [down, down_pos] = min (st(end/2+1:end));
          [up, up_pos] = max (st(1:end/2));
          % skip wake-up peak
          if up > 0.3 && up_pos < 5 ...
              && ~all (base_plus (up_pos:length(base)/2) > 0.9 * up)
            [up1, up_pos1] = max (st(5:end/2));
            if up1 > 0.3
              up = up1;
              up_pos = 4 + up_pos1;
            end
          end

          up_pos = ceil ((up_pos + 3) / 2);
          down_pos = floor ((down_pos + length (base)) / 2);
          step_mx = max (up, -down);
          step_mn = min (up, -down);
          if step_mn > 0.3 && down_pos - up_pos > 2 ...
             && all (base_plus (up_pos:down_pos) > 0.9 * step_mx) ...
             && max (base_plus (up_pos:down_pos)) - min (base_plus (up_pos:down_pos)) < 2 * step_mx
            cor = zeros(size(base));
            cor(up_pos:down_pos) = (step_mn + step_mx) / -2;
            b = b + cor;
            b_plus = b_plus + cor;
            inside      = bsxfun (@plus, inside,      cor);
            inside_plus = bsxfun (@plus, inside_plus, cor);
          end
          TVratio(count,1:2) = log (sum (abs (diff (b, 1, 3)), 3) ...
                                 ./ sum (abs (diff (b_plus, 1, 3)), 3));
          TVratio_inr(count,1) = log (sum (abs (diff (b(:,1,6:end-6), 1, 3)), 3) ...
                                   ./ sum (abs (diff (b_plus(:,1,6:end-6), 1, 3)), 3));
          TVratio_all(count,1) = log (sum (sum (abs (diff (inside, 1, 3)))) ...
                                   ./ sum (sum (abs (diff (inside_plus, 1, 3)))));
          vamp_ratio(count,1) = atan (min (vampires (i,rr(1,r):rr(2,r))) / min (b(:,1,:)));
          % If round-off, don't take log of negative.
          flatness(count,1) = log ((mean (inside_plus(:), 'omitnan') + 1e-3 - min (inside_plus(:), [], 'omitnan')) ...
                                ./ (mean (inside(:),      'omitnan') + 1e-3 - min (inside(:),      [], 'omitnan')));
          negatives(count,1) = sum (any (inside < 0, 3)) / size (inside, 2);

          if rr(1,r) > size (data, 2)/2 || rr(2,r) < size (data, 2)
            summerness(count,1) = abs (size (data, 2) - (rr(1,r) + rr(2,r))) / 2;
          else
            summerness(count,1) = mean (abs (size (data, 2) - rr(:,r)));
          end

          % TODO: Look for jumps (timed devices)
          % TODO: better "miss jump"
          % TODO: consider vampires.
          % TODO: Consider time of year: false positives more likely in winter or hot days
          if learn_fail_mode == 1
            miss_range = [min(s.daily_min(solar_user_range(i),rr(1,r):rr(2,r))), ...
                          max(s.daily_min(solar_user_range(i),rr(1,r):rr(2,r)))];
            miss_jump = (miss_range(1) + miss_range(2)) / 2 ...
                        - mean (s.daily_min(solar_user_range(i),surround));
            aaa(s.dark_end:s.dark_start, :) = squeeze (solar_correction(i, :, :))';
            aaa(48,1) = 0;
            bbb = aaa;
            bbb(:, rr(1,r):rr(2,r)) = 0;
            cv = (squeeze (data(su(i), :, :))');
            figure(100); imagesc (cv);
            figure(102); imagesc (max (0, cv + bbb));
            figure(103); imagesc (max (0, cv + aaa));
            figure(101); plot ([s.daily_min(solar_user_range(i),:)']);
            figure(104); plot (cv(:, rr(1,r):rr(2,r)));
            fprintf ('with %g  without %g  ratio %g  r %d\n', ...
                     ks_with, ks_without, ks_with / ks_without, r);
            disp (rr);
            ny = 'NY';
            my_guess = ((miss_range(2)-miss_range(1)) / miss_jump < 0.2);
            p = 'p';
            while p == 'p' || p == 'P'
              p = input (['Solar missing? (Y/N/U/P) [' ny(my_guess+1) '] '], 's');
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
          end
%{
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
%}
        end
      end
    end
    
    mr = min (diff (miss_range, 1, 2), diff (miss_range2, 1, 2));
    factor1 = log (ks_ww(:,1)) ./ (log (ks_ww(:,2)) - 3);
    factor2 = log (mr ./ max (0.01, miss_jump'));
    XX = [factor1, factor2, TVratio, TVratio_inr, TVratio_all, excursion, flatness, atan((summerness - 60) / 30), vamp_ratio, negatives, starts == 1]; 
    XX(isnan(XX)) = 0;
    XX = min (1e3, max (-1e3, XX));

    if learn_fail_mode == 1
      fclose (solar_fails);
    elseif learn_fail_mode == 2
      fails  = (solar_fails(:,5) == 1);
      nfails = (solar_fails(:,5) == 0);
      YY = solar_fails(:, 5);
      a = fails | nfails;   % those not "unknown"
      X = XX(a, :);
      Y = YY(a);
      SVMmodel = fitcsvm (X, Y);
      beta = SVMmodel.Beta;
      ff  = find ( fails);
      fnf = find (nfails);
      [~,  idx] = sort ( XX(ff, :) * beta);
      [~, nidx] = sort (-XX(fnf,:) * beta);
      j = 1;
      %ii = fnf(nidx(j));
      ii = ff(idx(j));
      i = su(solar_fails(ii, 2));
      cv = squeeze (data (i, :, :))';
      ccv = cv;
      ccv(s.dark_end:s.dark_start, :) = ccv(s.dark_end:s.dark_start,:) + squeeze (solar_correction(i, :, :))';
      figure(100); imagesc (squeeze (data(i,:,:))')
      figure(101); plot (1:365, s.daily_min (solar_user_range(i),:))
      figure(104); plot (cv(:, solar_fails(ii,3):solar_fails(ii,4)))
      figure(105); plot (ccv(:,solar_fails(ii,3):solar_fails(ii,4)))
      XX(ii,:) .* beta'
      a = [solar_fails(ii,3), solar_fails(ii,4)]

      figure(3);
      % convert fails and nfails to indices into X
      foo = zeros (size (fails));
      foo( fails) = 1;
      foo(nfails) = 2;
      foo = foo (foo > 0);
      ffails = find (foo == 1);
      fnfails = find (foo == 2);
      % Trim Infs.
      ss = X * beta;
      ss(ss > 500) = max (ss(ss <= 500));
      % plot them
      hold off; plot (excursion( fails), ss( ffails, :), 'x');
      hold all; plot (excursion(nfails), ss(fnfails, :), 'o');
      keyboard
    else
      valid = find (XX * beta > thresh);
      for i = valid(:)'
        discon(customers(i), starts(i):ends(i)) = true;
      end
    end

    last_solar_day = meta.Days - last_solar_day + 1;
    for i = find (first_solar_day ~= 1)'
      solar_correction(i,1:first_solar_day(i)-1,:) = 0;
    end
    for i = find (last_solar_day ~= meta.Days)'
      solar_correction(i,last_solar_day(i)+1:end,:) = 0;
    end
    for i = find (any (discon, 2))'
      cvt_cor = data(i, :, :);
      cvt_cor(:, :, s.dark_end:s.dark_start) = cvt_cor(:, :, s.dark_end:s.dark_start) + solar_correction(i,:,:);
      figure(100); imagesc (squeeze (cvt_cor)');
      figure(104); plot (squeeze (cvt_cor (:, discon(i,:), :))');
      % Should we save these corrections in case of mis-classification?
      solar_correction(i, discon(i,:), :) = 0;
      cvt_cor = data(i, :, :);
      cvt_cor(:, :, s.dark_end:s.dark_start) = cvt_cor(:, :, s.dark_end:s.dark_start) + solar_correction(i,:,:);
      figure(101); imagesc (squeeze (cvt_cor)');
      figure(102); imagesc (squeeze (solar_correction (i,:,:))');
      figure(103); plot (squeeze (data(i, discon(i,:), :))');
      fprintf ('i = %d\n', i);
    end
    
    % Check for over-correction
    df = bsxfun (@minus, vampires/2, data(su, :, s.dark_end:s.dark_start));
    ratio = df ./ solar_correction;
    ratio(:,:) = medfilt1 (ratio(:,:), 3, [], 2);
    ratio (ratio > 1) = 1;
    factor = max (ratio(:,:), [], 2);
    factor(~isfinite (factor)) = 1;
    new_cor = bsxfun (@times, solar_correction, factor);
    solar_correction = new_cor;
%{
    peak = 8:22;
    tmp = data(su, :, s.dark_end:s.dark_start) + solar_correction;
    twi1 = tmp(:, :, [1:peak(1), peak(end):end]);
    twi2 = solar_correction(:, :, [1:peak(1), peak(end):end]);
    twi1 = twi1(:, :)';
    twi2 = twi2(:, :)';
    twi1 = bsxfun (@minus, twi1, mean (twi1, 1));
    twi2 = bsxfun (@minus, twi2, mean (twi2, 1));
    twi1 = bsxfun (@times, twi1, 1 ./ sqrt (sum (twi1.^2, 1, 'omitnan')));
    twi2 = bsxfun (@times, twi2, 1 ./ sqrt (sum (twi2.^2, 1, 'omitnan')));
    c_twi = sum (twi1 .* twi2, 1, 'omitnan');
    
    pk1 = tmp(:, :, peak);
    pk2 = solar_correction(:, :, peak);
    pk1 = pk1(:,:)';
    pk2 = pk2(:,:)';
    pk1 = bsxfun (@minus, pk1, mean (pk1, 1));
    pk2 = bsxfun (@minus, pk2, mean (pk2, 1));
    pk1 = bsxfun (@times, pk1, 1 ./ sqrt (sum (pk1.^2, 1, 'omitnan')));
    pk2 = bsxfun (@times, pk2, 1 ./ sqrt (sum (pk2.^2, 1, 'omitnan')));
    c_pk = sum (pk1 .* pk2, 1, 'omitnan');
    
    over_est = max (0, min (c_twi,c_pk));
    new_cor = bsxfun (@times, 1 - over_est(:), solar_correction);
    too_little = sum (data(su, :, s.dark:s.dark_start) - new_cor < 0, 2);
%}
    
    tmp = zeros (length (su), meta.Days, meta.SamPerDay);
    tmp(:, :, s.dark_end:s.dark_start) = max (solar_correction,0);
    % any demand less than vampires is due to solar
    % (state.vampires is the diff, to make sparse matrix smaller)
    tmp = max (tmp, -(bsxfun (@minus, data(su,:,:), ...
                cumsum (single (full (s.vampires(s_su,:))),2))));
    cor_solar(su,:) = sparse (double (tmp(:,:)));
  end
end
