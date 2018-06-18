function [pump, cor_pump, power, profile, pattern] = has_untimed_pump (userData, vampires, meta, cust, issolar, state)
    record = [];
    pump = 0;
    cor_pump = sparse(size(userData,1), size(userData,2));
    power = 0;
    profile = 0;
    pattern = 0;

    c = (squeeze(userData)');


%    t = find_idle (c);
%    tt = reshape (t, size (c));
%%%    ttt = bwmorph (bwmorph (tt, 'close'), 'open');
%%    figure(1); imagesc (c);
%%    figure(2); imagesc (tt);
%%    figure(3); hist(abs(diff(c(t))), 40);
%%%    figure(3); imagesc (ttt);
%%    figure(4); plot (c(t));
%%%    figure(5); plot (c(ttt(:)));
%%    figure(6); plot (sum(tt'));
%%
%%%    [sum(abs(diff(c(ttt&isfinite(c))))), sum(abs(diff(c(t))))]
%%    spikiness = mean(abs(diff(c(t))));
%%    fprintf ('%g\n', spikiness);
%    
%    %spikiness = mean(abs(diff(c(1:9,:)(:))));
%    spikiness = sum (abs (diff (c(t))) > 0.4) / length (c(t));
%    if spikiness > 0.1
%      figure(1); imagesc (c);
%      figure(2); imagesc (tt);
%      figure(3); hist(abs(diff(c(t))), 40);
%      figure(4); plot (c(t));
%      figure(5); hist(abs(diff(c(1:9,:)(:))), 40);
%      figure(6); plot(c(1:9,:)(:));
%
%      fprintf('\nspikiness %g\n', spikiness);
%      keyboard
%    end
%    
%    power = spikiness;
%    if isempty (pump)
%      pump = 0;
%    end
%    return



    validDays = find(~isnan(c(1,:)));
    if 1
        if length(validDays) > 40
            [pump, power, cor_pump] = find_flat (userData, meta, state, validDays, cust);
        end
    else
        c = stripSpikes (squeeze(userData)');
        %c(:,isnan(vacation)) = NaN;
    
        validDays = find(~isnan(c(1,:)));
    %if 1
        if length(validDays) > 40
            [pump, power, cor_pump] = find_flat (userData, meta, state, ...
                                           validDays, cust);
        end
    %else
    %    c = stripSpikes (squeeze(userData)');
    %    %c(:,isnan(vacation)) = NaN;
    %
    %    validDays = find(~isnan(c(1,:)));
    %    cc = c (:,validDays);
    %
    %    %vampires = find_vampires(cc');
    %    vampires = vampires(validDays);
    %    cv = cc - repmat(vampires,  [size(cc,1), 1]);
    end

    
end

% Find a few features:
%      - jump size
%      - correlation with hot days
%      - correlation with cold days
%      - similarity of turn-on times
%      - similarity of durations
%      - variance during "on" period
% Classify
% Filter bursts not fitting profile
% Recalculate jump size


function [base, spikes] = strip_spikes(c, peakRatio)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wdw = 3;   % 1.5 hours
    base = -rolling_min(-rolling_min(mce, wdw),   wdw);
    if nargout > 1
        spikes = c -base;
    end
end

function [base, spikes] = stripSpikes_single(c, peakRatio)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 2
        peakRatio = 5;
    end
    idx = find(c(2:end-1) - c(1:end-2) > peakRatio * abs(c(3:end) - c(1:end-2)));
    base = c;
    base(idx+1) = (c(idx) + c(idx+2))/2;
    if nargout > 1
        spikes = c -base;
    end
end

function idle_times = find_idle (c)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % find times such that
  % a) one slot before us and one slot after are both low, or
  % b) we are low and one slot before and one slot after are both in class a.
  len = length(c(:));
  after  = [2:len, len];
  before = [1, 1:len-1];

  cp = c(c(:) > 0);
  if length (cp) > 2000
    h = cp(1:13:end);
  elseif length (cp) > 500
    h = cp(1:5:end);
  else
    h = cp;
  end
  if ~isempty (h)
    thresh = quantile (h, 0.2);
  else
    thresh = -10;
  end
  %thresh = 0.9 * min (c(:)) + 0.1 * max (c(:));

  class_a = (c(before) < thresh) & (c(after) < thresh);
  class_b = (c(:)' < thresh) & (class_a(before)  | class_a(after));

  idle_times = class_a | class_b;

  s = [sum(c(:)<thresh)/sum(c(:)< 1000), sum(idle_times)];
end


function [pump, power, cor] = find_flat (data, meta, state, validDays, cust)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %function [base, spikes] = find_flat (data, peakRatio)
    power = 0;
    jump = 0;
    cor = zeros (1, numel(data));
    ccc = data';
    ccc(abs(diff(ccc(:))) > 0.1) = NaN;
    [a, b] = hist(ccc(:),40);
    [p, c, s, r] = secondPeak (a, 3);
%figure(1); log_hist(ccc(:), 40)

%figure(2); imagesc(data');
%figure(10); plot(data(330:end,:)');
%%figure(4); plot(data(150:200,:)');
%%figure(5); plot(ccc(any(isfinite(ccc')),:));
%%b(p)

    % Find "big" peak that isn't from small switches etc
    if abs(p(2) - p(1)) < max(b(p(isfinite(p))))/3
        big = 3;
    else
        big = 2;
    end

    aircon = 1;
    heater = 2;
    pool_pump = 3;
    other = 4;
    type_names = { 'aircon', 'heater', 'pool_pump', 'other' };
    type = other;

    % In decreasing order, peel off "separable" thresholds
    [~, idx] = sort(p);
    dd = data';
    for n = length(idx):-1:2
        m = idx(n);
%fprintf ('Separation %g\n', s(m));
        if s(m) > 0 * 0.5
            if s(m) < 1
              % Should do some checking that the peaks are sufficiently
              % well separated
            end

            thresh = b(r(m,2));       % end of "rectangle" between peaks

%input (sprintf ('%d of %d  thresh %g  ', length(idx)-n+1, length(idx)-1, thresh), 's');
            % skip cases where the "normal case" peak wasn't the first peak
            if (   thresh < 0.5*skipNaN (@mean,data(:)) ...
                && thresh < 0.5*skipNaN (@mmedian, data(:)))
              continue;
            end

            dd(~isfinite(dd)) = 0;
            ddd = dd;
            ddd(:) = -rolling_min(-rolling_min(ddd(:), 3), 3);

            jump = NaN;
            [fj, jt, sj] = find_jumps(ddd);
            if ~any(fj > thresh)
                ddd = dd;
                [fj, jt, sj] = find_jumps(ddd);
            end
            % quick-and-dirty estimate for jump size
            if any (fj > thresh)
                jump = median (fj(fj > thresh));    % Assumes jump >> typical power
            else
                jump = thresh;
            end

            d = dd;
            d(:) = (ddd > thresh);

            %% Better(?) estimate for jump size:
            rr = ranges (find(d(:)));
            % find sj near quick-and-dirty estimate, near start/end of range
            starts = zeros(1,size(rr,2));
            ends   = zeros(1,size(rr,2));
            for i = size(rr,2):-1:1
%if (floor (rr(1,i) / 48) == 69);       %% misses period between two gaps in A/C
%keyboard;
%end;
                st = find (abs (jt - rr(1,i)) <= 3);
                if length (st) > 1
                    [mismatch, m] = min (abs(sj(st) - jump));
                    if mismatch < 0.5 * jump
                        st = st(m);
                    else
                        st = [];
                    end
                end
                if isempty(st)
                    starts(i) = [];
                else
                    starts(i) = st;
                end
                en = find (abs (jt - rr(2,i)) <= 3);
                if length (en) > 1
                    [mismatch, m] = min (abs (sj(en) + jump));
                    if mismatch < 0.5 * jump
                        en = en(m);
                    else
                        en = [];
                    end
                end
                if isempty(en)
                    ends(i) = [];
                else
                    ends(i) = en;
                end
            end
            if isempty (starts)
              continue;
            end

            jump = median (fj([starts, ends]));
            power (n) = jump;
%fprintf ('cust %d jump %g\n', cust, jump);
            %  (broadcast comparison is faster.  Perhaps trim out-of-range
            %  jumps first to conserve memory)
            % jump = mean (selected jumps);
            %
            %% Only match an interval if one end has jump near enough to "jump"
            %% or there is a run of values near the histogram peak.
            % ???
            ccc = data';
            ccc = ccc(:);
            ccc(abs(diff(ccc)) > 0.1) = NaN;
            ccd = [0; diff(ccc)];
            this_cor = zeros (size (dd));
            st_en = [];
            for i = size(rr,2):-1:1
                st = find (abs (jt - rr(1,i)) <= 3);
                if length (st) >= 1
                    [mismatch, m] = min (abs(sj(st) - jump));
                    if mismatch < 0.3 * jump
                        st = st(m);
                    else
                        st = [];
                    end
                end
                en = find (abs (jt - rr(2,i)) <= 3);
                if length (en) >= 1
                    [mismatch, m] = min (abs(sj(en) + jump));
                    if mismatch < 0.3 * jump
                        en = en(m);
                    else
                        en = [];
                    end
                end
                mid = find (any(isfinite (ccd(rr(1,i):rr(2,i)))));
                if length ([st', en', mid']) >= 2
                    %if isempty (st)
                        st = rr(1,i);
                    %else
                    %    st = jt(st);
                    %end
                    %if isempty (en)
                        en = rr(2,i);
                    %else
                    %    en = jt(en);
                    %end
                    st_en = [st_en; st, en];
                    d1 (st:en) = dd (st:en) - jump;
                    this_cor (st:en) = -jump;
                end
            end

            %dd(ddd>thresh) = max(0, dd(ddd>thresh) - jump);

            if isempty (st_en)
              continue;
            end

            if (jump < 0.4)
              continue;
            end

           [cwh, cwc, cwm, cwo] = corr_with_temp (this_cor, meta, state);
           sd  = similarity_of_durations (st_en);
           sot = similarity_of_on_times  (st_en, meta);
           vop = variance_of_on_power (st_en, d1);
           nd  = not_dinner (this_cor);
           [du, dus] = days_used (st_en, meta);
           du_norm  = (du / length (validDays)) / 6;
           dus_norm = (dus/ sum (isfinite (dd (1, meta.summer)))) / 3;
           mean_duration = mean (st_en(:,2) - st_en (:,1));

           % Try to find pool pump first
           factor(1) = -4 * (jump - 0.5) .^ 2;
           factor(2) = 2*(du_norm + dus_norm);
           factor(3) = -4 * max (du_norm - dus_norm, 0);
           factor(4) = min (du_norm - 0.3, 0);
           factor(5) = sd;
           factor(6) = min (0, mean_duration - 12) / 10;
           factor(7) = nd/8;
           factor(8) = -vop;
           factor(9) = cwm / 5;
           factor(10) = -cwh / 10;
           factor(11) = -cwc / 10;

                % > 0 for pool pump
           pool_ranking = sum (factor);

           if (pool_ranking > -0.2)  % check those that just miss
             type = pool_pump;
           else
               % Classify as Air con., heater, pool pump, other
               if (cwh > 0.7 && cwc < 0.3 && jump > 0.7)
                 type = aircon;
               elseif (du > 0.8 * length(validDays) && jump < 0.8 && jump > 0.1)
                 type = pool_pump;
               elseif (cwc > 0.7)
                 type = heater;
               else
                 type = other;
               end
           end

           %fprintf('Type is %s\n', type_names{type});

           % Identify days consistent with estimated type
           switch (type)
             case aircon
                % keep except days between 125 and 275
                toKeep = (  st_en(:,1) < 125 * meta.SamPerDay ...
                          | st_en(:,1) > 275 * meta.SamPerDay);
             case heater
                % keep except days before 75 and after 325
                toKeep = (  st_en(:,1) >  75 * meta.SamPerDay ...
                          & st_en(:,1) < 325 * meta.SamPerDay);
             case pool_pump
                % keep times at least 3 hours?  (Consider how often on?)
                toKeep = (st_en(:,2)-st_en(:,1) >= 6);
             case other
                toKeep = true (1, size (st_en, 1));
           end

           % if very few samples kept, probably noise
           % Should be more selective: if clustered, or start times clustered
           % more likely to be signal.
           if (sum (st_en(toKeep,2)) - sum (st_en(toKeep,1)) < 8)
             toKeep = toKeep & false;
           end

           this_cor = zeros (size (dd));
           discard_cor = this_cor;
           discarded = st_en(~toKeep,:);
           for i = 1:size (discarded,1)
               discard_cor (discarded(i,1):discarded(i,2)) = -jump;
           end
           discard_cor(end,end) = -jump;        % make scale consistent

           % Consider only consistent 
           st_en = st_en(toKeep,:);
           for i = 1:size (st_en,1)
               this_cor (st_en(i,1):st_en(i,2)) = -jump;
           end

           this_cor = max (this_cor, -dd);      % Remove negative bits
           % TODO: Interpolate turn-on / turn-off

           dd = dd + this_cor;

           dd_old = dd;
           cor (n,:) = this_cor(:);
        end
    end
    power = (type == pool_pump) * jump;
    pump = 6 * (type == pool_pump);       % TODO:  What is '6'??
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dummy(s, data)
    [sep, big] = max(s);          % most "separable" peak
    if sep > 1
        thresh = b(r(big,2));       % end of "rectangle" between peaks

        d = data;
        d(:) = (data > thresh);

        dd = data;
        dd(data>thresh) = NaN;

figure(2); imagesc(data');
figure(3); imagesc(d');
figure(4); imagesc(dd');
keyboard;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cwh, cwc, cwm, cwo] = corr_with_temp (this_cor, meta, state)
  used_hot  = sum (any (this_cor (:,meta.hotDays)));
  used_cold = sum (any (this_cor (:,meta.coldDays)));
  used_mild = sum (any (this_cor (:,state.mildDays)));

  other_summer = false (1,meta.Days);
  other_summer (meta.summer) = 1;
  other_summer (meta.hotDays) = 0;
  used_other_summer = sum (any (this_cor (:,other_summer)));

  other_winter = false (1,meta.Days);
  other_winter (meta.winter) = 1;
  other_winter (meta.coldDays) = 0;
  used_other_winter = sum (any (this_cor (:,other_winter)));

  other_days = true (1,meta.Days);
  other_days(meta.hotDays) = 0;
  other_days(meta.coldDays) = 0;
  other_days(state.mildDays) = 0;
  used_other= sum (any (this_cor (:, other_days)));

  all_bursts = used_hot + used_cold + used_mild + used_other;

  cwh = used_hot  / min (length (meta.hotDays),  all_bursts);
  cwh = cwh + (1-cwh) * (used_other_summer ...
                         / min (sum (other_summer), all_bursts));
  cwc = used_cold / min (length (meta.coldDays), all_bursts);
  cwc = cwc + (1-cwc) * (used_other_winter ...
                         / min (sum (other_winter), all_bursts));
  cwm = used_mild / all_bursts / (length (state.mildDays)  / meta.Days);
  cwo = used_other/ all_bursts / (sum (other_days) / meta.Days);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sd = similarity_of_durations (st_en)
    %[a, b] = hist(st_en(:,2) - st_en(:,1), 10);
    %figure(5); hist(st_en(:,2) - st_en(:,1), 10);
    
    len = size (st_en, 1);
    sd = 0.1 / (0.1 + var (st_en(:,2) - st_en(:,1)) / sqrt (len));
    if len <= 3
      sd = sd / (4-len);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sot = similarity_of_on_times (st_en, meta)
    sot = 1 / var (mod (st_en(:,1), meta.SamPerDay));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean variance of power while "on", weighted by the duration of the on period
function vop = variance_of_on_power (st_en, dd)
  c = arrayfun (@(s, e)(var (dd(s:e))), st_en(:,1), st_en(:,2));
  weights = st_en(:,2) - st_en(:,1);
  weights = weights / sum (weights);
  vop = sum (c .* weights);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Days used, and Days used in summer
function [du, dus] = days_used (st_en, meta)
    days = ceil (st_en(:,1) / meta.SamPerDay);
    mask = zeros (1, meta.Days);
    mask(days) = 1;
    du = sum (mask);
    dus = sum (mask (meta.summer));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nd = not_dinner (cor_this)
    profile = sum (cor_this, 2);
    % Should avoid computation for each call.
    %nd = sum (profile ([1:35, 47:48])) / sum (profile);
    nd = sum (profile ([1:35/48*end, end-1:end])) / sum (profile);
end
