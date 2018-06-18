function state = phase3 (s, meta, data)

% Provides:
%    (See variables in "isempty(s)")
% Needs:
%    meta.Days
%    meta.peakhours

%pragma warning on Octave:syntax-extension

global iif;
global fout_s;

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
    
    state.solarUsed = sparse(1,1);
    %%%% temporaries for this phase
    
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

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif isempty(data)            % Tidy up at end
    % Cluster by aircon power
    
    % Cluster by heating power

    fclose(fout_s);
    fclose(state.fout_hw);
    fclose(state.fout_pool);
    fclose(state.fout_solar);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else                            % Main process, called once per pass
    ind = meta.UserOffset + (1:size(data,1));
    
    % Calculate solar correction
    [su, cor_solar] = calc_cor_solar (data, s, meta, ind);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate pool pumps
    fprintf('Detecting pool pumps...\n');
    pump = zeros(1,size(data,1));
    records = zeros(length(ind),4);
    cor_pool = sparse(size(data,1),size(data(1,:),2));
    for i = 1:size(data,1)
        corrected = squeeze(data(i,:,:));
        corrected(:) = corrected(:)' + s.cor_hw(ind(i),:) + cor_solar(i,:);
        vamp = cumsum(s.vampires(ind(i),:));
        try
          [p, cor, power] = has_pool_pump(corrected, vamp, meta, s.NMIs{ind(i)}, any(su == i));
          pump(i) = power;
          if p>0
              cor_pool(i,:) = cor(:)';  % Slow?  Build cor_pool at end?
              fprintf('%d:%g\t', ind(i), power);
          end
        catch
          vv = cumsum(s.vampires(ind(i),:));
          save (['bug_',state.NMIs{ind(i)}], 'corrected', 'vv', 'meta', 'ind', 'i', 'su');
          fprintf ('Error detecting pool pump for user %s', state.NMIs{ind(i)});
          keyboard;
        end

%        [power, cor2] = has_untimed_pump (corrected + cor, vamp, meta, ind(i), any(su == i), state);
    end
    state.poolPump(ind) = pump;
%for i = find(pump(:))'
%  figure(1); imagesc (squeeze(data(i,:,:))');
%  figure(2); imagesc (reshape (cor_pool(i,:), [365, 48])');
%  keyboard
%end
    fprintf('\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Estimate airconditioner peak power, and power used on cooling
    fprintf('cooling... ');
    
    % Power consumption on hot days ("cooling") and hot non-vacations ("cooling_v")
    [i2, i1] = meshgrid(1:meta.SamPerDay, meta.hotDays);
    s2i = sub2ind ([meta.Days meta.SamPerDay], i1(:), i2(:));
    hot_no_hw_solar = reshape(data(:,s2i) + s.cor_hw(ind,s2i) + cor_solar(:,s2i), [size(data,1), length(s2i)/meta.SamPerDay, meta.SamPerDay]);
    cooling    = squeezeNaN(@mean, hot_no_hw_solar, 2);
    cooling_v  = squeezeNaN(@mean, hot_no_hw_solar + repmat(full(s.vacation(ind,meta.hotDays)),[1,1,meta.SamPerDay]), 2);
    clear hot_no_hw_solar;
    
    % Subtract power consumption on mild days ("cooling") and mild non-vacations ("cooling_v"), to estimate the actual cooling power
    [i2m, i1m] = meshgrid(1:meta.SamPerDay, s.mildDays);
    s2im = sub2ind ([meta.Days meta.SamPerDay], i1m(:), i2m(:));
    mild_no_hw_solar = reshape(data(:,s2im) + s.cor_hw(ind,s2im) + cor_solar(:,s2im), [size(data,1), length(s2im)/meta.SamPerDay, meta.SamPerDay]);
    mild   = squeezeNaN(@mean, mild_no_hw_solar, 2);
    mild_v = squeezeNaN(@mean, mild_no_hw_solar + repmat(full(s.vacation(ind,s.mildDays)),[1,1,meta.SamPerDay]), 2);
    clear mild_no_hw_solar;
    
    %state.hot(ind,:) = cooling;
    cooling   = cooling   - mild;
    cooling_v = cooling_v - mild_v;
    
        fprintf('Done\n');
    airconPower_(:,1) = meta.SamPerDay/24*squeezeNaN(@max,cooling_v,2);
    % Avoid undersized array if for loop skips users
    airconPower_2 = zeros(size(airconPower_));
    
    % Find steps
    % use loop to avoid creating huge array
    % Ignore 11pm spikes from hot water
    %meta.peakhours = (4*meta.SamPerDay/24-1):(23*meta.SamPerDay/24);
    ph = length(meta.peakhours);
    for i = 1:size(data,1)
        d = squeeze(data(i,:,:));
        d(:) = d(:) + cor_solar(i,:)';
        userData = 2*(squeeze(d(meta.hotDays,meta.peakhours))'+ones(ph,1)*s.vacation(i,meta.hotDays));
        userData = min(userData(1:end-1), userData(2:end));
        if all(isnan(userData))
            continue
        end
        userData = [userData userData(end)];    % make up to full size
        
        [dev, use, alt, conf, jumps, bins] = bigDevices(userData);
        airconPower_ (i,3) = alt(1,1);
        airconPower_2(i,3) = alt(1,2);
        airconPower_ (i,2) = alt(2,1);
        airconPower_2(i,2) = alt(2,2);
        binsj = bins{1};
        binsl = bins{2};
        
        %fprintf('Estimated AC power of %s: %f L %f(%f) J %f(%f)\n', s.NMIs{i}, airconPower_(i,1), airconPower_(i,2), airconPower_2(i,2), airconPower_(i,3), airconPower_2(i,3));
%keyboard
        
    end         % for
    
    state.airconPower(ind) = squeezeNaN(@mmedian,(airconPower_+airconPower_2),2);
    state.mildDay(ind,:) = mild;
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Estimate power used on heating
    fprintf('heating... ');
    
    % Power consumption on cold days ("heating") and hot non-vacations ("heating_v")
    [i2, i1] = meshgrid(1:meta.SamPerDay, meta.coldDays);
    s2i = sub2ind ([meta.Days meta.SamPerDay], i1(:), i2(:));
    cold_no_hw_solar = reshape(data(:,s2i) + s.cor_hw(ind,s2i)...
                                           + cor_solar(:,s2i),...
                               [size(data,1), length(s2i)/meta.SamPerDay,...
                                meta.SamPerDay]);
    heating    = squeezeNaN(@mean, cold_no_hw_solar, 2);
    heating_v  = squeezeNaN(@mean, cold_no_hw_solar ...
                 + repmat(full(s.vacation(ind,meta.coldDays)),...
                          [1,1,meta.SamPerDay]), 2);
    clear cold_no_hw_solar;
    
    % Subtract power consumption on mild days ("heating") and mild non-vacations ("heating_v"), to estimate the actual heating power
    %state.cold (ind,:) = reshape(heating, [size(heating,1), meta.SamPerDay]);
    heating   = heating - mild;
    heating_v = heating_v - mild_v;
    
    fprintf('Done\n');
    heatPower_(:,1) = 2*squeezeNaN(@max,heating_v,2);   % '2*' to convert to kW
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % estimate stay-at-home
    days = intersect(meta.weekends,s.mildDays);
    [i2, i1] = meshgrid(1:meta.SamPerDay, days);
    s2i = sub2ind ([meta.Days meta.SamPerDay], i1(:), i2(:));
    WEmildPatterns = reshape(data(:,s2i) + s.cor_hw(ind,s2i) + cor_solar(:,s2i), [size(data,1), length(s2i)/meta.SamPerDay, meta.SamPerDay]);
    WEmildPatterns = squeezeNaN(@mean, WEmildPatterns,2);
    
    days = intersect(meta.weekdays,s.mildDays);
    [i2, i1] = meshgrid(1:meta.SamPerDay, days);
    s2i = sub2ind ([meta.Days meta.SamPerDay], i1(:), i2(:));
    WDmildPatterns = reshape(data(:,s2i) + s.cor_hw(ind,s2i) + cor_solar(:,s2i), [size(data,1), length(s2i)/meta.SamPerDay, meta.SamPerDay]);
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
    bb = max(WDamPeak ./ business, 1);
    
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
    state.solarUsed(ind(slr)) = sum(cor_solar(slr,:) + min(0,data(slr,:)), 2);


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
    
    fprintf ('Done\n');
end
 
% vi: set expandtab
