function state = phase1 (s, meta, data)
%pragma warning on Octave:syntax-extension

% Provides:
%    (See variables in "isempty(s)")
% Needs:
%    meta.Days
%    meta.peakhours
%    meta. ...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2                  % Set up data structures
    %%%% data to pass to later phases
    state.NMIs = {};
    state.gridLoad      = zeros(meta.Days,1);
    state.gridLoadAll   = zeros(meta.Days,meta.SamPerDay);
    state.mildDays = uint16 ([]);

    state.daysValid = uint16 ([]);

    % Hot water
    state.hws_power   = sparse ([1;2]);     % force to column => sparser
    state.hws_energy  = sparse ([1;2]);
    state.has_hws_est = int32  ([]);
    state.hw_users    = int32  ([]);
    state.cor_hw      = sparse ([]);

    % Solar
    state.solar_power     = single ([]);
    %state.solar_users_est = int32  ([]);
    state.solar_users     = int32  ([]);
    state.solar_cap       = single ([]);
    state.solar_az        = single ([]);
    state.solar_ze        = single ([]);
    state.cor_solar       = sparse ([]);

    % Vacations
    state.vacation  = sparse ([]);
    state.emptyDays = uint16 ([]);

    % always-on devices.
    % Users-by-meta.Days array of *differences*, so cumsum(state.vampires,2)
    state.vampires = sparse(1,meta.Days);

    % Maximum use
    state.peak   = [];
    state.peak39 = [];  % Peak between 3pm and 9pm

    state.solar_start =  7*meta.SamPerDay/24;
    state.solar_end   = 17*meta.SamPerDay/24;
    state.dark_start  = 20*meta.SamPerDay/24;
    state.dark_end    =  5*meta.SamPerDay/24;

    %%%% temporaries for this phase
    state.gridLoadCount = zeros(meta.Days,1);
    state.gridLoadCountAll = zeros(meta.Days,meta.SamPerDay);
    %state.location;            % Used for solar.  real value set below
    %state.time;                % Used for solar.  real value set below

    %capFactor       = zeros(length(s.solar_users),meta.Days,state.dark_start-state.dark_end+1);
    state.capFactor = single ([]);

 %{
    location.latitude  = -37.813611;
    location.longitude = 144.963056;
    location.altitude  = 60;    % average altitude of 60m.  Won't matter much
    time.UTC = +10;
    time.year = 2013;
 %}
    location = meta.location;
    time = meta.time;
    time.min = 0;
    time.sec = 0;
    time.day = 21;

    fprintf('Calculating sun positions\n');
    % Calculate sun's angle half hourly at summer/winter solstice.
    clear SunPos_summer SunPos_winter;
    % Preallocate output structure array
    sunPos_summer.az((state.solar_end-state.solar_start)+1) = 0;
    sunPos_winter.az((state.solar_end-state.solar_start)+1) = 0;
    for i = (state.solar_end-state.solar_start):-1:0
        time.hour = (state.solar_start+i) * 24/meta.SamPerDay;
        time.month = 12; SunPos_summer (i+1) = sun_position(time, location);
        time.month = 6;  SunPos_winter (i+1) = sun_position(time, location);
        %       time.month = 9;  SunPos_equinox(i+1) = sun_position(time, location);
    end

    state.location = location;
    state.time = time;
    state.SunPos_summer = SunPos_summer;
    state.SunPos_winter = SunPos_winter;

    % Pre-compute location of the sun at each half hour of each day.
    time.month = 1;
    %for i = 1:meta.Days
    i = 1:meta.Days;
        time.day = i;           % Can handle "100th day of month 1" etc.
        for j = 0:(state.dark_start-state.dark_end)
            time.hour = (state.dark_end+j) * 24/meta.SamPerDay;

            % Record data to estimate generation at this time, assuming no cloud
            SunPos = sun_position_vector(time, location, true);
            state.s1(i,j+1) = cosd([SunPos.zenith]);
            state.s2(i,j+1) = sind([SunPos.zenith]);
            state.pp(i,j+1) = SunPos.azimuth;
            if (j >= state.solar_start - state.dark_end) && (j + state.dark_end <= state.solar_end)
                state.SunPos(i,j-state.solar_start+state.dark_end+1) = SunPos;
            end
        end
    %end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif isempty(data)            % Conclude phase 1
    s.gridLoad = s.gridLoad ./ s.gridLoadCount;
    s.gridLoadAll = s.gridLoadAll ./ s.gridLoadCountAll;
    [~, mildDays] = sort (s.gridLoad);
    s.mildDays = mildDays(1:150);

    %%%%%%%%
    % Calculate solar corrections
    % Ensure full size, even if isempty(u) for large k below.

    s = find_solar_by_pc (s, meta);

    state=rmfield(s, {'gridLoadCount', 'gridLoadCountAll'});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else 
    % Main process, called once per pass
    %s = find_solar_by_pc (s, meta);
    ind = meta.UserOffset + (1:size(data,1));
    state = s;
    state.NMIs = [s.NMIs, meta.NMIs];
    fprintf('Calculating Daily Use\n');
    state.daysValid(ind) = sum(any(~isnan(data),3),2);
    if length(meta.postcode) == length(ind) % mismatch for "labelled" data
        state.postcode(ind) = meta.postcode;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Identify and remove off-peak hot water usage
    %hws_power = 2*squeeze(skipNaN(@mean,max(data(:,:,47),data(:,:,48)) - data(:,:,46),2));
    fprintf('Identifying hot water parameters\n');
    % Identify off-peak HWS users by the presence of a "big" jump of a nearly
    % constant size at 11pm.  We don't know the size of the jump in advance,
    % and it will not occur on days that the HWS is turned off.
    % This misses some types of HWS, but detects some users that appear to
    % have HWS and are not flagged as such by  has_hws.
    hws_thresh = 0.5;
    hws_power = zeros(size(data,1),1);
    for i = 1:size(data,1)
        mx = max(data(i,:,end-1:end), [], 3);
        % big increase, starting When off peak starts
        fff = find((mx - data(i,:,end-2) > hws_thresh/2) ...
            .* (data(i,:,end-1) > data(i,:,end-2)));
        % big decrease starting when off peak starts
        ffg = find(mx - data(i,:,end-2) < -hws_thresh/2);
        % big increase starting just before off peak starts
        ffh = find(data(i,:,end-2) - data(i,:,end-3) > hws_thresh/2);
        m = skipNaN(@mean,mx(:,fff,:) - data(i,fff,end-2));
        v = skipNaN(@vvar,mx(:,fff,:) - data(i,fff,end-2),2);
        % (Customers identified as HWS, but not listed as such in  has_hws)
        %if ismember(i, [ 266   564   661   951]) keyboard end
        if length(fff) < 10 || m < hws_thresh || sqrt(v) > m || length(ffg) + length(ffh) > 0.75 * length(fff)
            hws_power(i) = 0;
            continue
        end

        ff = (data(i,fff,end) > data(i,fff,end-1));
        ffg = ff & (data(i,min(fff+1,size(data,2)),1) > 0.1 * data(i,fff,end))...
                 & (data(i,fff, end-2) < 1);
        if sum(ffg) > 5         % power at 48 more reliable if not turned off
                                % during that interval.  Use them if enough
            hp = 2*skipNaN(@mean,data(i,fff(ffg),end)-data(i,fff(ffg),end-2));
        elseif sum(ff) > 5      % power at 48 more reliable than max(47,48)
                                % so use it alone if there are enough
            hp = 2*skipNaN(@mean,data(i,fff(ff),end)-data(i,fff(ff),end-2));
        else
            hp = 2*skipNaN(@mean,mx(:,fff,:)-data(i,fff,end-2));
        end
        if ~isempty(hp)
            hws_power(i) = hp;
        end
    end
    state.hws_power = [s.hws_power; hws_power];
    has_hws_est = find(hws_power > hws_thresh);
    if (isfinite (meta.hw_users))
        hws_users = union (has_hws_est, meta.hw_users);
        state.hw_users  = [s.hw_users,  meta.hw_users+meta.UserOffset];
    else
        hws_users = has_hws_est';
        state.hw_users  = [s.hw_users,  hws_users+meta.UserOffset];
    end

    if ~isempty(has_hws_est)
        % Group off-peak times into contiguous blocks, instead of split over 2 days
        offpeak = setdiff(1:meta.SamPerDay, meta.peakhours(5:(end-1)));
        oph = length(offpeak);
        op = reshape(permute(data(has_hws_est,:,offpeak), [1 3 2]), ...
            [length(has_hws_est), meta.Days*oph]);

        % shift morning off-peak on 1 Jan to 31 Dec
        % to shift offpeak times to be contiguous
        offpeak_end = find(offpeak < meta.SamPerDay/2, 1, 'last');
        op = op(:, [offpeak_end+1:end 1:offpeak_end]);
        op = permute(reshape(op, [length(has_hws_est), oph, meta.Days]), [1 3 2]);
        op = op(:,1:meta.Days-1,:);     % ignore hybrid 1 Jan / 31 Dec day

        % turns off in middle of 30-min interval, so look at diff over 2 steps
        jumps2 = 2*(op(:,:,3:end)-op(:,:,1:(end-2)));
        % Ignore downward jumps in times after the power has dropped below
        % the estimated HWS power anyway
        cm = cummin(2*op(:,:,3:end-2),3);
        cm = cat(3,repmat(hws_power(has_hws_est),[1 meta.Days-1 2]), cm);
        jumps2 = jumps2 .* (bsxfun(@gt, cm, hws_power(has_hws_est)*0.8));

        % find what looks like the HWS turning off
        [~, turnsOff] = min(abs(bsxfun(@plus,jumps2, hws_power(has_hws_est))),[],3);
        % 3-D matrix, 1 except 0 in the first (turnsOff) els of each (i,j,:)
        afterTurnsOff = repmat(turnsOff, [1 1 oph]) < repmat(reshape(1:oph, [1 1 oph]), [length(has_hws_est), meta.Days-1,1]);
        turnOffInd = xor(afterTurnsOff(:,:,2:end-1), afterTurnsOff(:,:,1:(end-2)));

        %% The "turn-off" transition should be one or two downward jumps
        %% If hw_power_est is too small, sometimes the best match is an upward jump
        %% and a downward jump.  To correct those cases, take the largest two-step
        %% jump within +/- 1 of  turnsOff.
        %toi_left  = cat(3,turnOffInd(:,:,2:end, zeros(size(turnOffInd(:,:,1))));
        %toi_right = cat(3,zeros(size(turnOffInd(:,:,1))), turnOffInd(:,:,1:end-1));
        %[val deltaTurnsOff] = min(cat(4,jumps2.*toi_left, jumps2.*turnOffInd, jumps2.*toi_right));


        % Estimate power from downward jumps.
        jumps2 = reshape(jumps2, [length(has_hws_est) * (meta.Days-1), oph-2]);
        jumps2(turnsOff(:) < 3,:) = NaN;
        jumps2 = reshape(jumps2, [length(has_hws_est), (meta.Days-1), oph-2]);
        down_hws_power = -skipNaN(@mean, sum(jumps2 .*turnOffInd,3), 2);

        % Some HWS have an hourly "boost" of a fixed energy after turning off
        % These are in every second 30-min interval, starting at a point that
        % depends on when the main service turned off.

        % Look for hourly "boosts" after turn off.
        jumps = 2*(op(:,:,2:end)-op(:,:,1:(end-1)));
        peaks0 = (abs(op(:,:,1:(end-2)) - op(:,:,3:end)) < 2*jumps(:,:,1:end-1) ...
            .* (jumps(:,:,1:(end-1)) > 0.2)...
            .* (jumps(:,:,2:end) < 0)...
            .* afterTurnsOff(:,:,2:end-1)      );
        peaks0(:,:,2:end) = peaks0(:,:,2:end) & ~peaks0(:,:,1:(end-1));
        boost0 = sum(skipNaN(@sum,jumps(:,:,1:(end-1)) .* peaks0,2),3);
        boost0 = boost0 ./ sum(skipNaN(@sum,peaks0,2),3);

        % try again without out-liers
        peaks = peaks0 .* (bsxfun(@gt, jumps(:,:,1:(end-1)), 0.5*boost0))...
            .* (bsxfun(@lt, jumps(:,:,1:(end-1)),   2*boost0));
        boost = sum(skipNaN(@sum,max(jumps(:,:,1:(end-1)),-jumps(:,:,2:end)) .* peaks,2),3);
        boost = boost ./ sum(skipNaN(@sum,peaks,2),3);
        boost = boost .* (boost > 0.4);

        clear peaks0 boost0 turnOffInd;

        % Estimate heating in each interval
        % Normally, hws_power until turns off, except edge effects
        hw_usage = bsxfun(@times, hws_power(has_hws_est), ~afterTurnsOff(:,:,2:end));
        hw_usage = cat(3,zeros(size(hw_usage,1),meta.Days-1), hw_usage);
        for i = 1:size(has_hws_est,1)
            for j = 1:meta.Days-1        % don't use logical as j in args 2 & 3
                if turnsOff(i,j) > 1
                    %  TODO:  Handle cases where the two steps aren't both down
                    hw_usage(i,j,turnsOff(i,j)+1) = 2*(op(i,j,turnsOff(i,j)+1) - (op(i,j,turnsOff(i,j))-hws_power(has_hws_est(i))/2 + op(i,j,turnsOff(i,j)+2))/2);
                    %if hw_usage(i,j,turnsOff(i,j)) < 0 keyboard end
                end
            end
            % If power increases in two steps, assume HWS turned on part way through
            % sample, and so interpolate in proportion to observed usage
            ff = find(data(has_hws_est(i),1:meta.Days-1,end) > data(has_hws_est(i),1:meta.Days-1,end-1));
            hw_usage(i,ff,3) = min(hw_usage(i,ff,3), 2*data(has_hws_est(i),ff,end));
            hw_usage(i,ff,2) = 2*op(i,ff,2) - (2*op(i,ff,1) + 2*op(i,ff,3)-hw_usage(i,ff,3))/2;

            % clip to plausible range
            hw_usage(i,:,:) = min(hw_usage(i,:,:), hws_power(has_hws_est(i)));
            %hw_usage(i,:,1:3)   = min(hw_usage(i,:,1:3),   2*data(has_hws_est(i),ff,46:48));
            %hw_usage(i,:,4:end) = min(hw_usage(i,:,4:end), 2*data(has_hws_est(i),ff+1,1:oph-3));


            ffc = setdiff(1:meta.Days-1,ff);
            % Assume increase at the start is due to heating, but decrease is noise,
            % and increase beyond  hws_power  is also noise
            hw_usage(i,ffc,2)=min(hws_power(has_hws_est(i)),max(0,2*(op(i,ffc,2)-op(i,ffc,1))));

            % Model "boosts" that occur later
            peaks = squeeze(...
                (abs(op(i,:,1:(end-3)) - op(i,:,4:end)) < 0.2*boost(i)) ...
                .* (abs(jumps(i,:,1:(end-2)) - jumps(i,:,3:end) - boost(i)) < 0.2*boost(i)) ...
                .* (jumps(i,:,1:(end-2)) > 0)...
                .* (jumps(i,:,3:end) < 0.2*boost(i))...
                .* afterTurnsOff(i,:,2:end-2)...
                );
            peaks(:,2:end) = peaks(:,2:end) & ~peaks(:,1:(end-1));
            for j = 1:meta.Days-1
                %if j == d keyboard end
                %ptc = find(abs(peaks(1,j,:) - boost(i)) < 0.2*boost(i));
                ptc = find(peaks(j,:));
                portions = jumps(i,j,ptc) ./ (jumps(i,j,ptc) -jumps(i,j,ptc+2));
                %if portions < -0.1 || portions > 1.1 keyboard end
                portions = min(1,max(0,portions));
                hw_usage(i,j,ptc+1) = hw_usage(i,j,ptc+1) +    portions  * boost(i);
                hw_usage(i,j,ptc+2) = hw_usage(i,j,ptc+2) + (1-portions) * boost(i);
            end
        end

        clear afterTurnsOff jumps2;

        %if 0
        %    for i = 1:length(has_hws_est)
        %       validDays = find(all(~isnan(op(i,:,:)),3));
        %       figure(95); plot(squeeze(hw_usage(i,validDays,:))');
        %       figure(96); plot(squeeze(2*op(i,validDays,:) - hw_usage(i,validDays,:))');
        %       figure(97); hist(2*op(i,:),30)
        %       figure(98); plot (turnsOff(i,:))
        %       figure(99); plot (squeeze(2*op(i,validDays,:))')
        %       sleep(1);
        %       fprintf('User %d (%d of %d) has power %f (%f) and boost of %f (was %f %f)\n', has_hws_est(i), i, length(has_hws_est), hws_power(has_hws_est(i)), down_hws_power(i), boost(i), boost0(i), old_boost(i));
        %       %keyboard
        %    end
        %end

        % Convert back to a matrix of hot water usage per user / day / slot
        % suitable for
        % usage(has_hws_est,:,offpeak) = data(has_hws_est,:,offpeak) - hw_hours;"
        hw_hours = reshape(permute(hw_usage(:,:,:), [1 3 2]), ...
            [length(has_hws_est), (meta.Days-1)*oph]);
        hw_hours = [zeros(length(has_hws_est),offpeak_end), hw_hours,...
            zeros(length(has_hws_est),oph-offpeak_end)];
        hw_hours = permute(reshape(hw_hours, [length(has_hws_est), oph, meta.Days]), [1 3 2]) / 2;

        % Because of shifting, only 1:Days-1 were calculated.  Take ave for rest
        % Should ignore days for which data is NaN.
        hr = 0:(oph-offpeak_end+1);
        hw_hours(:,meta.Days,end-hr) = ...
            min(sum(hw_hours   (:, 1:meta.Days-1,end-hr),2) ...
            ./ max(1,sum(~isnan(data(has_hws_est,1:meta.Days-1,end-hr)),2)), ...
            data(has_hws_est,meta.Days,end-hr));
        hr = 1:offpeak_end;
        hw_hours(:,1,hr) = ...
            min(sum(hw_hours   (:, 1:meta.Days-1,hr),2) ...
            ./ max(1,sum(~isnan(data(has_hws_est,1:meta.Days-1,hr)),2)), ...
            data(has_hws_est,1,hr));

        tmp = zeros(size(has_hws_est,1),meta.Days,meta.SamPerDay);
        tmp(:,:,offpeak) = hw_hours;
        cor_hw = sparse(size(data,1), meta.Days*meta.SamPerDay);
        cor_hw(has_hws_est,:) = -tmp(:,:);
        cor_hw(has_hws_est,:) = min(0,max(cor_hw(has_hws_est,:), ...
                                    sparse(double(-data(has_hws_est,:)))));
        clear tmp hw_usage op turnsOff;
    else
        cor_hw = sparse(size(data,1), meta.Days*meta.SamPerDay);
    end

    state.cor_hw = [s.cor_hw; cor_hw];
    state.has_hws_est = [s.has_hws_est; meta.UserOffset+has_hws_est];
    state.hws_energy(meta.UserOffset+(1:size(data,1))) = -sum(cor_hw,2)/meta.Days;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solar capacities
    %if ~exist('solar_capacity', 'var')

    solar_users_est = (any(data(:,:) < 0, 2));
    if isnan(meta.solar_users)
        %solar_users = reshape(solar_users_est,[1,length(solar_users_est)]);
        solar_users = find( solar_users_est(:)');
        non_solar   = find(~solar_users_est(:)');
            
    elseif isempty(meta.solar_users)
        %solar_users = reshape(solar_users_est,[1,length(solar_users_est)]);
        solar_users = find( solar_users_est(:)');
        non_solar   = find(~solar_users_est(:)');
    else
        solar_users = meta.solar_users;
        non_solar = 1:size(data,1);
        non_solar(solar_users) = 0;
        non_solar = non_solar(non_solar ~= 0);
    end

    %solar_users_est = find(solar_users_est);
    %state.solar_users_est= [s.solar_users_est; solar_users_est+meta.UserOffset];
    state.solar_users = [s.solar_users, solar_users+meta.UserOffset];

    if ~isempty(solar_users)
	[cap, az, ze, seen, capFactor, daily_min, vampires] ...
			    = fit_all_solar (data(solar_users,:,:), meta, s);

	% range (within list of solar users) corresponding to *this* data batch
	solar_user_range = length(state.solar_ze)+(1:length(ze));

  state.solar_cap(solar_user_range) = cap;
  state.solar_az (solar_user_range) = az;
  state.solar_ze (solar_user_range) = ze;
  state.seen     (solar_user_range,:,:) = seen;
	state.seen     (solar_user_range,:,:) = seen;
	state.capFactor(solar_user_range,:,:) = capFactor;
	state.daily_min(solar_user_range,:) = daily_min;
                                                    % no sparse single type
        state.vampires (solar_users+meta.UserOffset,:) = double(vampires);
    end
    % Estimate vampires for non-solar users
    if length(solar_users) < length(ind)
	vampires = find_vampires(reshape(permute(data(non_solar,:,:), [2,1,3]), length(non_solar)*size(data,2), size(data,3)));
	vampires = reshape(vampires, [size(data,2), length(non_solar)])';
	state.vampires(non_solar+meta.UserOffset,:) ...  % no sparse single type
                       = double([vampires(:,1), diff(vampires,1,2)]);
    end




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set to NaN days that customers seem to be away
    % (currently lowest peak power days, with usage not correlated with
    % higher-power days.)  The main contribution is identifying how many low-power
    % days.  In future, should consider:
    %   peak power
    %   average power
    %   correlation with mild days
    %   jumps in correlation near jumps in power
    %   clustering of vacation days

    vac_customers = [];
    vac_days = [];
    % Find "step" in daily peak usage
    non_hws = setdiff(1:size(data,1),hws_users);
    if ~isempty(non_hws)
        fprintf('Calculating max daily for non-HWS users\n');
	% Don't need squeeze as trailing singleton removed anyway
	% Don't need skipNaN as each day is either all NaN or none.
        %maxDaily(non_hws,:) = squeeze(skipNaN(@mmax,data(non_hws,:,:),3));
        maxDaily(non_hws,:) = max(data(non_hws,:,:), [], 3);
    end
    % Daily max, excluding HWS peaks
    if ~isempty(hws_users)
        fprintf('Calculating max daily for HWS users\n');
        %maxDaily(hws_users,:) = squeeze(skipNaN(@mmax,data(hws_users,:,meta.peakhours),3));
        maxDaily(hws_users,:) = max(data(hws_users,:,meta.peakhours), [], 3);
    end
    maxDaily(isinf(maxDaily)) = Inf;    % get rid of -Inf

    fprintf('Identifying vacations\n');

    [sortDaily, sdi] = sort(maxDaily,2);
    sortDaily(~isfinite(sortDaily)) = NaN;

    p = nan(size(data,1),1);

    if 1
        offset = 2;     % +2 to reduce weight of tiny counts
        for i = 1:size(sortDaily,1)
            if exist('true_p', 'var') && ismember(true_p(i) , 0:meta.Days)
                continue;
            end

            userData = squeeze(data(i,:,:));
            validDays = find(isfinite(sortDaily(i,:)),1,'last');
            if isempty(validDays)
                continue
            end
            bins = max(2*sqrt(validDays), 6);
            done = 0;
            prev_min = meta.Days;
            new_pi = 0;
            while ~done
                noVacation = 0;
                h = hist(sortDaily(i,:),bins);
                r(1) = (offset+h(2)) / h(1);            j(1) = h(1);
                r(2) = (offset+h(3)) / max(h(1), h(2)); j(2) = h(1)+h(2);
                [a, b] = min(r);
                if a < 0.5 && j(b) > validDays / 50
                    colour = 'b';
                    % hack:  If first bar of histogram is very flat, use it
                    if b == 2 && h(1) > h(2) && ...
                            sortDaily(i,h(1)+1) - sortDaily(i,h(1)) > sortDaily(i,h(1)) - sortDaily(i,1)
                        b = 1;
                        colour = 'c';
                    end
                    % end hack
                    ppos = b;
                    %done = 1;
                else
                    [~, mid] = max(h(2:end));
                    % Give extra weight to "early" jumps,
                    % as vacations usually short
                    cm = cummax(h(1:mid));
                    rm = h(1:mid);              % "recent max"
                    for j = 1:(bins/20)
                        rm(j+1:end) = max(rm(j+1:end), h(1:mid-j));
                    end
                    %rm(3:end) = max(rm(3:end), h(1:mid-2));
                    weights = validDays ./ (validDays + cumsum(h(1:mid)));
                    [val,  pos]  = max(h(2:1+mid) ./ (cm+offset) .* weights);
                    [val1, pos1] = max ((rm+offset) ./ (h(1:mid)+offset) .* weights);
                    if val < val1
                        ppos = pos1;
                        colour = 'r';
                    else
                        ppos = pos;
                        colour = 'g';
                    end
                    if max(val,val1) < (1.2)^max(ppos,2) && h(1) < mean([h(2),h(3)])
                        noVacation = 1; % if fairly uniform, no vacation
                    end
                end
                % look for two consecutive zeros (a steep jump)
                zpos = max(ppos, ceil(bins/7));
                twinzeros = find((h(1:zpos-1) == 0) .* (h(2:zpos)==0), 1);
                if ~isempty(twinzeros)
                    p(i) = sum(h(1:twinzeros));
                    colour = 'm';
                    done = 1;
                    % hack:  If twin zero masks an early big jump up, use big jump
                    if p(i) < 5 && val > 4 && sum(h(1:pos)) < 20
                        p(i) = sum(h(1:pos));
                        colour = 'y';
                    end
                    % end hack
                elseif noVacation
                    p(i) = 0;
                    colour = 'w';
                    done = 1;
                else
                    new_pi = sum(h(1:ppos));
                    if new_pi < 20 || abs(new_pi - prev_min) < 0.2 * prev_min
                        done = 1;
                        p(i) = min(new_pi, prev_min);
                    else
                        prev_min = min(new_pi, prev_min);
                        bins = ceil(bins * 1.2);
                        %fprintf ('%d %c\t', new_pi, colour);
                    end
                end
                % If the start is very flat, we should detect *some* vacation,
                % but may miss it if there are too few bins
                if done && p(i) == 0 && sortDaily(i,6) - sortDaily(i,1) < 0.01 && bins < 100
                    done = 0;
                    bins = ceil(bins * 1.2);
                    %fprintf ('repeat\t');
                end
            end

            if ismember(i,hws_users)
                hr = meta.peakhours;
            else
                hr = 1:meta.SamPerDay;
            end

            % If not really vacation, usage probably correlated with 'normal' use
            % Reduce estimated number of vacation days until correlation < 0.3
            % Problem: If holiday house has fridge, its daily pattern is 'normal'
            %      and the vacant house has high correlation with that.
            hist_pi = p(i);
            if p(i)
                % Check for fridge in otherwise idle house
                % That causes higher-than-normal correlation
                if p(i) >= validDays/2;
                    st1 = p(i);
                    en1 = validDays;
                    st2 = 1;
                    en2 = validDays - p(i) + 1;
                else
                    st1 = 1+p(i);
                    en1 = 1+2*p(i);
                    st2 = 1;
                    en2 = 1+p(i);
                end
                [x, y]=deal(skipNaN(@mean,userData(sdi(i,st1  :2:en1),hr),1),...
                    skipNaN(@mean,userData(sdi(i,st1+1:2:en1),hr),1));
                busyCorr=corr (x',y');
                % Take same number of samples as above, so no extra smoothing
                [x, y]=deal(skipNaN(@mean,userData(sdi(i,st2  :2:en2),hr),1),...
                    skipNaN(@mean,userData(sdi(i,st2+1:2:en2),hr),1));
                idleCorr=corr (x',y');
                %fprintf ('busyCorr %f  idleCorr %f\n', busyCorr, idleCorr);

                % Unless most correlation due to devices on in low-power times...
                if busyCorr > 0.8*idleCorr || p(i) < 5
                    usual = skipNaN(@mean,userData(sdi(i,p(i):end),hr),1);
                    x = skipNaN(@mean,userData(sdi(i,1:p(i)),hr),1);
                    cc = corr (x', usual') * (1-min(x)/max(x));
                    while p(i) && cc > 0.12
                        %fprintf ('corr: %f for %d\n', cc, p(i));
                        p(i) = floor(p(i)/2);
                        if p(i) > validDays/2
                            usual = skipNaN(@mean,userData(sdi(i,p(i):end),hr),1);
                        end
                        if p(i)
                            x = skipNaN(@mean,userData(sdi(i,1:p(i)),hr),1);
                            cc = corr (x', usual') * (1-min(x)/max(x));
                        end
                    end
                end
            end
        end
    end
    f = find(p>1);              % those that seem to have *some* holiday

    for i = f(:)'
        vac_customers = [vac_customers, repmat(i,1,p(i))];
        vac_days = [vac_days, sdi(i,1:p(i))];
    end
    vacation = sparse(vac_customers, vac_days, NaN);
    if size(vacation,1) < size(data,1)
        vacation(size(data,1),1) = 0;
    end
    if size(vacation,2) < meta.Days
        vacation(1, meta.Days) = 0;
    end
    clear vac_customers vac_days;

    state.vacation = sparse([s.vacation; vacation]);
    state.emptyDays = [s.emptyDays; full(sum(isnan(vacation),2))];

    data = data + repmat(full(vacation), [1 1 meta.SamPerDay]);

    % take "grid load" over non-solar, non-hw users, if any.
    % Otherwise, take over all users
    other_users = setdiff(setdiff(1:size(data,1), hws_users), solar_users);
    if isempty(other_users)
        other_users = 1:size(data,1);
    end

    state.gridLoad = s.gridLoad + skipNaN(@sum, skipNaN(@sum, data(other_users,:,:), 3), 1)';
    % Exclude days that are entirely NaN
    state.gridLoadCount = s.gridLoadCount + sum((skipNaN(@sum, data(other_users,:,:), 3)>0),1)';
    
    a = size (data);
    state.gridLoadAll = s.gridLoadAll + reshape (skipNaN(@sum, data(other_users,:,:), 1), a(2), a(3));
    state.gridLoadCountAll = s.gridLoadCountAll + reshape (sum(~isnan (data(other_users,:,:)), 1), a(2), a(3));

end
