% flow:
% One pass of data, calculating:
% - Max solar to date per hour per postcode
% Global calculations
% - Load minus solar
% Second pass of data, calculating
% - Data sets to cluster
% - Initial centroids (based on many repetitions with first chunk)
% Global calculations
% - Output stats

%pragma warning on Octave:syntax-extension

tic;

% Melbourne
hemisphere = 'south';
location.latitude  = -37.813611;
location.longitude = 144.963056;
location.altitude  = 60;    % average altitude of 60m.  Won't matter much
time.UTC = +10;
time.min = 0;
time.sec = 0;
time.day = 21;

if 0>0
  if exist ('C:', 'dir')
    bulkDataPath = 'C:/Users/Lachlan Andrew/rsrch/NILM/PecanSt/';
  else
    bulkDataPath = '';
  end
  metaDataPath = bulkDataPath;
  Year = 2016;
  meta.SamPerDay = 24;
  meta.timestamp_uses_daylightsaving = false;
  meta.dataset_name = 'Austin';

  % Austin
  hemisphere = 'north';
  time.UTC = -6;
  location.latitude = +37;
  location.longitude = -87;
  location.altitude = 60;
elseif 1>0
  if exist ('C:', 'dir')
    bulkDataPath = 'G:/My Drive/rsrch/power/UnitedEnergy/';
    if ~exist (bulkDataPath, 'dir')
      bulkDataPath = 'C:/Users/seymo/Desktop/UnitedEnergy/';
    end
  else
    bulkDataPath = '~/UnitedEnergy/';
  end
  metaDataPath = './';
  Year = 2013;
  meta.SamPerDay = 48;
  meta.timestamp_uses_daylightsaving = true;
  meta.dataset_name = 'Melbourne';
else
  bulkDataPath = 'C:/Users/seymo/Desktop/UnitedEnergy/';
  metaDataPath = bulkDataPath;
  Year = 2014;
  meta.SamPerDay = 48;
  meta.timestamp_uses_daylightsaving = true;
end

if 0>0
    bulkDataPath = 'C:\Users\seymo\OneDrive\Documents\y5\FYP\Data\Swiss';
    metaDataPath = bulkDataPath;
    %Switzerland
    Year = 2015;
    hemisphere = 'north';
    time.UTC = 2;
    location.longitude = 7.519069;
    location.latitude = 47.50685;
    location.altitude = 400;
    meta.SamPerDay = 48;
    
end

time.year = Year;
meta.time = time;
meta.location = location;
meta.hemisphere = hemisphere;
meta.metaDataPath = metaDataPath;
meta.bulkDataPath = bulkDataPath;
meta.do_cooling = false;
meta.do_heating = false;

% Set to 1 to run on commercial and industrial customers, 0 for residential
DoCnI = 0;

% True if we only process users for which smart meter was present on day 1
full_years_only = true;

% Set to a non-zero value to process only  MaxUsers customers
MaxUsers = 100;
%MaxUsers = 300;
%MaxUsers = 1000;
%MaxUsers = 88000;

% Which phases to execute
startAt = 2;
endAt = 2;

% Set to 0 to run from an alternate set of files, and output to "short_..."
% Set to 1 to run from the main file.
DoLong = 0;

% Trust EnStandingData to tell us if customers are solar or HWS?
% If true (non-zero), only ESD customers are processed.
TrustESD = 1;

if exist ('OCTAVE_VERSION', 'builtin')
  batchUsers = 1500;
else
  batchUsers = 3000;
end
if MaxUsers ~= 0
  batchUsers = min (batchUsers, MaxUsers);
end
batchRows = 365 * batchUsers;		% number of lines to read at a time

dbstop if error;
%dbstop if naninf;        % currently crashes octave

meta.non_pools = [];	%non_pools;
meta.real_pools = [];	%real_pools;

meta.HotDayThresh = 32;
meta.ColdDayThresh = 15;
daysIn = @(Y)(365 + (mod(Y,4) == 0 & (mod(Y,100)~=0 | mod(Y,400) == 0)));
meta.Days = daysIn(Year);
DayOffset = 735234 + sum(daysIn(2014:Year));
meta.peakhours = (6*meta.SamPerDay/24-1):(23*meta.SamPerDay/24);
meta.ph = length(meta.peakhours);

addpath('scripts');
%addpath ('/usr/local/MATLAB/R2013b/toolbox/signal/signal')	% medfilt1
p = '/home/lachlan/.octave/signal-1.3.2/x86_64-pc-linux-gnu-api-v50+';
if exist (p, 'dir')
  addpath (p);
end

global iif;
iif = @(varargin) varargin{2+(varargin{1}==0)};

ESDfile = [metaDataPath, 'ESD.csv'];
fid = fopen (ESDfile, 'r');
if (fid == -1)
  fprintf ('Could not open file "%s" for reading\n', ESDfile);
  EnStandingDataNMIs = {};
  ESDhw = [];
  ESDresidential = [];
  pclist = [];
  TrustESD = false;
else
  rows = textscan(fid, '%s %u8 %u8 %*f32 %u8 %u16 %*f32 %*f32', 'delimiter', ',');
  fclose(fid);
  EnStandingDataNMIs = rows{1};
  ESDhw = rows{2};
  ESDresidential = rows{3};
  ESDsolar = rows{4};
  pclist = rows{5};
  clear rows
end

% postcode groups
try
  load ([metaDataPath, 'postcode_neighbours.txt']);
  meta.postcode_neighbours = postcode_neighbours;
  %pclist = EnStandingData(:,6);
  %pclist(isnan(pclist)) = -1;
  meta.pclist = unique(pclist);
  invalid_postcodes = setdiff(pclist, postcode_neighbours(:,1));
  if ~isempty(invalid_postcodes)
      fprintf('%d postcodes are invalid.  Assuming Melbourne''s demographic centre.\n', length(invalid_postcodes));
      [~,pos] = ismember(invalid_postcodes, pclist);
      EnStandingDataNMIs(pos);
      %pclist(pos) = 3146;
      %keyboard
  end
catch
  meta.pclist = 0;
  meta.postcode_neighbours = 1;
end

fprintf('Loading temperatures\n');
temperatures = load ([metaDataPath, sprintf('temperatures%d.txt',Year)]);
meta.temperatures = temperatures;
meta.max_temp = max (temperatures,[],2);
min_temp      = min (temperatures,[],2);
mean_temp     = mean(temperatures,2);

DayOfWeek = mod(1:meta.Days,7);
% In 2013, 5 and 6 of Jan were Saturday and Sunday.
meta.thursdays = find(DayOfWeek == mod(DayOffset+0,7));
meta.fridays   = find(DayOfWeek == mod(DayOffset+1,7));
meta.saturdays = find(DayOfWeek == mod(DayOffset+2,7));
meta.sundays   = find(DayOfWeek == mod(DayOffset+3,7));
meta.mondays   = find(DayOfWeek == mod(DayOffset+4,7));
meta.tuesdays  = find(DayOfWeek == mod(DayOffset+5,7));
meta.wednesdays= find(DayOfWeek == mod(DayOffset+6,7));
meta.weekends  = sort([meta.saturdays, meta.sundays]);
meta.weekdays  = setdiff(1:meta.Days, meta.weekends);

meta.holidays = [];

meta.January  = (1:31);
meta.February = (1:28+(meta.Days-365)) + max(meta.January);
meta.March    = (1:31) + max(meta.February);
meta.April    = (1:30) + max(meta.March);
meta.May      = (1:31) + max(meta.April);
meta.June     = (1:30) + max(meta.May);
meta.July     = (1:31) + max(meta.June);
meta.August   = (1:31) + max(meta.July);
meta.September= (1:30) + max(meta.August);
meta.October  = (1:31) + max(meta.September);
meta.November = (1:30) + max(meta.October);
meta.December = (1:31) + max(meta.November);

if isequal (hemisphere, 'south')
  meta.summer = [1:86, 331:meta.Days];
  meta.autumn = 87:171;
  meta.winter = 172:262;
  meta.spring = 263:330;
  %meta.daylight_saving(1) = min (meta.sundays(meta.sundays >= meta.April(1)));
  %meta.daylight_saving(2) = min (meta.sundays(meta.sundays >= meta.October(1)));
  meta.daylight_saving = (1:meta.Days) <  min (meta.sundays(meta.sundays >= meta.April(1))) ...
                       | (1:meta.Days) >= min (meta.sundays(meta.sundays >= meta.October(1)));
else
  meta.winter = [1:86, 331:meta.Days];
  meta.spring = 87:171;
  meta.summer = 172:262;
  meta.autumn = 263:330;
  %meta.daylight_saving(2) = min (meta.sundays(meta.sundays >= meta.April(1)));
  %meta.daylight_saving(1) = min (meta.sundays(meta.sundays >= meta.October(1)));
  meta.daylight_saving = (1:meta.Days) >= min (meta.sundays(meta.sundays >= meta.April(1))) ...
                       & (1:meta.Days) <  min (meta.sundays(meta.sundays >= meta.October(1)));
end


meta.coldDays = find(meta.max_temp <= meta.ColdDayThresh);
meta.hotDays  = find(meta.max_temp >= meta.HotDayThresh);
%mildDays = find((meta.max_temp > 20) .* (meta.max_temp < 25));

%t = 0.5:0.5:24;
t = (1:meta.SamPerDay) * 24/meta.SamPerDay;

% Load data in, if not already in memory

if DoLong ~= 0
    filename = {[bulkDataPath, 'data_only.csv'],...
                [bulkDataPath, 'LabelledUnified.csv']};
    phfiles= {[bulkDataPath, 'phase1.mat'], ...
              [bulkDataPath, 'phase2.mat'], ...
              [bulkDataPath, 'phase3.mat']};
    outPrefix = bulkDataPath;
else
    filename = {[bulkDataPath, 'solar_only.csv']};
%     filename = {[bulkDataPath, 'data_only.csv']};
     % filename = {[bulkDataPath, 'PecanSt_cooling.csv']};
     %filename = {[bulkDataPath, 'PecanSt_pool.csv']};
     %filename = {[bulkDataPath, 'PecanSt_solar.csv']};
%                [bulkDataPath, 'LabelledUnified.csv']};
%    filename = {[bulkDataPath, 'LabelledUnified.csv']};
%    filename = {[bulkDataPath, 'have_pools1.csv'],...
%                [bulkDataPath, 'LabelledUnified.csv']};
%        [bulkDataPath, 'LabelledUnified.csv']};
    phfiles= {[bulkDataPath, 'phase1_start3000.mat'], ...
              [bulkDataPath, 'phase2_start3000.mat'], ...
              [bulkDataPath, 'phase3_start3000.mat']};
    outPrefix = [bulkDataPath, 'short_'];
end
if ~exist(filename{1}, 'file')	% Month 13 may not exist.
    error('Cannot find input file %s\n', filename{1});
end

phases = {@phase1, @phase2, @phase3};
% To read the first tranche of data, and then exit:
%phases{1} = @stop_if_three_args; phases{2} = phases{1};

if DoCnI == 1
    meta.outCSV = [outPrefix, 'CnIusers.CSV'];
    meta.clusterCSV = [outPrefix, 'CnIclusters.CSV'];
    meta.out_solarCSV = [outPrefix, 'CnIsolar.CSV'];
    meta.out_poolCSV = [outPrefix, 'CnIpool.CSV'];
    meta.out_hwCSV = [outPrefix, 'CnIhw.CSV'];
    featureFile = [outPrefix, 'CnIfeatures.CSV'];
else
    meta.outCSV = [outPrefix, 'resUsers.CSV'];
    meta.clusterCSV = [outPrefix, 'resClusters.CSV'];
    meta.out_solarCSV = [outPrefix, 'resSolar.CSV'];
    meta.out_poolCSV = [outPrefix, 'resPool.CSV'];
    meta.out_hwCSV = [outPrefix, 'resHW.CSV'];
    featureFile = [outPrefix, 'resFeatures.CSV'];
end

% State is saved at the end of each phase.  If the previous phases have
% been debugged and work, load their output instead of regenerating.
%endAt = length(phases);
state = [];
% Find latest saved state before startAt, and start from
% following phase
fprintf('Loading state...');
for i = startAt-1:-1:1
    try
        load(phfiles{i});
        break;
    catch
        startAt = startAt-1;
        state = [];
    end
end
fprintf(' Done.  Starting in phase %d\n', startAt);

% Run remaining phases
for ph = startAt:endAt
    phase = phases{ph};
    state = phase(state, meta);	% initialize
    offset = 0;
    meta.batchNumber = 0;
    UseLabelled = 0;		% Set to 1 to use file with demographic labels
    for i = 1:(1+UseLabelled)
        done = false;
        A = single ([]);
        NMIlist = {};
        %labelled = [];
        maxOfFirstFile = 0;
        fid = fopen(filename{i});
        if fid == -1
          fprintf ('Could not open file "%s" for reading\n', filename{i});
          return
        end
        %NMIs = zeros(1,batchUsers);
        NMIs = {}; NMIs{batchUsers} = '';
        n=0;
        data = nan(batchUsers, meta.Days, meta.SamPerDay, 'single');
        meta.hw_users    = [];
        meta.solar_users = [];
        meta.postcode  = [];
        format = ['%s', repmat(' %f32', [1, meta.SamPerDay+1])];
        while ~done
            fprintf('reading %d from %d ... ', length(NMIs), offset);
            %rows = textscan(fid, '%s %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32', batchRows, 'collectoutput', 1, 'delimiter', ',');
            rows = textscan(fid, format, batchRows, ...
                            'collectoutput', 1, 'delimiter', ',');
            fprintf('Done\n');
            if ~isempty(rows{1})
              NMIlist = [NMIlist; rows{1}];
              A = [A; rows{2}];
            end
            clear rows

            % Correct any days for which the meter is inverted
            % (Note: one meter recorded 0.006kWh solar generation at 10pm.)
            range = 1 + [1:9*meta.SamPerDay/48, 40*meta.SamPerDay/48:meta.SamPerDay];
            inverted = any (A(:, range) <  -0.01, 2) ...
                     & all (A(:, range) <= 0, 2);
            r = ranges(find (inverted));
            if ~isempty (r)
              A(inverted,2:end) = -A(inverted,2:end);
              % TODO: identify day on which meter was fixed
              %       flip only those hours before the fix.
              % TODO: Check out meter 1526 with continuous 0.45kW gen.
            end

            doneInner = false;
            users = unique(NMIlist);
            % (Include CnI and non-CnI)
            ptr = 1;

            while ~doneInner
                ends = find (~strcmp (NMIlist(1:end-1), NMIlist(2:end)));
                if feof(fid)
                    ends = [ends(:); length(NMIlist)]';
                end
                usrs = NMIlist(ends);

                ESD_idx = zeros(length(usrs),1);
                usr_in_ESD = false (length(usrs),1);
                for j = 1:length(usrs)
                   [E, u] = find (strcmp (EnStandingDataNMIs, usrs{j}));
                   usr_in_ESD(j) = ~isempty (u);
                   if usr_in_ESD(j)
                       ESD_idx(j) = E;
                   end
                end     % for
                ESD_usrs = ESD_idx(usr_in_ESD);

                %% *** zeros -> only do users in EnStandingData
                %%     ones  -> do users not in EnStandingData (bug below)
                if TrustESD ~= 0
                    doUser = false (length(usrs), 1);
                else
                    doUser = true (length(usrs), 1);
                end
                if i == 1
                    % doUser <- (not in EnStandingData) || (matching (DoCnI))
                    doUser(usr_in_ESD) = ((ESDresidential(ESD_usrs)==5)...
                                          == (DoCnI==0));
                else
                    fprintf ('Probably buggy but only used for demographics\n');
                    [in, pos] = ismember(usrs, state.NMIs);
                    doUser = ~(in & pos < NumInFirstFile);
                end     % if

                % Allow us to discard the users with
                % smart meters not installed at the start of the year.
                if full_years_only
                  starts = [1; ends(:)+1]';
                  if starts(end) > length (NMIlist)
                    starts = starts (1:end-1);
                  end
                  start_days = A(starts, 1);
                  doUser(start_days > DayOffset + 1) = false;

                  % Also delete those with many missing days later
                  days_per_entry = diff (starts);
                  doUser(days_per_entry < 0.9 * meta.Days) = false;
                end

                isLast = false;     % last textscan in file?
                last = find (doUser, batchUsers - n);
                if length(last) < batchUsers-n
                                          % Add all these to data, then continue
                   doneInner = true;
                   last = find (doUser, 1, 'last');
                   if feof(fid)
                       isLast = true;
                       done = true;
                       next_ptr = size(A,1)+1;
                   elseif isempty (ends)
                       break;           % no complete user, so next textscan
                   else
                       next_ptr = ends(end)+1;  % first usr to read next pass
                   end
                else                % Add to data, call phase, then continue
                   last = last(end);
                   next_ptr = ends(last)+1;
                end
                count = sum(doUser(1:last));

                if (MaxUsers > 0 && offset+n+count >= MaxUsers)
                    doneInner = true;
                    done = true;
                    count = MaxUsers - (offset+n);
                    last = find (doUser, count);
                    last = last(end);
                elseif (isLast && feof(fid))
                    doneInner = true;
                    done = true;
                end

                [NMIs{n+(1:count)}] = NMIlist{ends(doUser(1:last))};

        % usr_in_ESD:    logical (in usrs) of users in EnStandingData
        % doUser:        logical (in usrs) of users we will process eventually
        % doUser(1:last) logical (in usrs) of users we will process now
        % doUser_idx:    index into 1:count of usr_idx
        % ESD_idx(i):    index into ESD of 'i'th user within usrs, or zero
        % last:          position in doUser of last processed now
        % count:         number processed now
                                  %%% Assumes doUser initialized to zero
                new_postcodes  = 3146 * ones (last, 1);
                doUser_idx       = usr_in_ESD (1:last) & doUser(1:last);
                doUser_ESD_idx   = ESD_idx(doUser_idx);
                new_postcodes(doUser_idx) = pclist(doUser_ESD_idx);
                if TrustESD ~= 0
                  new_solar_usrs   = n+find(ESDsolar(doUser_ESD_idx) == 1)';
                  new_hw_usrs      = n+find(ESDhw(doUser_ESD_idx) == 2)';
                  meta.hw_users    = [meta.hw_users,    new_hw_usrs];
                  meta.solar_users = [meta.solar_users, new_solar_usrs];
                else
                  meta.hw_users    = NaN;
                  meta.solar_users = NaN;
                end

                meta.postcode(n+(1:count)) = new_postcodes(doUser(1:last));

                doUser = [find(doUser); 2]; % 2, as dummy '1' for last iteration
                data(n+(1:count),:,:) = nan(count,meta.Days,meta.SamPerDay);
                if ~doUser(1)
                  ptr = ends(doUser(1) - 1) + 1;
                end
                for k = 1:count
                    n = n+1;
                    days = A(ptr:ends(doUser(k)), 1) - DayOffset;
                    dr = ranges (days);
                    if size (dr, 2) == 1 || length(unique(days)) == length(days)
                        data(n,days,:) = A(ptr:ends(doUser(k)), ...
                                           (1:meta.SamPerDay)+1);
                    else
                        % slow path.  Could be optimized, but only occurs for
                        % if days have multiple records, such as for solar or HW
                        fprintf ('-');
                        for j = ptr:ends(doUser(k))
                            if isnan(data(n,A(j,1)-DayOffset,1))
                                data(n,A(j,1)-DayOffset,:) = ...
                                               A(j,(1:meta.SamPerDay)+1);
                            else
                                data(n,A(j,1)-DayOffset,:) ...
                                    = data(n,A(j,1) -DayOffset,:) ...
                                      + reshape(A(j,(1:meta.SamPerDay)+1), ...
                                                [1,1,meta.SamPerDay]);
                            end
                        end
                        fprintf ('>');
                    end
                    ptr = ends(doUser(k + 1) - 1) + 1;
                end             % for k
                ptr = next_ptr;


                if n == batchUsers || done == true
                    % Free some memory for use in "phase" call
                    A = A(ptr:end,:);	% leave part-read user for next time
                    NMIlist = NMIlist(ptr:end,:);
                    ptr = 1;            % prevent re-truncation after the loop

                    % Set up meta, then call phase to process these users
                    NMIs = NMIs(1:n);
                    meta.NMIs = NMIs;
                    meta.UserOffset = offset;
                    if i == 1
                        %[yesno, meta.solar_users] = ismember(meta.has_solar, NMIs);
                        %meta.solar_users = meta.solar_users(yesno)';

                        %[yesno, meta.hw_users] = ismember(meta.has_hws, NMIs);
                        %meta.hw_users = meta.hw_users(yesno)';
                    else
                        meta.solar_users = NaN;
                        meta.hw_users = [];
                    end
                    meta.batchNumber = meta.batchNumber + 1;
                    % trim extra NaNs
                    data = data(1:length(NMIs),:,:);
                    meta.phase = ph;
                    fprintf('Starting phase %d: %d from NMI %s\n', ...
                            ph, size(data,1), NMIs{1});
                    state = phase(state, meta, data);
                    offset = offset + size(data,1);

                    if ~done
                        clear NMIs;
                        NMIs{batchUsers} = '';
                        data = nan(batchUsers, meta.Days, meta.SamPerDay);
                        n=0;
                        meta.hw_users    = [];
                        meta.solar_users = [];
                        meta.postcode    = [];
                    end % if
                end     % if
            end         % while ~doneInner
            A = A(ptr:end,:);		  % leave part-read user for next time
            NMIlist = NMIlist(ptr:end,:); % leave part-read user for next time
        end
        fprintf('\n');

%            while ~doneInner
%                usr = NMIlist{ptr};
%                nextUsr = ptr + find(~strcmp(usr, NMIlist(ptr+1:end)), 1, 'first');
%                % If partial user, read more before processing
%                isLast = false;
%                if isempty(nextUsr)
%                    doneInner = true;
%                    if ~feof(fid)
%                        break
%                    else
%                        isLast = true;
%                        done = true;
%                        nextUsr = size(A,1)+1;
%                    end
%                end
%%fprintf('usr %s nextUsr %d n %d\n', usr, nextUsr, n);
%
%                idx = strcmp(EnStandingDataNMIs, usr);
%                isResidential = (EnStandingData(idx, 3) == 5);
%                if i == 1
%                    skipUser = (isResidential ~= (DoCnI==0));
%                else
%                    [in, pos] = ismember(usr, state.NMIs);
%                    skipUser = (in && pos < NumInFirstFile);
%                end
%                if isempty(skipUser) || ~skipUser
%                    n = n+1;
%                    NMIs{n} = usr;
%                    %remap(usr) = n;
%
%                    % Read solar/hw/postcode data from EnStandingData
%                    if EnStandingData(idx, 5) == 1
%                        meta.solar_users = [meta.solar_users, n];
%                    end
%                    if EnStandingData(idx, 2) == 2
%                        meta.hw_users = [meta.hw_users, n];
%                    end
%                    if any(idx)
%                        meta.postcode(n) = EnStandingData(idx,6);
%                    else
%                        meta.postcode(n) = 3146;    % demographic centre
%                    end
%
%                    % Construct data structure
%                    data(n,:,:) = nan(1,meta.Days,meta.SamPerDay);
%                    days = A(ptr:nextUsr-1, 1) - DayOffset;
%                    dr = ranges (days);
%                    if size (dr, 2) == 1 || length(unique(days)) = length(days)
%                        data(n,days,:)=A(ptr:(nextUsr-1), (1:meta.SamPerDay)+1);
%                    else
%                        % slow path.  Could be optimized, but only occurs for
%                        % if days have multiple records, such as for solar or HW
%                        for j = ptr:(nextUsr-1)
%                            if isnan(data(n,A(j,1)-DayOffset,1))
%                                data(n,A(j,1)-DayOffset,:) = ...
%                                               A(j,(1:meta.SamPerDay)+1);
%                            else
%                                data(n,A(j,1)-DayOffset,:) ...
%                                    = data(n,A(j,1) -DayOffset,:) ...
%                                      + reshape(A(j,(1:meta.SamPerDay)+1), ...
%                                                [1,1,meta.SamPerDay]);
%                            end
%                        end
%                    end
%
%                    if offset+n == MaxUsers || (feof(fid) && isLast)
%                        doneInner = true;
%                        done = true;
%                    end
%                end
%
%                ptr = nextUsr;
%                if n == batchUsers || done == true
%                    % Free some memory for use in "phase" call
%                    A = A(ptr:end,:);	% leave part-read user for next time
%                    NMIlist = NMIlist(ptr:end,:);
%                    ptr = 1;
%
%                    % Set up meta, then call phase to process these users
%                    NMIs = NMIs(1:n);
%                    meta.NMIs = NMIs;
%                    meta.UserOffset = offset;
%                    if i == 1
%                        %[yesno, meta.solar_users] = ismember(meta.has_solar, NMIs);
%                        %meta.solar_users = meta.solar_users(yesno)';
%
%                        %[yesno, meta.hw_users] = ismember(meta.has_hws, NMIs);
%                        %meta.hw_users = meta.hw_users(yesno)';
%                    else
%                        meta.solar_users = NaN;
%                        meta.hw_users = [];
%                    end
%                    meta.batchNumber = meta.batchNumber + 1;
%                    % trim extra NaNs
%                    data = data(1:length(NMIs),:,:);
%                    fprintf('Starting phase %d NMI %s\n', ph, NMIs{1});
%%if offset >= 66000
%                    state = phase(state, meta, data);
%%end
%                    offset = offset + size(data,1);
%
%                    if ~done
%                        NMIs = {}; NMIs{batchUsers} = '';
%                        data = nan(batchUsers, meta.Days, meta.SamPerDay);
%                        n=0;
%                        meta.hw_users    = [];
%                        meta.solar_users = [];
%                        meta.postcode    = [];
%                    end
%                end
%            end
%            A = A(ptr:end,:);		% leave part-read user for next time
%            NMIlist = NMIlist(ptr:end,:);		% leave part-read user for next time
%        end
        %maxOfFirstFile = max(state.NMIs);
        if UseLabelled && i == 1 && ph == 1
            NumInFirstFile = length(state.NMIs);
        end
    end
    %data = data(1:length(NMIs),:,:);

    state = phase(state, meta, []);		% perform cleanup
    try
      save (phfiles{ph}, 'state');
    catch
      save (phfiles{ph}, '-v7.3', 'state');
    end
end

function arg1 = stop_if_three_args(arg1, arg2, arg3)
  if nargin == 3
    filename = evalin('base', 'filename');
    filename = filename{1};
    prefix = [evalin('base', 'bulkDataPath'), 'PecanSt_'];
    if length (filename) > length (prefix) ...
        && isequal (filename(1:length (prefix)), prefix) ...
        && isequal (filename(end-3:end), '.csv')
      var_name = filename(length (prefix)+1:end-4);
      eval ([var_name ' = arg3;']);
    else
      disp (arg3);
      keyboard
    end

    ground_NMIs = evalin('base', 'NMIs');

    save ([evalin('base', 'bulkDataPath'), var_name], var_name, 'ground_NMIs');
    keyboard
  end
end