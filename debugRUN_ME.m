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

tic;

bulkDataPath = '~/rsrch/power/UnitedEnergy/';
metaDataPath = './';

% Set to 1 to run on commercial and industrial customers, 0 for residential
DoCnI = 0;

% Set to a non-zero value to process only  MaxUsers customers
%  (Overwritten below if  DoLong = 0.)
MaxUsers = 0;
%MaxUsers = 88000;

% Which phases to execute
startAt = 1;
endAt = 2;

% Set to 0 to run from an alternate set of files, and output to "short_..."
% Set to 1 to run from the main file.
DoLong = 0;
    
if ~exist('OCTAVE_VERSION', 'builtin')	% if Matlab...
    batchRows = 1000000;		% number of lines to read at a time
    batchUsers = 3000;
    dbstop if error;
else				% else if Octave...
    batchRows = 130000;		% number of lines to read at a time
    batchUsers = 300;
    debug_on_error (1);
end

Year = 2013;

meta.non_pools = [];	%non_pools;
meta.real_pools = [];	%real_pools;

meta.HotDayThresh = 32;
meta.ColdDayThresh = 15;
daysIn = @(Y)(365 + (mod(Y,4) == 0 & (mod(Y,100)~=0 | mod(Y,400) == 0)));
meta.Days = daysIn(Year);
meta.SamPerDay = 48;
DayOffset = 735234 + sum(daysIn(2014:Year));
meta.peakhours = (6*meta.SamPerDay/24-1):(23*meta.SamPerDay/24);
meta.ph = length(meta.peakhours);

addpath('scripts');
addpath ('/usr/local/MATLAB/R2013b/toolbox/signal/signal')	% medfilt1

global iif;
iif = @(varargin) varargin{2+(varargin{1}==0)};

fid = fopen ([metaDataPath, 'ESD.csv'], 'r');
rows = textscan(fid, '%s %f %f %f %f %f %f %f\n', 'collectoutput', 1, 'delimiter', ',');
EnStandingDataNMIs = rows{1};
EnStandingData = [zeros(size(rows{2},1), 1), rows{2}];
fclose(fid);


% postcode groups
load ([metaDataPath, 'postcode_neighbours.txt']);
meta.postcode_neighbours = postcode_neighbours;
pclist = EnStandingData(:,6);
pclist(isnan(pclist)) = -1;
meta.pclist = unique(pclist);
invalid_postcodes = setdiff(pclist, postcode_neighbours(:,1));
if ~isempty(invalid_postcodes)
    fprintf('%d postcodes are invalid.  Assuming Melbourne''s demographic centre.\n', length(invalid_postcodes));
    [~,pos] = ismember(invalid_postcodes, pclist);
    EnStandingDataNMIs(pos);
    %pclist(pos) = 3146;
    %keyboard
end

fprintf('Loading temperatures\n');
temperatures = load ([metaDataPath, sprintf('temperatures%d.txt',Year)]);
meta.max_temp = max (temperatures,[],2);
min_temp      = min (temperatures,[],2);
mean_temp     = mean(temperatures,2);

DayOfWeek = mod(1:meta.Days,7);
% In 2013, 5 and 6 of Jan were Saturday and Sunday.
meta.wednesdays= find(DayOfWeek == mod(DayOffset+0,7));
meta.thursdays = find(DayOfWeek == mod(DayOffset+1,7));
meta.fridays   = find(DayOfWeek == mod(DayOffset+2,7));
meta.saturdays = find(DayOfWeek == mod(DayOffset+3,7));
meta.sundays   = find(DayOfWeek == mod(DayOffset+4,7));
meta.mondays   = find(DayOfWeek == mod(DayOffset+5,7));
meta.tuesdays  = find(DayOfWeek == mod(DayOffset+6,7));
meta.weekends  = sort([meta.saturdays, meta.sundays]);
meta.weekdays  = setdiff(1:meta.Days, meta.weekends);

meta.holidays = [];
meta.summer = [1:86, 331:meta.Days];
meta.autumn = 87:171;
meta.winter = 172:262;
meta.spring = 263:330;

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

meta.coldDays = find(meta.max_temp <= meta.ColdDayThresh);
meta.hotDays  = find(meta.max_temp >= meta.HotDayThresh);
%mildDays = find((meta.max_temp > 20) .* (meta.max_temp < 25));

t = 0.5:0.5:24;

% Load data in, if not already in memory

if DoLong
    filename = {[bulkDataPath, 'data_only.csv'],...
                [bulkDataPath, 'LabelledUnified.csv']};
    phfiles= {[bulkDataPath, 'phase1.mat'], ...
              [bulkDataPath, 'phase2.mat']};
    outPrefix = '';
else
    filename = {[bulkDataPath, 'solar_only.csv'],...
                [bulkDataPath, 'LabelledUnified.csv']};
%    filename = {[bulkDataPath, 'have_pools1.csv'],...
%                [bulkDataPath, 'LabelledUnified.csv']};
%    filename = {[bulkDataPath, 'solar_only.csv'],...
%        [bulkDataPath, 'LabelledUnified.csv']};
    phfiles= {[bulkDataPath, 'phase1_start10000.mat'], ...
              [bulkDataPath, 'phase2_start10000.mat']};
    outPrefix = 'short_';
    
    MaxUsers = 0;
MaxUsers = 1000;
end
if ~exist(filename{1}, 'file')	% Month 13 may not exist.
    error('Cannot find input file %s\n', filename{1});
end

phases = {@phase1, @phase2};

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
        A = [];
        NMIlist = {};
        %labelled = [];
        maxOfFirstFile = 0;
        fid = fopen(filename{i});
        %NMIs = zeros(1,batchUsers);
        NMIs = {}; NMIs{batchUsers} = '';
        n=0;
        data = nan(batchUsers, meta.Days, meta.SamPerDay);
        meta.hw_users    = [];
        meta.solar_users = [];
        meta.postcode  = [];
        while ~done
            fprintf('reading from %d (%d)... ', offset, length(NMIs));
            rows = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', batchRows, 'collectoutput', 1, 'delimiter', ',');
            fprintf('Done\n');
	    if ~isempty(rows{1})
		NMIlist = [NMIlist; rows{1}];
		A = [A; rows{2}];
	    end
            
            doneInner = false;
            users = unique(NMIlist);
            % (Include CnI and non-CnI)
            ptr = 1;
tic
            while ~doneInner
                usr = NMIlist{ptr};
                nextUsr = ptr + find(~strcmp(usr, NMIlist(ptr+1:end)), 1, 'first');
                % If partial user, read more before processing
                isLast = false;
                if isempty(nextUsr)
                    doneInner = true;
                    if ~feof(fid)
                        break
                    else
                        isLast = true;
                        done = true;
                        nextUsr = size(A,1)+1;
                    end
                end
%fprintf('usr %s nextUsr %d n %d\n', usr, nextUsr, n);
                
                % If this is the type of user we are interested in
                %if(usr(1)=='"')
                %    numeric = str2num(usr(2:end-1));
                %else
                %    numeric = str2num(usr);
                %end
                %idx = (EnStandingData(:,1) == numeric);
                idx = strcmp(EnStandingDataNMIs, usr);
                isResidential = (EnStandingData(idx, 3) == 5);
                if i == 1
                    skipUser = (isResidential ~= (DoCnI==0));
                else
                    [in, pos] = ismember(usr, state.NMIs);
                    skipUser = (in && pos < NumInFirstFile);
                end
                if isempty(skipUser) || ~skipUser
                    n = n+1;
                    NMIs{n} = usr;
                    %remap(usr) = n;

                    % Read solar/hw/postcode data from EnStandingData
                    if EnStandingData(idx, 5) == 1
                        meta.solar_users = [meta.solar_users, n];
                    end
                    if EnStandingData(idx, 2) == 2
                        meta.hw_users = [meta.hw_users, n];
                    end
                    if any(idx)
                        meta.postcode(n) = EnStandingData(idx,6);
                    else
                        meta.postcode(n) = 3146;    % demographic centre
                    end
                    
                    % Construct data structure
                    data(n,:,:) = nan(1,meta.Days,meta.SamPerDay);
                    for j = ptr:(nextUsr-1)
                        if isnan(data(n,A(j,1)-DayOffset,1))
                            data(n,A(j,1)-DayOffset,:) = A(j,(1:meta.SamPerDay)+1);
                        else
                            data(n,A(j,1)-DayOffset,:) = data(n,A(j,1)-DayOffset,:) + reshape(A(j,(1:meta.SamPerDay)+1), [1,1,meta.SamPerDay]);
                        end
                    end
                    
                    if offset+n == MaxUsers || (feof(fid) && isLast)
                        doneInner = true;
                        done = true;
                    end
                end
                
                ptr = nextUsr;
                if n == batchUsers || done == true
toc
                    % Free some memory for use in "phase" call
                    A = A(ptr:end,:);	% leave part-read user for next time
                    NMIlist = NMIlist(ptr:end,:);
                    ptr = 1;

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
                    fprintf('Starting phase %d NMI %s\n', ph, NMIs{1});
%if offset >= 66000
save ([bulkDataPath, NMIs{1}, "_dbg2"], 'data')
                    state = phase(state, meta, data);
%end                    
                    offset = offset + size(data,1);
                    
                    if ~done
                        NMIs = {}; NMIs{batchUsers} = '';
                        data = nan(batchUsers, meta.Days, meta.SamPerDay);
                        n=0;
                        meta.hw_users    = [];
                        meta.solar_users = [];
                        meta.postcode    = [];
                    end
tic
                end
            end
            A = A(ptr:end,:);		% leave part-read user for next time
            NMIlist = NMIlist(ptr:end,:);		% leave part-read user for next time
        end
        %maxOfFirstFile = max(state.NMIs);
        if UseLabelled && i == 1 && ph == 1
            NumInFirstFile = length(state.NMIs);
        end
    end
    %data = data(1:length(NMIs),:,:);
    
    state = phase(state, meta, []);		% perform cleanup
    save (phfiles{ph}, 'state');
end
