% [jumps times, signed_jumps] = find_jumps(userData)
% [jumps times, signed_jumps] = find_jumps(userData, offset, n)
% Finds jumps in the underlying continuous time signal approximated by
% sums over intervals
% Output jumps is a vector of the absolute sizes of the jumps
% Output times is a vector of the corresponding locations within userData
% Output signed_jumps is a vector of the upward jumps
% Increasing input offset increases the minimum size of jump detected
% Input n is the number of time steps before/after the jump for which the
% signal should be roughly constant.
% Considers the fact that a jump part way through an interval will cause
% small increases in two consecutive intervals
function [jumps times signed_jumps] = find_jumps(userData, offset, n)
    if nargin < 3
        n = 2;
        if nargin < 2
            offset = 0;
        end
    end
    diffs = diff(userData(:));
    diffs(~isfinite(diffs)) = 0;
    ab = cumsum(abs(diffs));

    % declare a "jump" if change this time step is as big as the sum of
    % the changes in the  n  previous and  n  subsequent time steps.
    diffs1 = abs(diffs);
    times1 = find(diffs1((n+1):(end-n)) > offset + 0.5*(ab((2*n+1):end)-ab(1:end-(2*n)))) + n;
    jumps1 = diffs1(times1);
    signed_jumps1 = diffs(times1);

    % declare a "jump" if the change over *two* adjacent time steps is as big
    % as the change in the  n  previous  and  n  subsequent time steps
    diffs2 = abs(diffs(1:end-1) + diffs(2:end));
    times2 = find(diffs2(2:end-1) >= diffs2(3:end) ...
               &  diffs2(2:end-1) >= diffs2(1:end-2)) + 1;
    times2 = setdiff (times2, times2 + 1);	% remove double counting
    before = ab(max(1, times2-n));
    after = ab(min(max(times2), times2+n+1));
    t = find(diffs2(times2) > offset + 0.5*(after-before));
    times2 = times2(t);
    jumps2 = diffs2(times2);
    signed_jumps2 = (diffs(times2) + diffs(times2+1));

        % Give "fractional" time for partial jumps,
        % but retain integers for old calls with fewer than 3 arguments
    if nargin >= 3
        times2 = times2 + (diffs(times2+1)./(diffs(times2)+diffs(times2+1)));
    end

%figure(1);
%hold off;
%lenHotDays = 36;
%plot(reshape(userData, [lenHotDays, 21]), 'g');
%hold on;
%
%t2 = kron(times2, [1 0 0]');
%t2(2:3:end) = t2(1:3:end)+2;
%t2(3:3:end) = 1;
%j2 = userData(t2);
%j2(3:3:end) = NaN;
%    % remove wrapped lines
%t2 = mod(t2, lenHotDays); toskip = 3*find(t2(2:3:end) < 2)-1; t2(toskip) = NaN;
%plot(t2, j2, 'b');
%
%t1 = kron(times1, [1 0 0]');
%t1(2:3:end) = t1(1:3:end)+1;
%t1(3:3:end) = 1;
%j1 = userData(t1);
%j1(3:3:end) = NaN;
%t1 = mod(t1, lenHotDays);  t1(t1 == 0) = NaN;	% remove wrapped lines
%plot(t1, j1, 'r');
%
%figure(2)
%hist(jumps1);
%
%figure(3)
%hist(jumps2);

%keyboard

    % Merge lists.
    % Eliminate entries in jumps1 that are one of the two steps of jumps2
    [~, t] = setdiff(times1, [floor(times2); floor(times2)+1]);
    [times idx] = sort([times2; times1(t)]);
    jumps = [jumps2; jumps1(t)];
    jumps = jumps(idx);
    signed_jumps = [signed_jumps2; signed_jumps1(t)];
    signed_jumps = signed_jumps(idx);
end
