function [time, height, my_jp, mate_jp] = match_jump (x, my_row, jp)
% Given a time-series  x,  find the jump that matches the one at  my_row,
% whose nominal size is  jp (used to ensure the sign of the match is correct).
% If  my_row  is negative, match a jump of size  jp  and disregard the
% actual step at time |my_row|.

  % Force to column vector
  x = x (:);

%    % estimate time and size of jumps that may occur over two time slots
%    d = diff([x(end); x; x(1:2)]);
%    pr = d .* [d(2:end); d(1)];
%    eps = min(abs(pr))/2;
%    [~, mx] = max([[pr(end); pr(1:end-1)], eps*ones(length(d),1), pr], [], 2);
%    mx (pr == 0) = 2;
%    idx = (mx ~= 2);
%    idx(1)   = (1==0);
%    idx(end) = (1==0);
%
%    new_d = d;
%    times = [0:length(x), 1]';
%    new_d(idx) = d(idx) + d(find(idx) + mx(idx)-2);
%    new_d(new_d == 0) = realmin;
%    times(idx) = times(idx) + (mx(idx)-2) .* (1 - d(idx) ./ new_d(idx));
%
%    times = times(2:end-1);	% fractional time slot of estimated jump
%    d     = new_d(2:end-1);	% size of estimated jump

  % Internally, a "jump" is before the transition,
  % but other code uses the time after the transition.
  me = 1 + mod (abs (my_row)-2, length (x));

  firstEdge = -1;
  for i = 1:2         % Allow second attempt if first match is "poor"

    % estimate time and size of jumps that may occur over two time slots
    d = diff([x; x(1:2)]);
    idx = (sign(d(1:end-1)) == sign(d(2:end)));
    new_d = d(1:end-1);
    new_d(idx) = d(idx) + d(find(idx)+1);
    times = (1:length(x))';
    delta = 1 - (d(idx) ./ new_d(idx));
    delta (isnan(delta)) = 0;
    times(idx) = times(idx) + delta;

    top_of_step = (1:length(x))';
    top_of_step(new_d > 0) = top_of_step(new_d > 0) + 1;
    % if my row is in the middle of a bigger transition, consider whole
    if top_of_step(me) == top_of_step(abs (my_row)) ...
       && top_of_step(me) ~= abs (my_row)
      me = abs (my_row);
    end
    % Avoid overflow later
    top_of_step(idx & (new_d>0)) = top_of_step(idx & (new_d>0)) + 1;
    
    d = new_d;


    % Rotate to make "me" first, and copy it to the end
    x2 = [x(me:end); x(1:me)];
    d2 = [d(me:end); d(1:me)];
    times2 = [times(me:end); times(1:me)];
    top_of_step2 = [top_of_step(me:end); top_of_step(1:me)+length(top_of_step)] - me+1;
    top_of_step2(end-1:end) = min(top_of_step2(end-1:end), length(x));

  %    d2(1)
    if my_row < 0
      d2(1) = jp;
    elseif sign(d2(1)) ~= sign(jp)
      d2(1) = eps * sign(jp);
    end

    mismatch_edge = (0.1 + abs(d2 + d2(1)));

    if d2(1) > 0           % matching down jump after
      body_height = cummin([x2(2); x2(2:end)]);
    else
      body_height = cummin(x2(end:-1:1));
      body_height = body_height(end:-1:1);
    end

    match_body = body_height(top_of_step2);
    idx = (match_body > 0);
    match_body = (0.1 + match_body(idx));

    [~, edge] = min(mismatch_edge(idx) ./ match_body);

    if ~isempty (edge)              % re-map back to position in original
      edge = find (idx, edge);
      edge = edge(end);
    else
      edge = 2;   % TODO: find sensible default
    end

    % Customized for pool pumps: If this is "expected", only do first pass
    if i == 1
      height = body_height(top_of_step2(edge));
      if (d2(1) > 0 && mod(edge,length(x)) > 2 && mod(edge,length(x)) < 30)...
        || ...
         (d2(1) < 0 &&mod(-edge,length(x)) > 2 &&mod(-edge,length(x)) < 30)...
        || ...
         (height > 0.1 && height < 1)
           break
      else
        % If no good match for this edge, look for one in a
        % "smoothed" version, if that doesn't change this edge too much
        %    TODO: avoid taking "diff" of whole vector
        firstEdge = edge;
        new_x = -rolling_min(-rolling_min(x(:),5),5);
        new_tmp = diff ([new_x; new_x(1:2)]);
        tmp = diff ([x; x(1:2)]);
        if abs(new_tmp(me) - tmp(me)) > 0.5 * abs(tmp(me))
          %keyboard
          break
        else
          x = new_x;
        end
      end
    end
  end


%    [mismatch_edge(edge)  match_body(top_of_step2(edge)) match_body(edge)]
  height = body_height(top_of_step2(edge));
  my_jp  = d2(1);
  mate_jp= d2(edge);
  time   = times2(edge);
  time   = 1 + mod(time, length(x));	% convert to time  after  jump
%if my_row == 39  keyboard  end
%    keyboard
end


%function [time, height] = match_jump_old (my_jump, jp, times, x)
%    	% convert all to column vectors
%    jp    = jp   (:);
%    times = times(:);
%    x     = x    (:);
%
%    x2     = x([round(times(my_jump)):length(x), 1:round(times(my_jump))]);
%    jp2    = jp   ([my_jump:length(jp), 1:my_jump]);
%    times2 = times([my_jump:length(jp), 1:my_jump]);
%
%    mismatch_edge = (0.1 + abs(jp2 (1:end-1)+ jp2(1)));
%
%    if jp2(1) > 0           % matching down jump after
%        body_height = cummin([x2(2); x2(2:end)]);
%    else
%        body_height = cummin(x2(end:-1:1));
%        body_height = body_height(end:-1:1);
%    end
%    body_height = body_height(mod(round(times2-times2(1)), length(x))+1);
%    match_body = (0.1 + body_height(1:end-1));
%
%    [~, edge] = min(mismatch_edge ./ match_body);
%    height = body_height(edge);
%    time   = times2     (edge);
%end
