function [morn, eve] = interp_edges (data)
% [morn, eve] = INTERP_EDGE (data)
%
% Input data is a matrix with two "1" entries in each column
% indicating the discrete locations of edges
% (transitions from shade to light or light to shade).
%
% Outputs morn and eve are vectors of width equal to data's
% containing interpolated times of the edges,
% linearly interpolated between changes in the discrete time.

  idx = (sum (data) == 2);
  if ~all (idx)
      fprintf ('Not two edges on every day.  Guessing edges...\n')
  end

  [present, morn] = max (data);
  morn = resolve (data, morn);
  morn = interp (morn);
  morn (~logical (present)) = NaN;

  d = data(end:-1:1, :);
  [present, eve] = max (d);
  eve = resolve (d, eve);
  eve = interp (size (d, 1) + 1 - eve);
  morn (~logical (present)) = NaN;
end

function times = interp (slots)
    changes = find ((diff (slots) ~= 0));
    segment_slots = [slots(1), slots(changes), slots(end)];
    directions = 0.5 * sign (diff (segment_slots));
    times = zeros (1, length (slots));
    %disp (changes)
    %disp (segment_slots)
    %disp (directions)
    for i = 1:length (changes) - 1
        pre_time  = slots(changes(i)+1) - directions(i+1);
        post_time = slots(changes(i)+1) + directions(i+2);
        %disp ([i, slots(i), pre_time, post_time])
        len = changes(i+1) - changes(i);
        if len == 1
            times (changes(i)+1) = (pre_time + post_time)/2;
        else
            times (changes(i)+1:changes(i)+len) ...
                = pre_time + (post_time - pre_time) * (1:len)/(len+1);
        end
    end
    times(1:changes(1)) = slots(1) + directions(2);
    times(changes(end)+1:end) = slots(end) - directions(end);
end


function edge = resolve (data, edge)
% First approximation to choosing better edge when there are multiple
% Currently just greedily minimizes distance to previous point
% Ideal would be shortest path, but too expensive
% Perhaps better would be 
    points = sum (data);

    choices = find (points > 1);
    if isempty (choices)
        return
    end

    if choices(1) == 1
        choices = choices (2:end);
    end

    for i = choices
        a = find (data(:, i));
        [~, pos] = min (abs (a - edge(i-1)));
        edge(i) = a(pos);
    end
end
