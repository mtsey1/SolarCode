function show_rectangles (rectangles, cv, days_to_hide)
% SHOW_RECTANGLES (RECTANGLES, CV)
% SHOW_RECTANGLES (RECTANGLES, CV, DAYS_TO_HIDE)
  tmp = cv;
  samPerDay = size (cv, 1);
  show_missed = false;

  if nargin == 2
    days_to_hide = [];
  elseif isequal (days_to_hide, 'missed')
    show_missed = true;
  end
  for i = days_to_hide
    tmp(:,i:7:end) = 0;
  end

  % Display power 'heatmap', preserving any zoom
  xx = get (gca, 'XLim');
  yy = get (gca, 'YLim');
  imagesc (tmp);
  if all (xx ~= 0) && all (yy ~= 0)
    set (gca, 'XLim', xx);
    set (gca, 'YLim', yy);
  end

  on_off = [rectangles(:).on_off];
  if (any (on_off(:) ~= round (on_off(:))))
    for r = 1:length (rectangles)
      if rectangles(r).on_off(1) ~= round (rectangles(r).on_off(1))
        %fprintf ('Error: rectangles(%d).on_off(1) = %g\n', r, rectangles(r).on_off(1));
        rectangles(r).on_off([1,3]) = round (rectangles(r).on_off([1,3]));
      end
    end
  end

  if ~isfield (rectangles, 'alt_hrs')
    draw_plain_rectangles (rectangles, samPerDay);
  else
    idx = cellfun (@length, {rectangles.alt_hrs}) < 2;
    draw_plain_rectangles (rectangles(idx), samPerDay);
    for r = find (~idx)
      draw_rectangle_DoW (rectangles(r), samPerDay);
    end
  end

  show_missed = false;
  if show_missed && isfield(rectangles, 'missed')
    oo = [rectangles(:).on_off];
    for r = 1:length (rectangles)
      m = oo(1,r) + ranges (find (rectangles(r).missed)) - 1;
      if ~isempty (m)
        x = [m(1,:); m(1,:); m(2,:)-1; m(2,:)-1; m(1,:)];
        y = repmat ([oo(2,r); oo(4,r); oo(4,r); oo(2,r); oo(2,r)] - 0.5, ...
                    [1, size(m,2)]);
        line (x, y, 'color', 'r');
      end
    end
  end
end

function draw_plain_rectangles (rr, samPerDay)
  if isempty (rr)
    return
  end
  oo = [rr(:).on_off];
  oo([1,3],:) = oo([1,3],:) - 0.5;
  % Ensure consistency of before/after midnight settings
  oo([2,4],:) = mod1 (oo([2,4],:), samPerDay);

  if isfield(rr, 'colour')
    colour = { rr.colour }';
    a = cellfun (@isempty, colour);
    colour(a,:) = repmat ({[1 1 1]'}, 1, sum (a));
  else
    colour = repmat ({[1 1 1]'}, 1, length (rr));
  end
  if isfield (rr, 'drift')
    oo_raw_drift = bsxfun (@plus, oo([2,4],:), [rr.drift] .* (oo(3,:) - oo(1,:)));
    oo_drift = oo_raw_drift;
    oo_drift(oo_drift < 1) = oo_drift(oo_drift < 1) + samPerDay;
    idx = (oo_drift > samPerDay + 1);
    oo_drift(idx) = oo_drift(idx) - samPerDay;
  else
    oo_raw_drift = oo([2,4],:);
    oo_drift = oo([2,4],:);
  end

  %s = (min (oo(4,:), oo_drift(2,:)) > max (oo(2,:), oo_drift(1,:)));   % "simple" case
  s = all (oo_drift == oo_raw_drift) & (oo(4,:) > oo(2,:));
  if any (s)
    set (gca, 'colorOrder', [colour{s}]');
    x = [oo(1,s); oo(1,s); oo(3,s); oo(3,s); oo(1,s)];
    y = [oo(2,s); oo(4,s); oo_drift(2,s); oo_drift(1,s); oo(2,s)]-0.5;
    line (x, y);
  end

  s = ~s;   % "split over days" case
  if any (s)
    % next-easiest case: individual edges not split
    t = (max (oo(4,:), oo_drift(2,:)) < min (oo(2,:), oo_drift(1,:))) & s;
    if any (t)
      set (gca, 'colorOrder', [colour{t}]');
      one = ones (1, sum (t));
      x = [oo(1,t); oo(1,t); oo(3,t); oo(3,t)];
      hrs_1 = 1 + samPerDay;
      y1 = [hrs_1*one; oo(2,t); oo_drift(1,t); hrs_1*one];
      y2 = [one; oo(4,t); oo_drift(2,t); one];
      line ([x, x+1], [y1, y2] - 0.5);
    end

    if any (s & ~t)
      % Messy case: individual edge is split
      for t = find (s & ~t)
        if isfield (rr, 'drift')
          drift = rr(t).drift;
        else
          drift = 0;
        end
        % TODO: Fix this!!
        set (gca, 'colorOrder', [colour{t}]');
        hrs_1 = 1 + samPerDay;
        x1 = [oo(1,t); oo(1,t); oo(3,t); oo(3,t)];
        x2 = x1 + 1;
        y1 = [hrs_1; oo(2,t); oo_drift(1,t); hrs_1];
        y2 = [1; oo(4,t); oo_drift(2,t); 1];

        x_t = [];
        x_b = [];
        y_t = [];
        y_b = [];
        if drift > 0
          if oo_drift(1,t) - oo(2,t) < 0      % cross midnight
            x_t = oo(1,t) + [0; 0; 0; (hrs_1-oo(2,t))/drift];
            y_t = [hrs_1; hrs_1; oo(2,t); hrs_1];
            x1 = [x_t(end); x_t(end); x1(3:4)] + 1;
            y1 = [1; 1; oo_drift(1,t); oo_drift(2,t)];
            if oo_drift(2,t) - oo(4,t) >= 0   % only top crosses
              y2(4) = oo_drift(2,t);
            else
              y_t(1:2) = oo(4,t);
            end
          end
          if oo_drift(2,t) - oo(4,t) < 0      % cross midnight
            x_b = [oo(1,t)+(hrs_1-oo(4,t))/drift; oo([3,3,3],t)] + 1;
            y_b = [1; oo_drift(2,t); 1; 1];
            x2 = [x2(1:2); x_b(1); x_b(1)] - 1;
            y2 = [oo(2,t); oo(4,t); hrs_1; hrs_1];
            if oo_drift(1,t) - oo(2,t) >= 0   % only bottom crosses
              y1(1) = oo(4,t);
            else
              y_b(3:4) = oo_drift(1,t);
            end
          end
        else
          if oo_drift(1,t) - oo(2,t) > 0      % cross midnight
            x_t = oo(3,t) + [0; 0; 0; (hrs_1-oo_drift(1,t))/drift];
            y_t = [hrs_1; hrs_1; oo_drift(1,t); hrs_1];
            x1 = [x1(1:2); x_t(end); x_t(end)];
            y1 = [oo(4,t); oo(2,t); 1; 1];
            if oo_drift(2,t) - oo(4,t) <= 0   % only top crosses
              y2(1) = oo(2,t);
            else
              y_t(1:2) = oo_drift(2,t);
            end
          end
          if oo_drift(2,t) - oo(4,t) > 0      % cross midnight
            x_b = [oo([1,1],t); oo(3,t)+(hrs_1-oo_drift(2,t))/drift; oo(1,t)];
            y_b = [1; oo(4,t); 1; 1];
            x2 = [x_b(3); x_b(3); x2(3:4)] - 1;
            y2 = [hrs_1; hrs_1; oo_drift(2,t); oo_drift(1,t)];
            if oo_drift(1,t) - oo(2,t) <= 0   % only bottom crosses
              y1(4) = oo_drift(2,t);
            else
              y_b(1:2) = oo(4,t);
            end
          end
        end

        line ([x1, x2, x_t, x_b], [y1, y2, y_t, y_b] - 0.5);
      end
    end
  end
end

function draw_rectangle_DoW (r, samPerDay)
  if ~isfield (r, 'drift') || r.drift == 0
    draw_rectangle_DoW_no_drift (r, samPerDay)
    return;
  end

  oo = [r.on_off; r.alt_hrs];
  oo_raw_drift = oo([2, 4, 5, 6]) + r.drift .* (oo(3) - oo(1));
  oo_drift = oo_raw_drift;
  oo_drift(oo_drift < 1) = oo_drift(oo_drift < 1) + samPerDay;
  idx = (oo_drift > samPerDay + 1);
  oo_drift(idx) = oo_drift(idx) - samPerDay;

  % If we cross midnight, it's too hard.  Just draw a horizontal rectangle.
  % Perhaps we could divide it into weeks or something.
  % s = all (oo_drift == oo_raw_drift) && (oo(4) > oo(2)) && (oo(6) > oo(5));
  % if ~s
  %   draw_rectangle_DoW_no_drift (r, samPerDay)
  %   return;
  % end

  week = r;
  last = floor ((oo(3) - oo(1)) / 7) - 1;
  ends = [r.on_off(1)+7*(0:last), r.on_off(3)];
  for i = 0:last
    week.on_off = r.on_off;
    week.on_off([2,4]) = week.on_off([2,4]) + 7*i*week.drift;
    week.on_off([1,3]) = ends([i+1, i+2]);
    week.alt_hrs = r.alt_hrs + 7*i*week.drift;
    omit_edges = [i~=0, i~=last];
    draw_rectangle_DoW_no_drift (week, samPerDay, omit_edges);
  end
end


function draw_rectangle_DoW_no_drift (r, samPerDay, omit_edges)
  if nargin < 3
    omit_edges = [0, 0];
  end

  if isfield(r, 'colour') && ~isempty (r.colour)
    colour = r.colour;
  else
    colour = 'w';
  end

  if ~isfield (r, 'alt_hrs') || length(r.alt_hrs) < 2
    line (r.on_off([1,3]-0.5), r.on_off([2,2])-0.5, 'Color', colour);
    line (r.on_off([1,3]-0.5), r.on_off([4,4])-0.5, 'Color', colour);

    % Find on and off times to use for the vertical edges.
    % Same both sides, since neither is an alt time
    on_l  = r.on_off(2);
    off_l = r.on_off(4);
    on_r  = on_l;
    off_r = off_l;
  else

    range = r.on_off(1):r.on_off(3);
    alt_days_idx = (ismember (1 + mod (range - 1, 7), r.alt_days));

    if ~any (alt_days_idx)
      draw_plain_rectangles (r, samPerDay);
      return;
    elseif all (alt_days_idx)
      rr.on_off = r.on_off;
      rr.on_off([2,4]) = r.alt_hrs;
      draw_plain_rectangles (rr, samPerDay);
      return;
    end

    alt_days = ranges (range (alt_days_idx));
    if ~isempty (alt_days)
      alt_days(1,:) = alt_days(1,:) - 0.5;
      if length (alt_days) >= 2
        alt_days(2,:) = alt_days(2,:) + 0.5;
      else
        alt_days(2,:) = alt_days(1,:) + 0.5;
      end
    end

    normal_days = ~alt_days_idx;
    normal_days = ranges (range (normal_days(1:end-1)));
    if ~isempty (normal_days)
      normal_days(1,:) = normal_days(1,:) - 0.5;
      normal_days(2,:) = normal_days(2,:) + 0.5;
    end
    if r.on_off(2) < r.on_off(4)
      normal_ends = normal_days;
    else
      normal_ends = min (normal_days + 1, r.on_off(3) - 0.5);
    end

    % horizontal lines
    if ~isempty (normal_days)
      line (normal_days, r.on_off([2,2])-0.5, 'Color', 'w');
      line (normal_ends, r.on_off([4,4])-0.5, 'Color', 'w');
    end
    if ~isempty (alt_days)
      if r.alt_hrs(1) < r.alt_hrs(2)
        alt_ends = alt_days;
      else
        alt_ends = alt_days + 1;
      end
      line (alt_days, r.alt_hrs([1,1])'-0.5, 'Color', 'm');
      line (alt_ends, r.alt_hrs([2,2])'-0.5, 'Color', 'm');
    end

    % Find on and off times to use for the end vertical edges
    if ~isempty (normal_days) && normal_days(1) < alt_days(1)
      s = 1;
      on_l  = r.on_off(2);
      off_l = r.on_off(4);
    else
      s = 0;
      on_l  = r.alt_hrs(1);
      off_l = r.alt_hrs(2);
    end
    if ~isempty (normal_days) && normal_days(end) > alt_days(end)
      e = 1;
      on_r  = r.on_off(2);
      off_r = r.on_off(4);
    else
      e = 0;
      on_r  = r.alt_hrs(1);
      off_r = r.alt_hrs(2);
    end

    % vertical lines between days
    if ~isempty (normal_days)
      normal_days = normal_days(1+s:end-e);
      X = [normal_days(:)'; normal_days(:)'];
      X_next = min (X + 1, r.on_off(3) - 0.5);
      all_hrs = [r.alt_hrs(1), r.alt_hrs(2); r.on_off(2), r.on_off(4)];
      if max (all_hrs(:)) - min (all_hrs(:)) < samPerDay / 2    % all on one day
        line (X, [r.on_off(2); r.alt_hrs(1)] - 0.5, 'Color', 'w');
        line (X, [r.on_off(4); r.alt_hrs(2)] - 0.5, 'Color', 'w');
      else
        d = max (all_hrs) - min (all_hrs);    % columnwise
        if d(1) < samPerDay / 2
          line (X, [r.on_off(2); r.alt_hrs(1)] - 0.5, 'Color', 'w');
        else
          line (X, [max(all_hrs(:,1))-0.5, samPerDay+0.5], 'Color', 'w');
          line (X_next, [0.5, min(all_hrs(:,1))-0.5], 'Color', 'w');
        end
        if d(2) < samPerDay / 2
          if r.on_off(4) < r.on_off(2)
            XX = X_next;
          else
            XX = X;
          end
          line (XX, [r.on_off(4); r.alt_hrs(2)] - 0.5, 'Color', 'w');
        else
          line (X, [max(all_hrs(:,2))-0.5, samPerDay+0.5], 'Color', 'w');
          line (X_next, [0.5, min(all_hrs(:,2))-0.5], 'Color', 'w');
        end
      end
    end
  end

  % vertical lines at the ends
  if ~omit_edges(1)
    if on_l < off_l
      line (r.on_off([1,1])-0.5, [on_l, off_l]-0.5, 'Color', colour);
    elseif on_l ~= off_l
      line (r.on_off([1,1])-0.5, [on_l-0.5, 48.5],  'Color', colour);
      line (r.on_off([1,1])+0.5, [0.5,  off_l-0.5], 'Color', colour);
    end
  end

  if ~omit_edges(2)
    if on_r < off_r
      line (r.on_off([3,3])-0.5, [on_r, off_r]-0.5, 'Color', colour);
    elseif on_r ~= off_r
      line (r.on_off([3,3])-0.5, [on_r-0.5, 48.5],  'Color', colour);
      line (r.on_off([3,3])+0.5, [0.5,  off_r-0.5], 'Color', colour);
    end
  end
end
