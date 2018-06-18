function [aircon, ok, confidence] = aircon_ground_truth (cv, aircon, house_no, total_cor, confidence, meta)
  dims = get (0, 'ScreenSize');
  if isempty (timer)
    dims = dims/2;
  end

  % TODO:
  % Align dates of aircon not overlapping in time
  % Allow timer to be set to produce "non-aircon"
  % "Snap" to correct start / end date of aircon
  % Change alt_day of overlapping aircon
  % Make changing alt_day easier
  % Insert gap
  %
  % BUG FIXES:
  %



  f = figure ('Visible', 'off', 'Position', [0, 0.5*dims(4), 0.8*dims(3), 0.8*dims(4)]);

  ha = axes ('Units', 'pixels', 'Units', 'normalized', ...
             'Position', [0.05, 0.05, 0.7, 0.9]);

  imagesc (cv);
  tic;

  ok = true;
  if all (confidence == 0)
    confidence = ones (size (cv));
  end
  undo_stack = {};
  selected_ac = [];
  ac_colour = [0, 0, 0];
  background = cv;
  bg = 'cv';
  weekends = [];
  temperatures = meta.temperatures';

  uicontrol ('Style', 'Text', 'String', num2str (house_no), ...
             'Units', 'normalized', ...
             'Position', [0.8, 0.96, 0.15, 0.04]);

  %showHide =
  uicontrol ('Style', 'pushbutton', 'String', 'Bkgrnd', ...
             'Units', 'normalized', ...
             'Position', [0.8, 0.92, 0.075, 0.04], ...
             'Callback', {@choose_background});

  uicontrol ('Style', 'pushbutton', 'String', 'Show a/c', ...
             'Units', 'normalized', ...
             'Position', [0.875, 0.92, 0.075, 0.04], ...
             'Callback', {@cycle_ac_colour});

  uicontrol ('Style', 'pushbutton', 'String', 'Weekends',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.84, 0.15, 0.04], ...
             'Callback', {@toggle_weekends});

  uicontrol ('Style', 'pushbutton', 'String', '<--',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.75, 0.05, 0.04], ...
             'Callback', {@show_left});

  uicontrol ('Style', 'pushbutton', 'String', 'O',...
             'Units', 'normalized', ...
             'Position', [0.85, 0.75, 0.05, 0.04], ...
             'Callback', {@show_all});

  uicontrol ('Style', 'pushbutton', 'String', '-->',...
             'Units', 'normalized', ...
             'Position', [0.9, 0.75, 0.05, 0.04], ...
             'Callback', {@show_right});

  %del_button =
  uicontrol ('Style', 'pushbutton', 'String', 'Delete aircon',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.68, 0.15, 0.04], ...
             'Callback', {@del}, ...
             'ButtonDownFcn', {@del});

  %add_button =
  uicontrol ('Style', 'pushbutton', 'String', 'Add aircon',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.60, 0.15, 0.04], ...
             'Callback', {@add});

%{
  uicontrol ('Style', 'pushbutton', 'String', 'Merge',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.68, 0.15, 0.04], ...
             'Callback', {@merge_selected});

  uicontrol ('Style', 'pushbutton', 'String', 'Split aircon',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.60, 0.15, 0.04], ...
             'Callback', {@set_splitting});

  uicontrol ('Style', 'pushbutton', 'String', 'Alt Days',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.52, 0.15, 0.04], ...
             'Callback', {@toggle_alt_days});
%}

  uicontrol ('Style', 'Text', 'String', 'Confidence', ...
             'Units', 'normalized', ...
             'Position', [0.8, 0.44, 0.15, 0.04]);

  handle_conf = ...
  uicontrol ('Style', 'slider', 'String', 'Confidence',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.40, 0.15, 0.04], ...
             'Min', 0, 'Max', 10, 'SliderStep', [1, 1], 'Value', 10);
  if isfield (aircon, 'manual_confidence')
    set (handle_conf, 'Value', max ([aircon.manual_confidence]));
  end

  uicontrol ('Style', 'pushbutton', 'String', 'Undo',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.30, 0.15, 0.04], ...
             'Callback', {@undo});

  %done_button =
  uicontrol ('Style', 'pushbutton', 'String', 'Next',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.16, 0.15, 0.06], ...
             'Callback', {@done});

  uicontrol ('Style', 'pushbutton', 'String', 'Prev',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.08, 0.15, 0.04], ...
             'Callback', {@prev});

  uicontrol ('Style', 'pushbutton', 'String', 'Skip',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.02, 0.075, 0.04], ...
             'Callback', {@cancel}, ...
             'ButtonDownFcn', {@cancel});

  uicontrol ('Style', 'pushbutton', 'String', 'Quit',...
             'Units', 'normalized', ...
             'Position', [0.875, 0.02, 0.075, 0.04], ...
             'Callback', {@exit_loop});

  set (f, 'WindowButtonDownFcn', @start_select_or_move);
  set (f, 'WindowButtonUpFcn',   @final_select_or_move);

  set (f, 'KeyPressFcn', @key_press);

  % initialize for use below
  selected_ac = [];
  downtime = 0;
  coordinates = 0;

  if (~exist ('OCTAVE_VERSION', 'builtin'))
    movegui (f, 'center');
  end
  
  sanity_check;
  
  draw_aircon;     % will set aircon_visible = y
  set (f, 'Visible', 'on');
  
  uiwait (f);

  function deleted = sanity_check
    % Add single-slot cooling
    % Find mean of "cooling", mean of summer except cooling
    mean_on  = mean (cv(aircon ~= 0));
    mean_off = mean (cv(aircon == 0));
    sd_off = sqrt (var (cv(aircon == 0)));
    missed = cv;
    missed = missed - mean_off;
    missed(missed < (3 * mean_on - mean_off) / 4 + sd_off) = 0;
    missed(temperatures < 25) = 0;
    % When in place:
    %     smoothed(missed) = me_year(missed)
    %     do maximum_likelihood
    missed(aircon ~= 0) = 0;
    figure(14); imagesc (missed);
    figure(12); draw_aircon;
    aircon = aircon + missed;

    % Remove cooling from cold times.
    max_temp = max (temperatures);
    prev_max = [max_temp(1), max_temp(1:end-1)];
    cold_days = (max_temp < 15 ...
              | (max_temp < 20 & prev_max < 25));
    cool_days = (max_temp < 25 & max_temp + prev_max < 45);
    cooling_on_cold = aircon(:, (cold_days | cool_days));
    deleted = {};     % message to be displayed
    if any (cooling_on_cold(:))
      cold_hrs = temperatures < 18;
      cold_hrs = cold_hrs & bsxfun (@and, cool_days, true (size (cv, 1), 1));
      cold_hrs = cold_hrs | bsxfun (@and, cold_days, true (size (cv, 1), 1));
      range = find (cold_hrs & aircon);
      if ~isempty (range)
        add_undo ({'replace', range, aircon(range)});
        aircon(range) = 0;

        [hrs, days] = ind2sub (size (cv), ranges(range));
        msg = sprintf ('(%d %d) - (%d %d)\n', [days(:)'; hrs(:)']);
        msg = strsplit (msg, '\n');
        deleted = {'sanity check deleted', msg{:}};
      else
        fprintf ('cold day cooling\n');
        % keyboard;
      end
    end

    if ~any (confidence)    % Don't overwrite manual tweaks
      % Tweak ends
      % Trim by one or extend by several.
      runs = ranges(find (aircon));
      if ~isempty (runs)
        hrs = size (cv, 1);
        offsets = [-hrs; -hrs-1; -1; hrs-1; hrs];

        % Extend runs, on slot at a time
        changed = true;
        old_aircon = aircon;
        figure(12); draw_aircon;
        while changed
          changed = false;
          incr = -1;
          for i = 1:2
            next = runs(i,:) + incr;
            next(end) = min (next(end), numel (cv));
            next(1) = max (next(1), 1);
            ref_val = cv(runs(i,:)) + 0.05;
            prev_day = runs(i,:);
            next_day = prev_day;

            % ignore if smaller than previous day non-aircon neighbour
            idx = next > hrs;
            prev_day(idx) = next(idx) - hrs;
            prev_ok = (cv(prev_day) < cv(next) | aircon(prev_day) > 0);
            % If aircon neighbour, use it when checking drop in power is small
            idx = aircon(prev_day) > 0;
            ref_val(idx) = max (ref_val(idx), cv(prev_day(idx)));

            % ignore if smaller than next day non-aircon neighbour
            idx = next <= numel (cv) - hrs;
            next_day(idx) = next(idx) + hrs;
            next_ok = (cv(next_day) < cv(next) | aircon(next_day) > 0);
            % If aircon neighbour, use it when checking drop in power is small
            idx = aircon(next_day) > 0;
            ref_val(idx) = max (ref_val(idx), cv(next_day(idx)));

            ratio = cv(next) ./ ref_val;

            idx = (ratio >= 0.5 & prev_ok & next_ok) ...
                 |(ratio >= 0.75 & (prev_ok | next_ok)) ...
                 |(ratio >= 1);
            [a, b] = ind2sub (size (cv), runs(i,idx));
            [a(:)'; b(:)'; ratio(idx); cv(next(idx)); ref_val(idx)]
            if any (idx)
              next = next(idx);
              prev = runs(i,idx);
              aircon(next) = aircon(prev) * (cv(next) / cv(prev));
              runs(i,idx) = next;
              changed = true;
            end

            incr = 1;
          end
        end
        figure(13); draw_aircon;

        % Try to remove excess from ends
        for i = 1:2
          neighbours = bsxfun(@plus, runs(i,:), offsets);
          valid_neighbours = neighbours;
          % Point out-of-bounds neighbours to a point with aircon ~= 0
          valid_neighbours(neighbours < 1 | neighbours > numel (aircon)) = runs(1);
          % Set all invalid neighbours (including out-of-bounds) to NaN
          valid_neighbours(aircon(valid_neighbours) ~= 0) = NaN;
          idx = isfinite (valid_neighbours);
          % Replace index by power value
          valid_neighbours(idx) = cv(valid_neighbours(idx));

          % Find mean and variance of neighbours
          means = mean (valid_neighbours, 1, 'omitnan');
          deviations = bsxfun (@minus, valid_neighbours, means);
          deviations(isnan (deviations)) = 0;
          vars = diag (deviations' * deviations)';
          vars = vars ./ sum (isfinite (valid_neighbours));

          % ends that are indistinguishable from noise
          small_ends = (cv(runs(i,:)) < 1.1 * means + sqrt (vars));
          two_step = (runs(2,:) > runs(1,:));
          small_ends(two_step) = small_ends(two_step) ...
                             &   aircon(runs(i,two_step) - offsets(3)) ...
                               > aircon(runs(i,two_step)) * 2;

          % Truncate ends indistinguishable from noise
          if any (small_ends)
            r = runs(i, small_ends);
            add_undo ({'replace', r, aircon(r)});
            figure(10); draw_aircon;
            aircon(r) = 0;
            figure(11); draw_aircon;
          end

          % Second iteration uses neighbours later than itself
          offsets = -offsets;
        end
      end
    end
  end

  function choose_background (~, ~)
    if isequal (bg, 'cv')
      background = total_cor;
      bg = 'cor';
    elseif isequal (bg, 'cor')
      background = temperatures;
      bg = 'temperatures';
    elseif isequal (bg, 'temperatures')
      background = confidence;
      bg = 'conf';
%       background = ones (size (background));
%       background (:, evalin('base', 'meta.weekends')) = 1.5;
%       background(1:2) = [0, 3];    % Force imagesc's scaling
%       bg = 'DoW';
    else
      background = cv;
      bg = 'cv';
    end
    
    draw_aircon;
  end

  function toggle_weekends (~, ~)
    if isempty (weekends)
      weekends = meta.weekends;
    else
      weekends = [];
    end
    draw_aircon;
  end

  function cycle_ac_colour (~, ~)
    if isempty (ac_colour)
      ac_colour = [0, 0, 0];
    elseif ac_colour(1) == 0
      ac_colour = [1, 1, 1];
    else
      ac_colour = [];
    end
    draw_aircon;
  end

  function show_left (~, ~)
    set (gca, 'XLim', [0.5, 120]);
  end

  function show_right (~, ~)
    set (gca, 'XLim', [265, size(cv,2)+0.5]);
  end

  function show_all (~, ~)
    set (gca, 'XLim', [0.5, size(cv,2)+0.5]);
    set (gca, 'YLim', [1, size(cv,1)]);
  end

  function draw_aircon (~, ~)
    show_aircon (aircon, background, ac_colour, weekends);
    draw_as_selected (selected_ac);
  end

%{
  function set_splitting (~, ~)
    set (f, 'WindowButtonUpFcn', @split_aircon);
  end

  function split_aircon (~,~)
    % split a aircon into two, leaving a gap of size determined by
    % the mouse movement.  The split is horizontal if most mouse movement
    % (as a fraction of the size of the aircon) is vertical and vice
    % versa.
    set (f, 'WindowButtonUpFcn', @final_select_or_move); % restore default

    old_coords = coordinates;
    coordinates_now = get (ha, 'CurrentPoint');
    coordinates_now = coordinates_now (1, 1:2);
    move = coordinates_now - old_coords;
    [r, found] = get_aircon_from_coords ();  % clobbers 'coordinates'
    if ~found
      return
    end
    aircon(end+1) = aircon (r);
    add_undo ({ 'replace', r, aircon(r) });
    add_undo ({ 'delete', length(aircon), [] });

    width = aircon(r).on_off(3) - aircon(r).on_off(1);
    height = aircon(r).on_off(4) - aircon(r).on_off(2);
    if move(1) / width >  move(2) / height
      % vertical split
      aircon(end).on_off(4) = old_coords (2);
      aircon(r).on_off(2) = coordinates_now (2);
    else
      % horizontal split
      [en, st] = align_with_existing (old_coords, coordinates_now);
      aircon(end).on_off(3) = round (en);
      aircon(r).on_off(1) = round (st);
    end
    draw_aircon;
  end
%}

  function del (object_handle, ~)
    if nargin > 0 && ~isempty (object_handle)   % If used as callback, not explicitly
      modifiers = get_modifiers (object_handle);
      if ismember('control',modifiers)
        selected_ac = ranges (find (aircon));
      end
    end

    if isempty (selected_ac)
      return;
    end
    e = length (undo_stack);
    undo_stack{e+length (selected_ac)} = [];   % pre-allocate to avoid Matlab warning
    for i = 1:size (selected_ac, 2)
      r = selected_ac(:,i);
      add_undo ({'replace', r(1):r(2), aircon(r(1):r(2))});
      aircon(r(1):r(2)) = 0;
    end
    selected_ac = [];

    draw_aircon;
  end

  function yes = inorder (st, mid, en)
    % Copied from has_pool_pump
    % True if time  mid  is in the interval  st:en, mod  x  for some x>max(en,st)
    mask = (st < en);
    if length (mid) == 1 && length (mask) > 1
      mid = repmat (mid, size (mask));
    end
    yes(mask)  = st(mask)  < mid(mask)  & mid(mask)  < en(mask);
    yes(~mask) = st(~mask) < mid(~mask) | mid(~mask) < en(~mask);
  end

  function select_aircon (object_handle, ~)
    size (aircon)
    [r, found] = get_aircon_from_coords;

    if found
      modifiers = get_modifiers (object_handle);
      if ismember('shift',modifiers) || ismember('control', modifiers)
        selected_ac(:,end+1) = r;
      else
        draw_aircon;
        selected_ac = r;
      end
      draw_as_selected (r);
    else
      selected_ac = [];
      draw_aircon;
    end
  end

  function draw_as_selected (r)
    [y, x] = ind2sub (size (cv), r);
    for i = 1:size (x, 2)
      if x(1,i) == x(2, i)
        line(x, bsxfun (@plus, y, [-0.5; 0.5]), 'Color', [1, 0, 0]);
      else
        line ([x(1, i), x(1, i)], [y(1, i)-0.5, size(cv, 1)+0.5], 'Color', [1, 0, 0]);
        line ([x(2, i), x(2, i)], [0.5, y(2, i)+0.5], 'Color', [1, 0, 0]);
        for j = x(1, i)+1 : x(2, i)-1
          line ([j, j], [0.5, size(cv, 1)+0.5], 'Color', [1, 0, 0]);
        end
      end
    end
  end

  function [r, found] = get_aircon_from_coords
    % Find point identified as an airconditioner near mouse cursor.
    % found = true if exact selected point has aircon, -1 if nearby point.
    coordinates = get (ha, 'CurrentPoint');
    coordinates = coordinates (1, 1:2);
    x = coordinates (1, 1);
    y = coordinates (1, 2);
    found = false;
    r = [];

    if isfield (aircon, 'uncertain')    % if cancelling
      return;
    end
    xx = round (x);
    yy = round (y);
    if xx > size (cv,2)+0.5 || yy >= size (cv,1)+0.5 || any ([xx yy] < 0.5)
      return;
    end

    if aircon (yy, xx) > 0
      found = true;
    else
      for radius = 1:10
        xxx = max (1, xx - radius) : min (size (cv,2), xx + radius);
        yyy = max (1, yy - floor(radius/5)) : min (size (cv,1), yy + radius/5);
        [xxx, yyy] = find (aircon (yyy, xxx), 1);
        if ~isempty (xxx)
          xx = xx + xxx-1;
          yy = yy + yyy-1;
          found = -radius;
          break;
        end
      end
    end
    
    if found
      r = (xx - 1) * size (cv, 1) + yy;
      r1 = sum (r - find (aircon (r:-1:1) == 0, 1)) + 2;  % sum([]) = 0
      r2 = r + find (aircon (r:end)  == 0, 1) - 2;
      if isempty (r2)
        r2 = length (aircon(:));
      end
      if r1 <= r2
        r = [r1; r2];
      else
        r = [];
        found = false;
      end
    end
  end

  function start_select_or_move (~, ~)
    coordinates = get (ha, 'CurrentPoint');
    coordinates = coordinates (1, 1:2);
    downtime = toc;
  end

  function modifiers = get_modifiers (h)
    s = '';
    % Find actual object that contains selection type
    while isempty (s)
      try
        s = get (h, 'selectiontype');
        fprintf ('got selectiontype\n');
      catch
        try
          h = get (h, 'Parent');
          fprintf ('got parent\n');
        catch
          s = 'oops';
        end
      end
    end

    if (strcmp (s, 'extend'))
      modifiers = {'shift'};
      display (modifiers);
    elseif (strcmp (s, 'alt'))
      modifiers = {'control'};
      display (modifiers);
    else
      modifiers = {};
    end
  end

  function final_select_or_move (object_handle, event_data)
    modifiers = get_modifiers (object_handle);

    coordinates_now = get (ha, 'CurrentPoint');
    coordinates_now = coordinates_now (1, 1:2);
    move = coordinates_now - coordinates;
    move(1) = round (move(1));
    scaled = [move(1), 5*move(2)];

    a = ranges (find (aircon(:)));
    if isempty (a)
      on_off = [];
    else
      [on_off1a, on_off1b] = ind2sub (size (cv), a(1,:));
      [on_off2a, on_off2b] = ind2sub (size (cv), a(2,:));
      on_off = [on_off1b(:)'; on_off1a(:)'; on_off2b(:)'; on_off2a(:)'];
    end
    mn = min (coordinates_now, coordinates);
    mx = max (coordinates_now, coordinates);

    % aircon entirely within dragged region
    r = [];
    if ~isempty (on_off)
      r = on_off(1,:) > mn(1) & on_off(2,:) > mn (2) ...
          & on_off(3,:) < mx(1) & on_off(4,:) < mx(2) ...
          & on_off(4,:) >= on_off(2,:);
    elseif ismember ('normal', modifiers)
      return
    end

    if toc - downtime < 0.5 && norm (scaled) < 3
      % assume a click, not a drag
      select_aircon (object_handle, event_data);
      if ismember ('control', modifiers)  %if right click
        auto_fix;
        draw_aircon;
      end
    elseif any (r)
      % select all aircon in the range
      if ismember ('shift', modifiers)
        selected_ac = unique ([selected_ac, a(:,r)]', 'rows')';
      else
        selected_ac = a(:, r);
      end
      draw_as_selected (selected_ac);
    elseif ismember ('shift', modifiers)
      coordinates_cache = coordinates;
      select_aircon (object_handle, event_data);
      coordinates = coordinates_cache;

   if isempty (selected_ac)
        add_aircon
      else
        x = round ((mn(1) + mx(1))/2);
        y = round (mn(2) + 0.5):round(mx(2) - 0.5);
        selected_ac = (x-1)*size(cv, 1) + round ([mn(2)+0.5; mx(2)-0.5]);
        draw_aircon;
     end
    else
      if ismember ('normal', modifiers)
        return
      end
      % Assume a drag.
      % Use both direction of drag and proximity to edges to work out
      % which edge was being dragged.

      % vertical move => horizontal edge
      seems_horiz  = (abs (scaled(1)) < abs (scaled(2)));
      no_way_horiz = (abs (scaled(1)) > 3 * abs (scaled(2)));

      if (isempty (on_off))
        if ismember ('control', modifiers)
          auto_fix;
        else
          add_aircon (object_handle, event_data);
        end
        draw_aircon;
        return
      end
%{

      vert_candidates = find (inorder(on_off(2,:) - 0.5, coordinates(2), ...
                                      on_off(4,:) - 0.5));
      edges = on_off([1,3], vert_candidates);
      [~, vert_closest] = min (abs (edges(:) - coordinates(1)));
      vert_closest = vert_closest / 2;
      if vert_closest == round (vert_closest)
        vert_closest = vert_candidates(vert_closest);
        vert_edge = 3;
      else
        vert_closest = vert_candidates(vert_closest + 0.5);
        vert_edge = 1;
      end
%}

      % If a horizontal drag ('vertical edge'), select the dragged runs
      if ~seems_horiz && (no_way_horiz || isempty (selected_ac))
        x1 = round (min (coordinates(1), coordinates_now(1)));
        x2 = round (max (coordinates(1), coordinates_now(1)));
        y = round ((coordinates(2) + coordinates_now(2)) / 2);

        if x1 > size (cv, 2) || x2 < 1
          return
        else
          if x2 > size (cv, 2)
            x2 = size (cv, 2);
          end
          if x1 < 1
            x1 = 1;
          end
        end

        x = (aircon(y, x1:x2) > 0);
        x = x1 + find (x) - 1;
        if ~isempty (x)
          pos = (x-1) * size (cv, 1) + y;
          a = ranges (find (aircon(:)));
          selected_ac = zeros (2, length (x));
          for i = 1:length (pos)
            col = a(1,:) <= pos(i) & a(2,:) >= pos(i);
            selected_ac(:,i) = a(:, col);
          end
          draw_as_selected (selected_ac);
        end
        return;
      end

      horiz_candidates = find (on_off(1,:) < coordinates(1) + 0.5 ...
                               & coordinates(1) < on_off(3, :) + 0.5);
      if isempty (horiz_candidates)
        return;
      end
      edges = [on_off(2, horiz_candidates) - 0.5, ...
               on_off(4, horiz_candidates) + 0.5];

      % If burst wraps midnight, it may have an on/off close to coord(2)
      % but be on the wrong day, so penalise mismatch.
      day_mismatch = [on_off(1,horiz_candidates)-round(coordinates(1)), ...
                      on_off(3,horiz_candidates)-round(coordinates(1))];

      [~, horiz_closest] = min (abs (edges - coordinates(2)) ...
                                + 10 * abs (day_mismatch));
      if horiz_closest > length (horiz_candidates)
        horiz_closest = horiz_candidates(horiz_closest ...
                                         - length (horiz_candidates));
        horiz_edge = 4;
      else
        horiz_closest = horiz_candidates(horiz_closest);
        horiz_edge = 2;
      end

      if seems_horiz && ~isempty (horiz_closest)
        closest = horiz_closest(1);
        edge = horiz_edge;
%{
      elseif ~isempty (vert_closest)
        closest = vert_closest(1);
        edge = vert_edge;
%}
      else
        return
      end

      if ismember ('control', modifiers)
      end

      % Move the edge by the amount of the mouse move
      if edge == 2
        if move(2) > 0
          interval = a(1, closest) : a(1, closest) + round (move(2)) - 1;
          add_undo ({'replace', interval, aircon(interval)});
          aircon(interval) = 0;
        else
          interval = a(1, closest) + round (move(2)) : a(1, closest) - 1;
          if ~isempty (interval)
            interval(aircon(interval) > 0) = [];  % don't change old vals
            add_undo ({'replace', interval, aircon(interval)});
            aircon(interval) = min (aircon(interval(end)+1), cv(interval));
          end
        end
      elseif edge == 4
        if move(2) > 0
          interval = a(2, closest) : a(2, closest) + round (move(2));
          interval(aircon(interval) > 0) = [];  % don't change old vals
          if ~isempty (interval)
            add_undo ({'replace', interval, aircon(interval)});
            aircon(interval) = min (aircon(interval(1)-1), cv(interval));
          end
        else
          interval = a(2, closest) + 1 + round (move(2)) : a(2, closest);
          add_undo ({'replace', interval, aircon(interval)});
          aircon(interval) = 0;
        end
      end

      draw_aircon;
    end
  end

  function add (~, ~)
    set (f, 'WindowButtonUpFcn', @add_aircon);
  end

%{
  function [st, en] = align_with_existing (down_pos, up_pos)
    % Find matching start and end times from existing aircon
    if isempty (aircon)
      match = [];
    else
      on_off = [aircon(:).on_off];
      match = find (abs (on_off(3,:) - down_pos(1)) < 2);
    end
    if ~isempty (match)
      st = on_off(3, match(1));
    else
      st = round (down_pos(1));
    end
    if ~isempty (aircon)
      match = find (abs (on_off(1,:) - up_pos(1)) < 2);
    end
    if ~isempty (match)
      en = on_off(1, match(1));
    else
      en = round (up_pos (1));
    end
  end
%}


  function add_aircon (~,~)
    coordinates_now = get (ha, 'CurrentPoint');
    coordinates_now = coordinates_now (1, 1:2);

    x1 = round (coordinates(1));
    x2 = x1;
    y1 = min (coordinates(2), coordinates_now(2)) + 0.5;
    y2 = max (coordinates(2), coordinates_now(2)) - 0.5;

    % 0.5 because aircon specifies midpoint of pixel.
    on_off = round ([x1; y1; x2; y2]);

    do_add_aircon (on_off);

    % Next time, go back to nudging aircon
    set (f, 'WindowButtonUpFcn', @final_select_or_move);
    draw_aircon;
  end

  function do_add_aircon (on_off)
    on_off([2,1,4,3]) = min (max (on_off([2,1,4,3]), 1), [size(aircon)'; size(aircon)']);
    r = sub2ind (size (aircon), on_off([2,4]), on_off([1,3]));
    range = r(1):r(2);
    range(aircon(range) > 0) = [];
    aircon(range) = min(cv(range));

    add_undo ({'delete', range, []});
  end

  function add_undo (action)
    undo_stack{end + 1} = action;
  end

  function undo (~, ~)
    if isempty (undo_stack)
      return
    end
    undo_action = undo_stack{end};
    undo_stack(end) = [];
    if isempty (undo_action)
      return;
    end
    r = undo_action{2};
    switch undo_action{1}
      case 'replace'
        aircon(r) = undo_action{3};
      case 'delete'
        aircon(r) = zeros (size (r));
      case 'conf'
        confidence(r) = undo_action{3};
        if isequal (bg, 'conf')
          background = confidence;
        end
%       case 'insert'
%         aircon(r+1:end+1) = aircon(r:end);
%         aircon(r) = undo_action{3};
    end
    draw_aircon;
  end

  function key_press (object_handle, event_data)
    switch event_data.Key
      case 'delete'
        if isempty (selected_ac)
          select_aircon (object_handle, event_data)
        end
        if ~isempty (selected_ac)
          del
        end
      case 'slash'
        if ~isempty(selected_ac)
          conf = get (handle_conf, 'Value');
          range = [];
          for i = 1:size (selected_ac, 2)
            % aircon(i).colour = [0.05; 0.05; 0.05] * (10+conf);
            % confidence(selected_ac(1,i):selected_ac(2,i)) = conf;
            range = [range, selected_ac(1,i):selected_ac(2,i)];
          end
          add_undo ({'conf', range, confidence(range)});
          confidence(range) = conf;
          selected_ac = [];
          if isequal (bg, 'conf')
            background = confidence;
          end
          draw_aircon;    % remove highlighting of selected_ac.
        end
      case 'm'
        merge_selected (object_handle, event_data);
      case 'u'
        undo (object_handle, event_data);
    end
  end

  function affected_acs = auto_fix
    % Guess what could be wrong.
    %  - close a gap
    %  - make timing match the other chain
    %  - missed-days
    %  - copy from other chain
    if coordinates(1,1) > size (cv, 2) || coordinates(1,2) > size (cv, 1)
      return
    end

    x = max (min (round (coordinates(1, 1)), size (cv, 2)), 1);
    y = max (min (round (coordinates(1, 2)), size (cv, 1)), 1);
    r = (x - 1) * size (cv, 1) + y;

%{
    if ~isempty (selected_ac) ...
        && selected_ac(1) <= r && selected_ac(end) >= r
      r = selected_ac(1);
      rr = aircon(r);
      selected_ac = [];
      
      return;
    end
%}

    if all (aircon == 0)
      [h, bins] = hist (cv(:), 20);
      peak = sort (secondPeak (h(:)', 2));
      bp = bins(peak);
      power = bp(2) - bp(1);
    else
      power = mean (aircon(aircon ~= 0));
    end

    selected_ac = [];
    if cv(y, x) > 0.5 * power
      jj = find_jumps (cv(:));
      % step = cv(2:end) ./ cv(1:end-1);
      % step = log (max (step, 1./step));
      % enough(r-1:-1:1) = cummin (cv(r-1:-1:1)) / cv(r);
      % enough(r+1:length (cv(:)))  = cummax (cv(r+1:end))  / cv(r);
      r1 = r - find (cv(r:-1:1) < 0.5 * cv(r), 1) + 2;
      r2 = r + find (cv(r:end)  < 0.5 * cv(r), 1) - 2;
      range = r1:r2;
      range(aircon(range) > 0) = [];  % only change values not already set
      add_undo ({'replace', range, aircon(range)});
      aircon(range) = min (cv(range), power);
      return
    end
  end

  function exit_loop (~, ~)
    ok = -2;
    delete (f);
    f = 0;
  end

  function prev (~, ~)
    ok = -1;
    delete (f);
    f = 0;
  end

  function done (~, ~)
    conf = get (handle_conf, 'Value');
    if conf ~= 10 && ~any (confidence(:) == conf) % If slider moved since last setting a confidence
      confidence(confidence == 10) = conf;
    elseif isempty (undo_stack)
      ok = false;       % skip saving if no changes
    end
    delete (f);
    f = 0;
  end

  function cancel (object_handle, ~)
    ok = false;
    modifiers = get_modifiers (object_handle);
    if ismember('control',modifiers)
      ok = 2;
    end
    done;
  end

end


% Functions copied from has_pool_pump.m
%{
function m = mod1 (value, modulus)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Modulo operator, returning values in [1, modulus], not [0, modulus-1].
  m = mod (value-1, modulus) + 1;
end
%}

function yes = between (A, mid, B, modulo_2)
  % Consider the shorter of the two "intervals" modulo 2*modulo_2
  % bounded by A and B.
  % Return true if mid in in that interval.
  if (abs (A - B) < modulo_2)     %usual case, not crossing midnight
    if A > B
      yes = (B < mid & mid < A);
    else
      yes = (A < mid & mid < B);
    end
  else
    if A < B
      yes = (B < mid | mid < A);
    else
      yes = (A < mid | mid < B);
    end
  end
end

%{
function yes = inorder (st, mid, en)
  % True if time  mid  is in the interval  st:en, mod  x  for some x>max(en,st)
  if st < en
    yes = st < mid && mid < en;
  else
    yes = st < mid || mid < en;
  end
end

function d = hr_diff (x, y, modulus)
 % hr_diff (x, y, modulus) is the difference, modulo modulus, between x and y
  d = abs (x - y);
  d = min (d, modulus - d);
end
%}
