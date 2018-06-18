function [timer, rectangles, ok, next] = pool_ground_truth (cv, timer, rectangles, house_no)
  dims = get (0, 'ScreenSize');
  if isempty (timer)
    dims = dims/2;
  end

  % TODO:
  % Align dates of rectangles not overlapping in time
  % Allow timer to be set to produce "non-rectangles"
  % "Snap" to correct start / end date of rectangle
  % Change alt_day of overlapping rectangles
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
  next = [];        % Don't specify next house, unless text field is edited
  undo_stack = {};
  redo_stack = {};
  days_to_hide = [];
  rectangles_visible = true;
  background = cv;
  bg = 'cv';

  sanitise;            % make sure rectangles are all valid
  draw_rectangles;     % will set rectangles_visible = y

  if nargin >= 4
    house_str = num2str (house_no);
  else
    house_str = '';
  end
  uicontrol ('Style', 'edit', 'String', house_str, ...
             'Units', 'normalized', ...
             'Position', [0.8, 0.96, 0.05, 0.04], ...
             'Callback', {@set_next_house});

  %showHide =
  uicontrol ('Style', 'pushbutton', 'String', 'Show rectangles', ...
             'Units', 'normalized', ...
             'Position', [0.8, 0.92, 0.07, 0.04], ...
             'Callback', {@draw_or_clear_rectangles});

  uicontrol ('Style', 'pushbutton', 'String', 'Background', ...
             'Units', 'normalized', ...
             'Position', [0.88, 0.92, 0.07, 0.04], ...
             'Callback', {@cycle_background});

  %del_button =
  uicontrol ('Style', 'pushbutton', 'String', 'Delete rectangle',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.84, 0.15, 0.04], ...
             'Callback', {@del}, ...
             'ButtonDownFcn', {@del});

  %add_button =
  uicontrol ('Style', 'pushbutton', 'String', 'Add rectangle',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.76, 0.15, 0.04], ...
             'Callback', {@add});

  uicontrol ('Style', 'pushbutton', 'String', 'Merge',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.68, 0.15, 0.04], ...
             'Callback', {@merge_selected}, ...
             'ButtonDownFcn', {@merge_selected});

  uicontrol ('Style', 'pushbutton', 'String', 'Split rectangle',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.60, 0.15, 0.04], ...
             'Callback', {@set_splitting});

  uicontrol ('Style', 'pushbutton', 'String', 'Alt Days',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.52, 0.07, 0.04], ...
             'Callback', {@toggle_alt_days});

  %add_button =
  uicontrol ('Style', 'pushbutton', 'String', 'Redo drift',...
             'Units', 'normalized', ...
             'Position', [0.88, 0.52, 0.07, 0.04], ...
             'Callback', {@redo_drift});

  uicontrol ('Style', 'Text', 'String', 'Confidence', ...
             'Units', 'normalized', ...
             'Position', [0.8, 0.44, 0.15, 0.04]);

  handle_conf = ...
  uicontrol ('Style', 'slider', 'String', 'Confidence',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.40, 0.15, 0.04], ...
             'Min', 0, 'Max', 10, 'SliderStep', [1, 1], 'Value', 10);
  if isfield (rectangles, 'manual_confidence')
    m = max ([rectangles.manual_confidence]);
    if ~isempty (m)
      set (handle_conf, 'Value', m);
      confidence_global = [];
    end
  end
  addlistener (handle_conf, 'Value', 'PostSet', @read_confidence);

  uicontrol ('Style', 'pushbutton', 'String', 'Undo',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.30, 0.10, 0.04], ...
             'Callback', {@undo});

  uicontrol ('Style', 'pushbutton', 'String', 'Redo',...
             'Units', 'normalized', ...
             'Position', [0.9, 0.30, 0.05, 0.04], ...
             'Callback', {@redo});

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
  selected_rect = [];
  confidence_global = [];
  conf_copied_from = [];
  downtime = 0;
  coordinates = 0;

  if (~exist ('OCTAVE_VERSION', 'builtin'))
    movegui (f, 'center');
  end
  set (f, 'Visible', 'on');

  uiwait (f);

  function draw_or_clear_rectangles (~, ~)
    if ~rectangles_visible
      draw_rectangles;
    else
      %plot_handle =
      if isequal (bg, 'diff')
        background = cv + rectangles_to_correction (rectangles, ...
                                        zeros (size (cv,2), size (cv,1)), ...
                                        1:size (cv, 2))';
      end
      tmp = background;
      for i = days_to_hide
        tmp(:,i:7:end) = 0;
      end
      imagesc (tmp);
      colorbar;
      rectangles_visible = false;
    end
  end

  function draw_rectangles (~, ~)
    if isequal (bg, 'diff')
      background = cv + rectangles_to_correction (rectangles, ...
                                      zeros (size (cv,2), size (cv,1)), ...
                                      1:size (cv, 2))';
    end
    show_rectangles (rectangles, background, days_to_hide);
    colorbar;
    rectangles_visible = true;
  end

  function cycle_background (~, ~)
    switch bg
      case 'cv'
        background = medfilt1 (cv, 5, [], 2);
        bg = 'med';
      case 'med'
        background = -rolling_min (-rolling_min (cv', 5), 5)';
        bg = 'roll_min';
      case 'roll_min'
        background = min (cv, 2);
        bg = 'min';
      case 'min'
        % background is set dynamically by draw_rectangles
        %background = cv + rectangles_to_correction (rectangles, zeros (size (cv,2), size (cv,1)), 1:size (cv, 2))';
        bg = 'diff';
      case 'diff'
        background = cv;
        bg = 'cv';
    end
    rectangles_visible = ~rectangles_visible;
    draw_or_clear_rectangles;
  end

  function redo_drift (~, ~)
    rectangles = find_drift (rectangles, cv, true);
    draw_rectangles;
  end

  function set_splitting (~, ~)
    set (f, 'WindowButtonUpFcn', @split_rectangle);
  end

  function split_rectangle (~,~)
    % split a rectangle into two, leaving a gap of size determined by
    % the mouse movement.  The split is horizontal if most mouse movement
    % (as a fraction of the size of the rectangle) is vertical and vice
    % versa.
    set (f, 'WindowButtonUpFcn', @final_select_or_move); % restore default

    old_coords = coordinates;
    coordinates_now = get (ha, 'CurrentPoint');
    coordinates_now = coordinates_now (1, 1:2);
    move = coordinates_now - old_coords;
    [r, found] = get_rectangle_from_coords ();  % clobbers 'coordinates'
    if ~found
      return
    end
    rectangles(end+1) = rectangles (r);
    add_undo ({ 'replace', r, rectangles(r) });
    add_undo ({ 'delete', length(rectangles), [] });

    width = rectangles(r).on_off(3) - rectangles(r).on_off(1);
    height = rectangles(r).on_off(4) - rectangles(r).on_off(2);
    if move(1) / width >  move(2) / height
      % vertical split
      rectangles(end).on_off(4) = old_coords (2);
      rectangles(r).on_off(2) = coordinates_now (2);
    else
      % horizontal split
      [en, st] = align_with_existing (old_coords, coordinates_now);
      rectangles(end).on_off(3) = round (en);
      if isfield (rectangles, 'drift')
        rectangles(r).on_off([2,4]) = rectangles(r).on_off([2,4]) ...
                                + rectangles(r).drift ...
                                  * (round (st) - rectangles(r).on_off(1));
      end
      rectangles(r).on_off(1) = round (st);
      rectangles(r).on_off = tweak_hrs (rectangles(r), cv, 2);
      rectangles(r).on_off = tweak_hrs (rectangles(r), cv, 4);
      rectangles(end).on_off = tweak_hrs (rectangles(end), cv, 2);
      rectangles(end).on_off = tweak_hrs (rectangles(end), cv, 4);
    end
    draw_rectangles;
  end

  function del (object_handle, ~)
    if nargin > 0 && ~isempty (object_handle)   % If used as callback, not explicitly
      modifiers = get_modifiers (object_handle);
      if ismember('control',modifiers)
        selected_rect = 1:length(rectangles);
      end
    end

    sr = -sort (-selected_rect);
    if isempty (sr)
      return;
    end
    keep_map = true (size (rectangles));
    for i = 1:length (sr)
      r = sr(i);
      add_undo ({'insert', r, rectangles(r)});
      keep_map(r) = false;
    end
    rectangles = rectangles(keep_map);
    selected_rect = [];
    confidence_global = [];

    draw_rectangles;
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

  function read_confidence (~, event)
    % Recommeded at [https://au.mathworks.com/matlabcentral/answers/264979-continuous-slider-callback-how-to-get-value-from-addlistener]
    % but gives error "Cannot find 'get' method for event.PropertyEvent class."
    %confidence_global = get (event, 'newValue');

    confidence_global = event.AffectedObject.Value;
    conf_copied_from = [];
  end

  function select_rectangle (object_handle, ~)
    [r, found] = get_rectangle_from_coords;

    if found
      modifiers = get_modifiers (object_handle);
      if ismember('shift',modifiers) || ismember('control', modifiers)
        selected_rect(end+1) = r;
      else
        draw_rectangles;
        selected_rect = r;

        % Display confidence level
        if isfield (rectangles(r), 'manual_confidence') ...
            && ~isempty (rectangles(r).manual_confidence) ...
            && isempty (confidence_global)
          set (handle_conf, 'Value', rectangles(r).manual_confidence);
          conf_copied_from = selected_rect;
        end
      end
      draw_as_selected (rectangles(r).on_off);
    elseif coordinates(1,1) <= size (cv, 2) + 1
      selected_rect = [];
      confidence_global = [];
      draw_rectangles;
    end
  end

  function draw_as_selected (oo)
    if isempty (oo)
      return
    end
    oo([1,3],:) = oo([1,3],:) - 0.5;

    s = (oo(4,:) > oo(2,:));   % "simple" case
    if any (s)
      x = [oo(1,s); oo(1,s); oo(3,s); oo(3,s)];
      y = [oo(2,s); oo(4,s); oo(4,s); oo(2,s)]-0.5;
      patch (x, y, 'w');
    end

    s = ~s;   % "split over days" case
    if any (s)
      one = ones (1, sum (s));
      x = [oo(1,s); oo(1,s); oo(3,s); oo(3,s)];
      y1 = [(size(cv,1)+1)*one; oo(2,s); oo(2,s); (size(cv,1)+1)*one];
      y2 = [one; oo(4,s); oo(4,s); one];
      patch ([x, x+1], [y1, y2] - 0.5, 'w');
    end

  end

  function [r, found] = get_rectangle_from_coords
    coordinates = get (ha, 'CurrentPoint');
    coordinates = coordinates (1, 1:2);
    x = coordinates (1, 1);
    y = coordinates (1, 2);
    found = false;

    if isfield (rectangles, 'uncertain')    % if cancelling
      r = 0;
      return;
    end

    % TODO: Vectorise
    for r = 1:length (rectangles)
      if rectangles(r).on_off(1)-0.5 <= x && rectangles(r).on_off(3)-0.5 >= x ...
          && inorder (rectangles(r).on_off(2)-0.5, y, rectangles(r).on_off(4)-0.5)
        found = true;
        break
      end
    end
  end

  function merge_selected (object_handle, ~)
    % Form one rectangle covering the extent covered by the highlighted
    % rectangles.

    % If right-click, set drift from the offsets between rectangles
    if length (selected_rect) >= 2
      if ismember ('control', get_modifiers (object_handle))
        % find drift
        oo = [rectangles.on_off];
        oo = oo(:, selected_rect);
        [~, idx] = sort (oo(1,:));
        oo = oo(:, idx);
        if isfield (rectangles, 'drift')
          drift = rectangles(selected_rect(1)).drift;
        else
          drift = 0;
        end
        changes = (oo(3,:) - oo(1,:)) * drift;
        ends = oo([2,4], :) + [changes; changes];
        jumps = oo([2,4], 2:end) - ends(:, 1:end-1);
        new_drift = (ends(2,end) - oo(4,1) + ends(1,end) - oo(2,1)) / (2 * (oo(3,end) - oo(1,1)));
        for i = 1:length (rectangles)
          add_undo ({ 'replace', i, rectangles(i) });
          rectangles(i).drift = new_drift;
          oo = rectangles(i).on_off;
          rectangles(i).on_off([2,4]) = mod1 (oo([2,4]) + 0.5 * (oo(3)-oo(1)) * (drift - new_drift), size (cv, 1));
        end
      end

      %actual merge
      s = sort (selected_rect);
      for i = length(s)-1:-1:1
        merge_rectangles (s(i), s(i+1));
        rectangles(s(1)).on_off = tweak_hrs (rectangles(s(1)), cv, 2);
        rectangles(s(1)).on_off = tweak_hrs (rectangles(s(1)), cv, 4);
      end
    end
    selected_rect = [];
    confidence_global = [];
    draw_rectangles;
  end

  function merge_rectangles (r1, r2)
    add_undo ({ 'replace', r1, rectangles(r1) });
    add_undo ({ 'insert', r2, rectangles(r2) });

    on_off1 = rectangles(r1).on_off;
    on_off2 = rectangles(r2).on_off;
    rectangles(r2) = [];

    rectangles(r1).on_off(1) = min (on_off1(1), on_off2(1));
    rectangles(r1).on_off(3) = max (on_off1(3), on_off2(3));
    if isfield (rectangles, 'drift')
      drift = rectangles(r1).drift * ([on_off1(1), on_off2(1)] - rectangles(r1).on_off(1));
    else
      drift = 0;
    end
    if on_off1(2) < on_off1(4)    % doesn't span midnight
      rectangles(r1).on_off(2) = min ([on_off1(2), on_off2(2)] - drift);
      rectangles(r1).on_off(4) = max ([on_off1(4), on_off2(4)] - drift);
    else
      rectangles(r1).on_off(2) = max ([on_off1(2), on_off2(2)] - drift);
      rectangles(r1).on_off(4) = min ([on_off1(4), on_off2(4)] - drift);
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
      catch
        try
          h = get (h, 'Parent');
        catch
          s = 'oops';
        end
      end
    end

    if (strcmp (s, 'extend'))
      modifiers = {'shift'};
    elseif (strcmp (s, 'alt'))
      modifiers = {'control'};
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

    on_off = [rectangles(:).on_off];
    mn = min (coordinates_now, coordinates);
    mx = max (coordinates_now, coordinates);

    % Rectangles entirely within dragged region
    r = [];
    if ~isempty (on_off)
      r = on_off(1,:) > mn(1) & on_off(2,:) > mn (2) ...
          & on_off(3,:) < mx(1) & on_off(4,:) < mx(2) ...
          & on_off(4,:) > on_off(2,:);
    elseif ismember ('normal', modifiers)
      return
    end

    if toc - downtime < 0.5 && norm (scaled) < 3
      % assume a click, not a drag
      select_rectangle (object_handle, event_data);
      if ismember ('control', modifiers)  %if right click
        auto_fix;
        draw_rectangles;
      end
    elseif any (r)
      % select all rectangles in the range
      if ismember ('shift', modifiers)
        selected_rect = unique ([selected_rect, find(r)]);
      else
        selected_rect = find (r);
      end
      draw_as_selected ([rectangles(selected_rect).on_off]);
    elseif ismember ('shift', modifiers)
      add_rectangle
    else
      if ismember ('normal', modifiers)
        return
      end
      % Assume a drag.
      % Use both direction of drag and proximity to edges to work out
      % which edge was being dragged.

      % vertical move => horizontal edge
      seems_horiz = (abs (scaled(1)) < abs (scaled(2)));

      if (isempty (on_off))
        if ismember ('control', modifiers)
          auto_fix;
        else
          add_rectangle (object_handle, event_data);
        end
        draw_rectangles;
        return
      end

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

      horiz_candidates = find (on_off(1,:) < coordinates(1) ...
                               & coordinates(1) < on_off(3, :));
      edges = [on_off(2, horiz_candidates), ...
               on_off(4, horiz_candidates)];
      if isfield (rectangles, 'alt_days')
        alt_cand = ~cellfun (@isempty, {rectangles(horiz_candidates).alt_days});
        alt_candidates = horiz_candidates(alt_cand);
        hrs = [rectangles(alt_candidates).alt_hrs]';
        edges = [edges, hrs(:)'];
      end

      [~, horiz_closest] = min (abs (edges - coordinates(2)));
      if horiz_closest > length (horiz_candidates)
        if horiz_closest <= 2*length (horiz_candidates)
          horiz_closest = horiz_candidates(horiz_closest ...
                                           - length (horiz_candidates));
          horiz_edge = 4;
        else
          horiz_closest = horiz_closest - 2*length (horiz_candidates);
          if horiz_closest > length (alt_candidates)
            horiz_closest = alt_candidates (horiz_closest...
                                            - length (alt_candidates));
            horiz_edge = -2;
          else
            horiz_closest = alt_candidates (horiz_closest);
            horiz_edge = -1;
          end
        end
      else
        horiz_closest = horiz_candidates(horiz_closest);
        horiz_edge = 2;
      end

      % Could decide here whether it really is horizontal or vertical.
      % For now, just take what it seems.
      if seems_horiz && ~isempty (horiz_closest)
        closest = horiz_closest(1);
        edge = horiz_edge;
      elseif ~isempty (vert_closest)
        closest = vert_closest(1);
        edge = vert_edge;
      else
        return
      end
      cl = rectangles(closest);

      add_undo ({'replace', closest, cl});

      if ismember ('control', modifiers)
        if (~isfield (rectangles, 'alt_days')) ...
            || isempty (rectangles(closest).alt_days)
              % TODO: reduce redundancy w.r.t. code below
          day = mod1 (round ((coordinates(1) + coordinates_now(1))/2), 7);
          rectangles(closest).alt_days = [day; mod1(day+1,7)];
          rectangles(closest).alt_hrs = rectangles(closest).on_off([2,4]);
        end
        if ~isempty (horiz_closest)
          closest = horiz_closest(1);
          cl = rectangles(closest);
        end
        A = 0;
        B = 0;
                % if seems moving a vertical edge
        if ~seems_horiz && ~isempty (cl.alt_hrs)
          A = between (cl.alt_hrs(1), coordinates(2),...
                       cl.on_off(2), size (cv, 2)/2);
          B = between (cl.alt_hrs(2), coordinates(2),...
                       cl.on_off(4), size (cv, 2)/2);
        end
        if A + B > 0
          edge = -3;
          % Find whether it is leading edge or trailing edge that is moved
          alt_edge = mod (coordinates(1), 7);
          if length (cl.alt_days) == 1
            alt_edge = 1 + (alt_edge > cl.alt_days);
          else
            % Take minimum modulo 7 distance
            dist = abs (alt_edge - cl.alt_days([1,end]));
            dist(dist > 3.5) = 7 - dist(dist > 3.5);
            [~, alt_edge] = min (dist);
          end
        else
          edge = horiz_edge;
          if horiz_edge > 0
            if isempty (cl.alt_days)
              % If this rectangle has no alt days, copy from an overlap
              for i = 1:horiz_candidates
                if ~isempty (rectangles(i).alt_days)
                  rectangles(closest).alt_days = rectangles(i).alt_days;
                  rectangles(closest).alt_hrs=cl.on_off([2,4]);
                  break
                end
              end
              % If no overlap has alt days either, take from mouse
              % TODO: reduce redundancy w.r.t. code above
              if isempty (cl.alt_days)
                day = mod1 (round ((coordinates(1) + coordinates_now(1))/2), 7);
                rectangles(closest).alt_days = [day; mod1(day+1,7)];
                rectangles(closest).alt_hrs = cl.on_off([2,4]);
              end
            end
            cl = rectangles(closest);
            dist = abs (cl.alt_hrs - coordinates(2));
            if length (dist) < 2 || dist(1) < dist(2)
              edge = -1;
            else
              edge = -2;
            end
          end
        end
      end

      % Move the edge by the amount of the mouse move
      if edge > 0       % a non-alt setting
        old_value = cl.on_off(edge);
        new = cl.on_off(edge) + move(1 + seems_horiz);
        if edge == 1 || edge == 2
          new = max (1, new);
        elseif edge == 3
          new = min (size (cv, 2) + 1, new);
        else
          new = min (size (cv, 1) + 1, new);
        end
        if edge == 1 && isfield (rectangles, 'drift') ...
                     && cl.drift ~= 0
          rectangles(closest).on_off([2,4]) = ...
            mod1 (cl.on_off([2,4]) + cl.drift * (new - cl.on_off(edge)),...
                  size (cv, 1));
        end
        rectangles(closest).on_off(edge) = new;
      else
        if edge > -3    % changing alt_hrs
          new_time = cl.alt_hrs(-edge) + move(1 + seems_horiz);
          [val, pos] = min (abs (new_time - cl.on_off([2,4])));
          if val < 1
            new_time = cl.on_off(2*pos);
          end
          rectangles(closest).alt_hrs(-edge) = new_time;
          cl.alt_hrs(-edge) = new_time;
                    % remove alt settings if they're not longer needed
          if all (cl.alt_hrs(:) == cl.on_off([2,4]))
            rectangles(closest).alt_hrs = [];
            rectangles(closest).alt_days = [];
            if all (cellfun (@isempty, {rectangles(:).alt_days}))
              rectangles = rmfield (rectangles, {'alt_days', 'alt_hrs'});
            end
          end
        else        % changing alt_days
          % TODO: Check if we should add/remove alt days, or shift them.
          alt_days = cl.alt_days;
          new_alt = alt_days + round(move(1));
          %[means, signs] = alt_day_means (cl, cv);
          done = false;
          if alt_edge == 1
            if move(1) < 0
              new_alt = [alt_days(1)-1; alt_days];
            elseif length (alt_days) > 1
              new_alt = alt_days(2:end);
            end
          else
            if move(1) > 0
              new_alt = [alt_days; alt_days(end)+1];
            elseif length (alt_days) > 1
              new_alt = alt_days(1:end-1);
            end
          end
          rectangles(closest).alt_days = mod1 (new_alt, 7);
        end
      end

      % After a horizontal move (changing start/end date), see if other
      % rectangles abutted this originally, and now overlap.
      % If so, adjust their edges too, unless the user is holding 'SHIFT'.
      if ~seems_horiz && edge > 0 % moving vertical edge -> horizontal move
        if ~ismember('shift',modifiers)
          cl = rectangles(closest);     % update, as on_off changed
          other_edge = 4 - edge;  % 1 -> 3 and 3 -> 1
          matches = find (abs (on_off (other_edge,:) - old_value) < 2);
          disp (matches);
          for i = matches(:)'
            if i == closest       % avoid issues with very narrow rectangles
              continue
            end
            if ~isempty (intersect (cl.burst, rectangles(i).burst))
              oo = rectangles(i).on_off;
              oo(other_edge) = cl.on_off(edge);
              if oo(3) > oo(1)    % if we haven't swapped left<->right
                add_undo ({'replace', i, rectangles(i)});
                rectangles(i).on_off(other_edge) = cl.on_off(edge);
              end
            end
          end
        end
      end

      draw_rectangles;
    end
  end

  function [days, confidences] = guess_alt_days (rect, cv)
    % Should be based on both on/off jumps and
    % differences between days.
    [means, signs] = alt_day_means (rect, cv);
    if all (means(:,1) == 0)
      means = means(:,2);
      signs = signs(:,2);
    elseif all (means(:,2) == 0)
      means = means(:,1);
      signs = signs(:,1);
    end
    jumps = diff ([means; means(1,:)]);
    [s, idx] = sort (abs (jumps));
    if all (idx(1,:) == idx(1,1) && idx(2,:) == idx(2,1))
      m = mean (s, 2);
      a = (m(1) > 0.5 * m(2) && m(1) < 2 * m(2));

      % Check sign is as expected
      % check magnitude is at least half of power
      % Check on period is at least 80% of power
      % What if it isn't?  Just reduce confidence?
    else
      if isfield (rect, 'alt_days')
        days = rect.alt_days;
        confidences = 0 * days;
      else
        days = [];
        confidences = 0;
      end
    end
  end

  function [means, signs] = alt_day_means (rect, cv)
    % return 7x2 array.
    % Top row is mean intensity of time between alt_hrs(1) and on_off(2)
    % on each of seven days.
    % Bottom row is the same for time between alt_hrs(2) and on_off(4).
    means = zeros (2, 7);
    signs = [0;0];
    if ~isfield (rect, 'alt_days')
      return;
    end
    % Find sets of rows to operate on
    diff_rows = cell{1,2};
    for i = 1:2
      a = rect.alt_hrs(i);
      b = rect.on_off(2*i);
      if a == b         % no alt-day behaviour on this edge
        continue;
      end
      if abs (a - b) > size (cv,1) / 2  % avoid wrap of midnight
        if a < b
          a = a + size (cv, 1);
        else
          b = b + size (cv, 1);
        end
      end
      signs(i) = a - b;
      tmp = max (a,b);
      a = min (a,b);
      b = tmp;

      % Omit partial rows, or take biggest if there is no non-partial row
      aa = ceil (a);
      bb = floor (b) - 1;
      if bb < aa
        if aa - a < bb+1 - b
          bb = aa;
        else
          aa = bb;
        end
      end

      aa = mod1 (aa:bb, size (cv, 1));
      vals = mean (cv (aa, rect.on_off(1):rect.on_off(3)-1), 1);
      vals = [nan*ones(mod(rect,on_off(1),7)), vals];
      len = length (vals);
      vals = [vals, nan*ones(7*(ceil(len/7)-len/7))];
      vals = reshape (vals, [7, length(vals)/7]);
      means(:,i) = mean (vals, 2);
    end
  end

  function add (~, ~)
    set (f, 'WindowButtonUpFcn', @add_rectangle);
  end

  function [st, en] = align_with_existing (down_pos, up_pos)
    % Find matching start and end times from existing rectangles
    if isempty (rectangles)
      match = [];
    else
      on_off = [rectangles(:).on_off];
      match = find (abs (on_off(3,:) - down_pos(1)) < 2);
    end
    if ~isempty (match)
      st = on_off(3, match(1));
    else
      st = round (down_pos(1));
    end
    if ~isempty (rectangles)
      match = find (abs (on_off(1,:) - up_pos(1)) < 2);
    end
    if ~isempty (match)
      en = on_off(1, match(1));
    else
      en = round (up_pos (1));
    end
  end


  function add_rectangle (~,~)
    coordinates_now = get (ha, 'CurrentPoint');
    coordinates_now = min (coordinates_now(1,1:2),[size(cv,2) size(cv,1)]+0.5);
    coordinates_now = max (coordinates_now, [1, 0.5]);

    [st, en] = align_with_existing (coordinates, coordinates_now);

    % even if drawn right-to-left, we need st < en
    % If drawn bottom to top, it is a rectange than spans midnight.
    if st > en
      [st, en] = deal (en, st);
    end

    % 0.5 and ceil (not round) as rectangle specifies midpoint of pixel
    on_off = [ceil(st); coordinates(2)+0.5; ceil(en); coordinates_now(2)+0.5];

    do_add_rectangle (on_off);

    % If we overlap another significantly, break it
    oo = [rectangles(:).on_off];
    overlap = on_off(3) <= oo(3,:) & on_off(1) >= oo(1,:) ...
              & on_off(2) <= oo(4,:) & on_off(4) >= oo(2,:);
    overlap = overlap(1:end-1);     % omit newly added rectangle
    if any (overlap)
      overlap = find (overlap);
      if length (overlap) > 1
        overlap = overlap(1);   % TODO: Find best match
      end

      drift = 0;
      if isfield (rectangles, 'drift')
        drift = rectangles(overlap).drift;
      end

      % Split to new rectangle, if this is in the middle
      if on_off(3) < rectangles(overlap).on_off(3) ...
          && on_off(1) > rectangles(overlap).on_off(1)
        new_rectangle = rectangles(overlap);
        new_rectangle.on_off([2,4]) = new_rectangle.on_off([2,4]) ...
                           + drift * (on_off(3) - new_rectangle.on_off(1));
        new_rectangle.on_off(1) = on_off(3);
        new_rectangle.on_off = tweak_hrs (new_rectangle, cv, 2);
        new_rectangle.on_off = tweak_hrs (new_rectangle, cv, 4);
        new_rectangle.missed = [];    % TODO: recalculate, or copy

        rectangles = [rectangles(:); new_rectangle];
        add_undo ({'delete', length(rectangles), []});
      end

      add_undo ({ 'replace', overlap, rectangles(overlap) });
      if isempty (rectangles(overlap).missed)
        rectangles(overlap).missed = false (1, diff(rectangles(overlap).on_off([1, 3])));
      end

      if on_off(1) > rectangles(overlap).on_off(1)
        rectangles(overlap).missed = rectangles(overlap)...
                  .missed(1:on_off(1) - rectangles(overlap).on_off(3)+1);
        rectangles(overlap).on_off(3) = on_off(1);
      else
        rectangles(overlap).missed = rectangles(overlap)...
                  .missed(on_off(3) - rectangles(overlap).on_off(1)+1:end);
        rectangles(overlap).on_off([2,4]) = rectangles(overlap).on_off([2,4]) ...
                                      + drift * (on_off(3) - on_off(1));
        rectangles(overlap).on_off(1) = on_off(3);
        rectangles(overlap).on_off = tweak_hrs (rectangles(overlap), cv, 2);
        rectangles(overlap).on_off = tweak_hrs (rectangles(overlap), cv, 4);
      end

      r = other_chain (on_off, oo, cv);
      if ~isempty (r)

      end
    end

    % Next time, go back to nudging rectangles
    set (f, 'WindowButtonUpFcn', @final_select_or_move);
    draw_rectangles;
  end

  function do_add_rectangle (on_off)
    new_rectangle.on_off = [];           %initialize
    for fn = fieldnames (rectangles)'
      new_rectangle.(fn{1}) = [];
    end
    if isfield (new_rectangle, 'drift')
      if ~isempty (rectangles)
        new_rectangle.drift = median ([rectangles.drift]);
        new_rectangle.power = mean ([rectangles(:).power]);
      else
        new_rectangle.drift = 0;
        new_rectangle.power = 0.4;    % TODO: estimate from on/off jumps
      end
    end

    new_rectangle.on_off = on_off;
    new_rectangle.missed = [];
    if on_off(2) < on_off(4)
      new_rectangle.burst = ceil(on_off(2)):floor(on_off(4));
    else
      new_rectangle.burst = [ceil(on_off(2)):48, 1:floor(on_off(4))];
    end
    new_rectangle.burst = round (new_rectangle.burst);
    new_rectangle.trust = mean ([rectangles(:).trust], 2);

    rectangles = [rectangles(:); new_rectangle];

    add_undo ({'delete', length(rectangles), []});
  end

  function toggle_alt_days (object_handle,~)
    if ~isempty (selected_rect)
      if isempty (rectangles (selected_rect).alt_days)
        days_to_hide = [];
      else
        days_to_hide = setdiff (1:7, rectangles (selected_rect).alt_days);
      end
    else
      modifiers = get_modifiers (object_handle);
      if ismember ('shift', modifiers)
        days_to_hide = [];
      else
        days_to_hide = setdiff (1:7, days_to_hide);
      end
    end
    selected_rect = [];
    confidence_global = [];
    draw_rectangles;
  end

  function add_undo (action)
    undo_stack{end + 1} = action;
    redo_stack = {};
  end

  function undo (~, ~)
    if isempty (undo_stack)
      return
    end

    % If one rectangle is selected, and we can identify it,
    % undo the most recent change to it, instead of the most recent of all
    found_match = false;
    if length (selected_rect) == 1
      for i = length (undo_stack):-1:1
        undo_action = undo_stack{i};
        if undo_action{2} == selected_rect
          if length (undo_action) == 3 && ~isempty (undo_action{3})
            old_rect = undo_action{3};
            if sum (old_rect.on_off == rectangles(selected_rect).on_off) >= 3
              found_match = true;
              break;
            end
          else
            found_match = true;
            break;
          end
        end
      end
      if ~found_match
        selected_rect = [];
        draw_rectangles;
        return
      end
    end

    if found_match
      undo_stack = undo_stack([1:i-1, i+1:end]);
    else
      undo_action = undo_stack{end};
      undo_stack(end) = [];
    end

    r = undo_action{2};
    if r <= length (rectangles)
      rec = rectangles(r);
    else
      rec = [];
    end
    redo_stack{end + 1} = {undo_action{1}, r, rec};
    if isstruct (undo_action{3}) && ~isequal (fields (rectangles), fields (undo_action{3}))
      rectangles = rmfield (rectangles, setdiff (fields (rectangles), fields (undo_action{3})));
      for i = setdiff (fields (undo_action{3}), fields (rectangles))'
        empty = [];
        if strcmp (i, 'drift')
          empty = 0;
        end
        rectangles(1).(i{1}) = empty;
      end
    end
    switch undo_action{1}
      case 'replace'
        try
          rectangles (r) = undo_action{3};
        catch
          keyboard
        end
      case 'delete'
        rectangles (r) = [];
      case 'insert'
        rectangles(r+1:end+1) = rectangles(r:end);
        rectangles(r) = undo_action{3};
    end
    draw_rectangles;
  end

  function redo (~, ~)
    if isempty (redo_stack)
      return
    end
    redo_action = redo_stack{end};
    redo_stack(end) = [];
    r = redo_action{2};
    if r <= length (rectangles)
      rec = rectangles(r);
    else
      rec = [];
    end
    undo_stack{end + 1} = {redo_action{1}, r, rec};
    switch redo_action{1}
      case 'replace'
        rectangles (r) = redo_action{3};
      case 'delete'
        rectangles(r+1:end+1) = rectangles(r:end);
        rectangles(r) = redo_action{3};
      case 'insert'
        rectangles(r) = [];
    end
    draw_rectangles;
  end

  function key_press (object_handle, event_data)
    switch event_data.Key
      case 'delete'
        if isempty (selected_rect)
          select_rectangle (object_handle, event_data)
        end
        if ~isempty (selected_rect)
          del
        end
      case 'slash'
        if ~isempty(selected_rect)
          confidence = get (handle_conf, 'Value');
          for i = selected_rect
            rectangles(i).colour = [0.05; 0.05; 0.05] * (10+confidence);
            rectangles(i).manual_confidence = confidence;
          end
          selected_rect = [];
          draw_rectangles;
        end
      case 'm'
        merge_selected (object_handle, event_data);
      case 'u'
        undo (object_handle, event_data);
      case 'y'
        redo (object_handle, event_data);
      case 'z'
        undo (object_handle, event_data);
    end
  end

  function affected_rects = auto_fix
    % Guess what could be wrong.
    %  - close a gap
    %  - make timing match the other chain
    %  - missed-days
    %  - copy from other chain
    if coordinates(1,1) > size (cv, 2) || coordinates(1,2) > size (cv, 1)
      return
    end

    oo = [rectangles.on_off];

    if ~isempty (selected_rect)
      r = selected_rect(1);
      rr = rectangles(r);
      selected_rect = [];
      confidence_global = [];

      % Try to re-adjust nearest edge
      % 1. choose nearest edge
      rr.on_off([2,4]) = mod1 (rr.on_off([2,4]), size (cv, 1));
      left_right = (rr.on_off(3) - coordinates(1,1)) ...
                    / (rr.on_off(3) - rr.on_off(1));
      if rr.on_off(2) < rr.on_off(4)
        top_bottom = (rr.on_off(4) - coordinates(1,2)) ...
                   / (rr.on_off(4) - rr.on_off(2));
      else
        % > 0.5 means top, < 0.5 means bottom
        if coordinates(1,2) < 0.5 * (rr.on_off(2) + rr.on_off(4))
          top_bottom = min (0.5, abs (rr.on_off(4) - coordinates(1,2)) ...
                            / (rr.on_off(4) - rr.on_off(2) + size (cv,1)));
        else
          top_bottom = 1 - min (0.5, abs (coordinates(1,2)- rr.on_off(2)) ...
                            / (rr.on_off(4) - rr.on_off(2) + size (cv,1)));
        end
      end

      % 2. If change vertical edge
      if abs (left_right - 0.5) > abs (top_bottom - 0.5)
          % TODO: should handle alt_days
        % First, shrink if the first/last day has many below-power values
        if ~isfield (rr, 'alt_days') || isempty (rr.alt_days)
          b = rr.burst (1:max (1, end-1));
          pwr = 0.8 * rr.power;
          st = rr.on_off(1);
          en = rr.on_off(3);
          days = st:en-1;

          if rr.on_off(2) < rr.on_off(4)
            patch = cv(b, days);
          else
            patch=[cv(ceil (rr.on_off(2)):end, days);
                   cv(1:floor (rr.on_off(4))-1, min (days+1,size (cv,2)))];
          end
          num_small = sum (patch < 0.9 * pwr, 1);

          done = false;
          if left_right > 0.5
            pos = 1;
            factor = 1;
          else
            pos = 3;
            factor = -1;
            num_small = num_small(end:-1:1);
          end
          for i = 1:min (length (num_small), 4)
            if num_small(i) < i
              break;
            end
            done = true;
          end
          if done
            add_undo ({ 'replace', r, rectangles(r) });
            rectangles(r).on_off(pos) = rr.on_off(pos) + factor * (i-1);
            return;
          end
        end

        % Extend to close a gap
        time_overlap = (between (oo(2,:), rr.on_off(4), oo(4,:),size (cv, 2)/2)...
                      | between (oo(2,:), rr.on_off(2), oo(4,:),size (cv, 2)/2) ...
                      | between (rr.on_off(2), oo(2,:), rr.on_off(4), size (cv, 2)/2) ...
                      | between (rr.on_off(2), oo(4,:), rr.on_off(4), size (cv, 2)/2) ...
                      | rr.on_off(2) == oo(2,:) | rr.on_off(4) == oo(4,:));
        %  First check right
        if coordinates(1) > rr.on_off(3) - 0.25 * (rr.on_off(3) - rr.on_off(1))
          % check if there is another rectangle slightly further to the right
          near = time_overlap & oo(1,:) > rr.on_off(3) - 10 ...
                              & oo(1,:) < rr.on_off(3) + 10;
          near(r) = false;
          near = find (near);

          if ~isempty (near)
            add_undo ({ 'replace', r, rectangles(r) });
            n = oo(2, near) < rr.on_off(4) & oo(4, near) > rr.on_off(2);
            if any (n)
                  % simple case, not straddling midnight
              rectangles(r).on_off(3) = min (oo(1,near(n)));
            elseif oo(2, near(1)) > rr.on_off(4)
                  % add 1 if only one straddles midnight, or 0 if both do
              rectangles(r).on_off(3) = min (oo(1,near)) ...
                                        + any (oo(4, near) > rr.on_off(2));
            else
              rectangles(r).on_off(3) = min (oo(1,near)) - 1;
            end
            affected_rects = [r, near];
            return
          end
        %  Next check left
        elseif coordinates(1) < rr.on_off(1) + 0.25 * (rr.on_off(3) - rr.on_off(1))
          % check if there is another rectangle slightly further to the right
          near = time_overlap & oo(3,:) < rr.on_off(1) + 10 ...
                              & oo(3,:) > rr.on_off(1) - 10;
          near(r) = false;
          near = find (near);

          if ~isempty (near)
            add_undo ({ 'replace', r, rectangles(r) });
            n = oo(2, near) < rr.on_off(4) & oo(4, near) > rr.on_off(2);
            if any (n)
                  % simple case, not straddling midnight
              rectangles(r).on_off(1) = max (oo(3,near(n)));
            elseif oo(2, near(1)) > rr.on_off(4)
              rectangles(r).on_off(1) = max (oo(3,near)) ...
                                        + any (oo(4, near) > rr.on_off(2));
            else
              rectangles(r).on_off(1) = max (oo(3,near)) - 1;
            end
            affected_rects = [r, near];
            return
          end
        end

        % big adjustment of start or end day
        [st, en] = timed_days (rectangles, r, cv, ...
                               [left_right>0.5, left_right<0.5]);
        if left_right > 0.5   % adjust left edge
          add_undo ({'replace', r, rectangles(r)});
          rectangles(r).on_off(1) = st;
        else
          add_undo ({'replace', r, rectangles(r)});
          rectangles(r).on_off(3) = en;
        end
      else
        % adjust start or end hour
        if top_bottom > 0.5
          edge = 2;
        else
          edge = 4;
        end
        [on_off, unchanged] = tweak_hrs (rectangles(r), cv, edge);
        if unchanged
          % look for big change
          centre = (rr.on_off([1,2]) + rr.on_off([3,4])) / 2;
          if rr.on_off(4) < rr.on_off(2)
            centre(2) = mod1 (centre(2) + size (cv, 1)/2, size (cv, 1));
          end
          % Find a box around centre.
          % It would be better to fix the sides, but this should do
          on_off = find_region (cv, rr.on_off, rr.power, centre);
          rr.on_off(edge) = on_off(edge);
          on_off = tweak_hrs (rr, cv, edge);
        end
        add_undo ({'replace', r, rectangles(r)});
        rectangles(r).on_off = on_off;
      end
      return
    end

    if isempty (oo)
      [h, bins] = hist (cv(:), 20);
      peak = sort (secondPeak (h(:)', 2));
      bp = bins(peak);
      power = bp(2) - bp(1);

      left = 1;
      right = size (cv, 2);
      a = find (cv(round(coordinates(1,2)), round(coordinates(1,1)):-1:left) < power, 1);
      if ~isempty (a)
        left = round (coordinates(1,1)) - a + 2;  % +2: First *of* rectangle
      end
      a = find (cv(round(coordinates(1,2)), round(coordinates(1,1):right)) < power, 1);
      if ~isempty (a)
        right = round(coordinates(1,1)) + a - 1;  % -1: First after rectangle
      end
      pair_on_off = [left; 1; right; 1];

      on_off = find_region (cv, pair_on_off, power, ...
                      coordinates(1,1:2), [], rectangles);
      if ~isempty (on_off)
        do_add_rectangle (on_off);
        rectangles(end).power = power;
        rectangles(end).on_off = tweak_hrs (rectangles(end), cv, 2);
        rectangles(end).on_off = tweak_hrs (rectangles(end), cv, 4);
      end
      return;
    end

    % neighbours = find (between (oo(2, :), coordinates(1,2), ...
    %                             oo(4, :), size (cv, 1) / 2));
    left  = max (oo(3, oo(3,:) < coordinates(1)));
    right = min (oo(1, oo(1,:) > coordinates(1)));

    if isempty (left)
      left = 1;
    end
    if isempty (right)
      right = size (cv, 2);
    end

    a = find (cv(round(coordinates(1,2)), round(coordinates(1,1)):-1:left) < mean ([rectangles.power]), 1);
    if ~isempty (a)
      left = round (coordinates(1,1)) - a + 2;  % +2: First *of* rectangle
    end
    a = find (cv(round(coordinates(1,2)), round(coordinates(1,1):right)) < mean ([rectangles.power]), 1);
    if ~isempty (a)
      right = round(coordinates(1,1)) + a - 1;  % -1: First after rectangle
    end

    same_day = find (oo(1,:) <= coordinates(1) & coordinates(1) <= oo(3,:));
    if length (same_day) > 1
      same_day (between (oo(2,same_day), coordinates(2), ...
                         oo(4,same_day), size (cv, 1) / 2)) = [];
    end
    if length (same_day) == 1
      pair_on_off = rectangles(same_day).on_off;
      pair_on_off(1) = max (pair_on_off(1), left);
      pair_on_off(3) = min (pair_on_off(3), right);
      on_off = find_region (cv, pair_on_off, mean ([rectangles.power]), ...
                            coordinates(1,1:2), [], rectangles);
      if ~isempty (on_off)
        add_or_alt_day (on_off);
        return
      end
    end

    if right > left
      on_off = find_region (cv, [left; 1; right; 1], mean ([rectangles.power]), ...
                            coordinates(1,1:2), [], rectangles);
      if ~isempty (on_off)
        add_or_alt_day (on_off);
        return
      end
    end
  end

  function [on_off, unchanged] = tweak_hrs (rect, cv, edge)
    % adjust on_off([2,4]) by at most one hour
    on_off = rect.on_off;
    if isfield (rectangles, 'drift')
      drift = sum (rect.drift);
    else
      drift = 0;
    end
    on_off = rect.on_off;
    if floor (on_off(edge)) ...
        == floor (on_off(edge) + drift * (on_off(3) - on_off(2)))
      row = mod1 (floor (min (on_off(edge) + [0 drift])) + [-1 0 1], size (cv,1));
      % [pre, mid, post] = edge_neighbours (on_off(edge) + drift, size (cv, 1));
      pr = cv (row(1), on_off(1):on_off(3)-1);
      md = cv (row(2), on_off(1):on_off(3)-1);
      po = cv (row(3),on_off(1):on_off(3)-1);
      m = median ((po - md) ./ (po - pr), 'omitnan');
      candidate = floor (on_off(edge)) + m;
      if on_off(edge) < floor (on_off(edge)) + 0.5
        offset = -1;
      else
        offset = 1;
      end
      old_jump = mean (po - pr);

      row = mod1 (row + offset, size (cv, 1));
      pr = cv (row(1), on_off(1):on_off(3)-1);
      md = cv (row(2), on_off(1):on_off(3)-1);
      po = cv (row(3),on_off(1):on_off(3)-1);
      new_jump = mean (po - pr);
      if edge == 4
        old_jump = -old_jump;
        new_jump = -new_jump;
      end
      p = rect.power;
      if (on_off(edge) == floor (on_off(edge)) || old_jump < 1.1 * p) ...
          && abs (new_jump - p) < abs (old_jump - p)
        m = median ((po - md) ./ (po - pr), 'omitnan');
        candidate = row(2) + m;
      end
      candidate = candidate - drift * (on_off(3) - on_off(1))/2;
      if rect.on_off(edge) == candidate
        unchanged = [min(abs ([new_jump old_jump] / p - 1)), ...
                     max(abs ([new_jump old_jump]))];
      else
        on_off(edge) = candidate;
        unchanged = false;
      end
    else
      % split over multiple rows
      unchanged = true;
    end
  end

  function add_or_alt_day (on_off)
    % either add a rectangle with on_off,
    % or set set the alt_day of a nearby rectangle if appropriate
        if on_off(3) - on_off(1) ~= 2
          do_add_rectangle (on_off);
        else
                  % Check for alt_days
          r = set_alt_day (on_off);
          if isempty (r)
            do_add_rectangle (on_off);
          end
        end
  end

  function r = set_alt_day (on_off)
    oo = [rectangles.on_off];
    d = find (on_off(3) > oo(1,:) & on_off(1) < oo(3,:));
    time_overlap = (between (oo(2,d), on_off(4), oo(4,d),size (cv, 2)/2)...
                  | between (oo(2,d), on_off(2), oo(4,d),size (cv, 2)/2) ...
                  | between (on_off(2), oo(2,d), on_off(4), size (cv, 2)/2) ...
                  | between (on_off(2), oo(4,d), on_off(4), size (cv, 2)/2) ...
                  | on_off(2) == oo(2,d) | on_off(4) == oo(4,d));
    r = find (time_overlap, 1);
    if isempty (r)
      return;
    end
    r = d(r);
    % Choose alt days
    alt_days = false(1, 7);
    alt_days([mod1(on_off(1), 7), mod1(on_off(1)+1,7)]) = true;

    on_off_2 = oo(2,r);
    on_off_4 = oo(4,r);
    known_on = round ((on_off_4 + on_off_2) / 2);
    if on_off_4 < on_off_2
      known_on = known_on + size (cv, 1) / 2;
    end

    range = oo(1,r):oo(3,r)-1;
    power = rectangles(r).power;

    [~, time_on,  conf1] = alt_jump (alt_days, [], on_off_2, known_on, -power, range, cv);
    [~, time_off, conf2] = alt_jump (alt_days, [], on_off_4, known_on,  power, range, cv);
    if conf1 < 0.5 && conf2 < 0.5
      r = [];
      return;
    end

    %wasnt_alt = ~isfield (rectangles, 'alt_days');
    add_undo ({ 'replace', r, rectangles(r) });
    rectangles(r).alt_days = find (alt_days);
    rectangles(r).alt_hrs = [time_on; time_off];
    if conf1 < 0.5 || hr_diff (on_off_2, rectangles(r).alt_hrs(1), size (cv, 1)) < 0.5
      rectangles(r).alt_hrs(1) = on_off_2;
    end
    if conf2 < 0.5 || hr_diff (on_off_4, rectangles(r).alt_hrs(2), size (cv, 1)) < 0.5
      rectangles(r).alt_hrs(2) = on_off_4;
    end

    % Find times on non-alt days

    % Check other rectangles
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

  function sanitise
    if ~isempty (rectangles)
      % restore sanity
      oo = [rectangles.on_off];

      % Trim to valid size
      oo(1:2,:) = max (oo(1:2,:), 1);
      oo([1,3],:) = min (oo([1,3],:), size (cv, 2) + 1);
      oo([2,4],:) = min (oo([2,4],:), size (cv, 1) + 1);
      oo([1,3],:) = sort (round (oo([1,3], :)));

      % Omit zero-width rectangles
      idx = (oo(1,:) ~= oo(3,:));
      oo = oo(:, idx);
      rectangles = rectangles(idx);

      % Put back in rectangles
      for i = 1:length (rectangles)
        rectangles(i).on_off = oo(:, i);
      end
      if isfield (rectangles, 'drift')
        for i = 1:length (rectangles)
          if ~isfinite (rectangles(i).drift)
            rectangles(i).drift = 0;
          end
        end
      end
      if isfield (rectangles, 'alt_days')
        for i = 1:length (rectangles)
          rectangles(i).alt_days = rectangles(i).alt_days(:);
          rectangles(i).alt_hrs  = rectangles(i).alt_hrs(:);
        end
      end
    end
  end

  function done (~, ~)
    sanitise;
    confidence = get (handle_conf, 'Value');
    if confidence ~= 10
      explicit_confidence = [];
      if ~isempty (rectangles)
        if ~isfield (rectangles, 'manual_confidence')
          rectangles(1).manual_confidence = [];
        else
          explicit_confidence = [rectangles.manual_confidence];
        end
      end

      if all (explicit_confidence ~= confidence) ...
          || ~isempty (confidence_global)  % if changed since last use
        for i = 1:length (rectangles)
          if isempty (rectangles(i).manual_confidence)
            rectangles(i).manual_confidence = confidence;
          end
        end
      end
    end
    delete (f);
    f = 0;
  end

  function cancel (object_handle, ~)
    ok = false;
    modifiers = get_modifiers (object_handle);
    if ismember('control',modifiers)
      if ~isempty (rectangles)
        rectangles(1).uncertain = true;
      end
      ok = true;
    end
    done;
  end

  function set_next_house (object_handle, ~)
    next = str2double (get (object_handle, 'String'));
    cancel (object_handle);
  end

end

% New functions
function pair = other_chain (me, oo, cv)
  % Find the number of the rectangle whose on_off is in  oo,  and
  % is the "match" for  me  in the other "chain".
  % Returns [] if none.
  hr = size (cv, 1) / 2;
  pair = find ((me(1) < oo(3,:) & me(3) > oo(1,:)) ...
               & ~between (oo(2,:), me(2), oo(4,:), hr) ...
               & ~between (oo(2,:), me(4), oo(4,:), hr));
  if length (pair) > 1
    pair = pair(1);   % TODO: Find "best" match
  end
end

function on_off = find_region (cv, pair_on_off, power, same, ~, ~) % Missing are 'other, rectangles'
  on_off = [];
  % Find rectangle matching one in the other chain
  st = pair_on_off(1);
  en = pair_on_off(3);
  % Find neighbours
  % Find jumps of top and bottom from neighbours.
  % Find would-be neighbours in other chain.
  % Find symmetric tops and bottoms for candidate

  pad = floor (size (cv, 1) / 2);
  m = min (cv (:, st:max (st, en-1)), [], 2);
  m_pre = min (cv(end-pad+1:end, max (1, st-1):max ([1, st-1 en-2])), [], 2);
  m_post= min (cv(1:pad, min (size (cv, 2), st+1:max(en, st+1))), [], 2);
  mm = [m_pre; m; m_post];
  r = ranges (find (mm > power));
  if isempty (r)
    r = ranges (find (mm > 0.75 * power));
    if isempty (r)
      return
    end
  end

  % If passed co-ordinates instead of a rectangle index,
  % find rectangle containing the point
  if length (same) == 2
    r = r - pad;    % cancel padding of mm
%    r = r (:, r(1,:) > 0 & r(2,:) <= length (m));

    holds_me = find (between (r(1,:)-0.5, same(2), r(2,:)+0.5, size (cv, 1) / 2));
    if isempty (holds_me)
      [gap, holds_me] = min (abs ([r(1,:)-same(2), r(2,:)-same(2)]));
      if gap > 2
        holds_me = [];
      elseif holds_me > size (r, 2)
        holds_me = holds_me - size (r, 2);
      end
    end
    if length (holds_me) == 1
      r3 = r + pad;
      step1 = mm(r3(1, holds_me))     - mm(r3(1, holds_me) - 1);
      step2 = mm(r3(1, holds_me) - 1) - mm(r3(1, holds_me) - 2);
      stepA = step1;
      if abs (step1 - power) > abs (step2 - power)
        r(1, holds_me) = 1 + mod (r(1, holds_me) - 2, size (cv, 1));
        stepA = step2;
      end

      step1 = mm(r3(2, holds_me) - 1) - mm(r3(2, holds_me));
      step2 = mm(r3(2, holds_me))     - mm(r3(2, holds_me) + 1);
      stepB = step1;
      if abs (step1 - power) > abs (step2 - power)
        r(2, holds_me) = 1 + mod (r(2, holds_me), size (cv, 1));
        stepB = step2;
      end

      % TODO: do biggest improvement first.
      if stepA < 1.1 * power
        [j, best_jump] = min (abs (diff (mm(r3(1,holds_me):r3(2,holds_me))) - power));
        if j < abs (power - stepB)
          r(1, holds_me) = 1 + mod (r(1, holds_me) + best_jump - 1, ...
                                    size (cv, 1));
        end
      end
      if stepB < 1.1 * power
        [j, best_jump] = min (abs (diff (mm(r3(1,holds_me):r3(2,holds_me))) + power));
        if j < abs (power - stepB)
          r(2, holds_me) = 1 + mod (r(1, holds_me) + best_jump - 1, ...
                                    size (cv, 1));
        end
      end

      r(:, holds_me) = mod1 (r (:, holds_me), size (cv, 1));
      on_off = [st; r(1, holds_me); en; r(2, holds_me)];
      return;
    end
  end

  %{
  short = r(:, r(1,:)+2 > r(2,:));
  r = r(:, r(1,:)+2 <= r(2,:));   % only consider rectangles of length > 2

  r(r <= 3) = r(r <= 3) + size (cv, 1);
  idx = (r > size (cv,1) + 3);
  r(idx) = r(idx) - size (cv, 1);
  if size(r,2) > 1 && all (r(:, end) == r(:,1))
    r = r(:, 1:end-1);
  end

  jp_up = max ([mm(r(1,:)), mm(r(1,:)+1), mm(r(1,:)+2)]') - min ([mm(r(1,:)-1), mm(r(1,:)-2)]');
  jp_dn = min ([mm(r(2,:)+1), mm(r(2,:)+2)]') - max ([mm(r(2,:)), mm(r(2,:)+1), mm(r(2,:)+2)]');
  jp_up = max (jp_up, 1e-10);
  jp_dn = min (jp_dn, 0);
  times = [r(1,:), r(2,:)] - 3;
  jp = [jp_up, jp_dn];
  u_jp = abs (jp);

  score = power ./ (0.2 * power + wabs (u_jp - power, -0.2));
  up = (jp > 0);
  score(up)  = score(up)  + 10 ./ (10 + abs (times(up)  - rectangles(i).on_off(2)));
  score(~up) = score(~up) + 10 ./ (10 + abs (times(~up) - rectangles(i).on_off(4)));
  pairs = bsxfun (@plus, score(up), score(~up)');
  pairs = pairs + 0.5 ./ (0.5 + bsxfun (@minus, times(~up)', times(up) + rectangles(i).on_off(4)-rectangles(i).on_off(2)));
  pairs = pairs + power ./ (0.2 * power + abs (bsxfun (@minus, u_jp(~up)', u_jp(up))));
  down = find (~up);
  up = find (up);

  oo = [rectangles.on_off];
  olap = [oo(:, oo(3,:) >= en-2 & oo(1,:) <= st-2 & i ~= 1:size(oo,2))];
  for n = 1:length (up)
    for j = 1:length (down)
      burst = ceil (times(up(n))):floor (times(down(j)));
      if isempty (burst)
        burst = [ceil(times(up(n))):size(cv,1), 1:floor(times(down(j)))];
      end
      if any (m(burst) < power)
        pairs(n,j) = 0;
      end
      % Skip if we substantially overlap an existing rectangle
      if any (max (abs ([burst(1)-olap(2,:) burst(end)-olap(4,:)+1])) <= 1)
        pairs(n,j) = 0;
      end
    end
  end
  %}
end

% Functions copied from has_pool_pump.m
function m = mod1 (value, modulus)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Modulo operator, returning values in [1, modulus], not [0, modulus-1].
  m = mod (value-1, modulus) + 1;
end

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

function [days, time, confidence] = alt_jump (no_jump, DoW, main_time, known_on, power, range, cv)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Find turn on/off time on days that don't match the standard time
 % no_jump = logical vector of which days don't have a big step at main_time
 % DoW = vector of mean jump size on each day of week
 % main_time = time of current estimate of jump
 % known_on = a time that the device is "known" to be on
 % range = list of days in this rectangle
 % Days are 1 = 1st,8th... day of valid_days, 7 = 7th,14th... day of valid_days
  confidence = 0;
  reduce_confidence_factor = 1; % If we "nudge" a choice, reduce this to (0,1)
  time = 0;
  days = find (no_jump);
  n_alt = length (days);
  if n_alt > 0
    if n_alt == 1
      tmp = DoW;
      tmp(no_jump) = max (DoW);
      [~, next_smallest] = min (tmp);
      if hr_diff (days, next_smallest, 7) == 1
        reduce_confidence_factor = 0.7;
        days = [days, next_smallest];
        n_alt = 2;
      end
    end
    match = zeros (size (days));
    for j = 1:length (days)
      try
        x = mean (cv (:, range(days(j):7:end)), 2);
      catch
        keyboard
      end
      match(j) = match_jump (x, -known_on, power);
    end
    % if 2 days with no jump at off_time have "better" jumps at similar times
    % to each other, but different from off_time, be probably alt timer
    if n_alt == 2 && any (diff (days) == [-1, -6, 1, 6])
      if hr_diff (match(1), match (2), size(cv,1)) < 2 ...
         &&  hr_diff ((match(1)+match(2))/2, main_time, size(cv,1)) > 1
          confidence = 1;
      else
        x = mean(cv(:,range([days(1):7:end,days(2):7:end])), 2);
        both_match = match_jump (x, -known_on, power);
        if min (abs (match - both_match)) < 1
          confidence = 0.5;
        end
      end
    % What other cases suggest alternate timer?
    else
      [~, most] = min (DoW + DoW ([2:end, 1]));
      diffs = hr_diff (match, main_time, size(cv,1));
      if sum (diffs > 1.5) == 2
        days = days(diffs > 1.5);
        if all (diff (days) ~= [-1, -6, 1, 6])
          days = [];
        elseif any (most == days) && any (mod1 (most+1, 7) == days)
          confidence = 0.9;
        else
          confidence = 0.5;
        end
      else
        m1 = find (days == most);
        m2 = find (days == mod1 (most+1, 7));
        if ~isempty (m1) && ~isempty (m2)
          if hr_diff (match (m1), match (m2), size (cv,1)) < 2
            days = [most, mod1(most+1,7)];
            if hr_diff ((match(1)+match(2))/2, main_time, size(cv,1)) > 1
              confidence = 0.9;
              % TODO -- find which other days also match
            end
          else
            %keyboard
            days = [];
          end
        else
          days = [];
          %keyboard
        end
      end
    end

    if ~isempty (days)
      % Match difference between normal and alt days
      % 1. difference of means
      y = 0;
      x = 0;
      for d = 1:min (7, length (range))
        if any (mod (range(d) - days, 7) == 0)
          x = x + mean (cv (:, range(d:7:end)), 2);
        else
          y = y + mean (cv (:, range(d:7:end)), 2);
        end
      end
      x = x / length (days);
      y = y / (7- length (days));

      % 2. First guess, to see whether alt times are earlier or later
      time = match_jump (x, -known_on, power);

      % 3. Match, if time seems to be different
      if hr_diff (time, main_time, size(cv,1)) > 1
        if (power < 0 && inorder (time, main_time, known_on)) ...
            || (power > 0 && inorder (known_on, main_time, time))
          time = match_jump (x - y, round (main_time), power);
        else
          time = match_jump (y - x, round (main_time), -power);
        end
        % TODO: If large difference between nominal and alt times, reduce conf
        confidence = confidence * reduce_confidence_factor;
      end
      days = mod1 (days + range(1) - 1, 7);
    end
  end
end
