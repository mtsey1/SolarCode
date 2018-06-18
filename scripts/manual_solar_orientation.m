function [az, ze, cap, abort, orientation_data, jump_to_unsaved, jump_to_house] = manual_solar_orientation (az, ze, cap, seen, sunPos, data, data_no_vamp, s, hot_days, cold_days, solar_start, solar_end, solar_range, meta, data_idx, orientation_data_cell, disconnected, irradiation, funcs)
  orientation_data = orientation_data_cell.orientation_data;
  connected = find (~disconnected);

  dims = get (0, 'ScreenSize');
  f = figure ('Visible', 'off', 'Position', [0, 0.5*dims(4), 0.8*dims(3), 0.8*dims(4)]);

  ha = axes ('Units', 'pixels', 'Units', 'normalized', ...
             'Position', [0.05, 0.05, 0.7, 0.9]);

  if(~isnan(data_idx))
%       if((size(orientation_data_cell.orientation_data,1) >= data_idx) && ~isequal(orientation_data_cell.orientation_data(data_idx,:),[0 0 0 0]))
      if((length(orientation_data_cell.orientation_data) >= data_idx) && ~isempty(orientation_data_cell.orientation_data{data_idx}))

          current_saved_flag = true;

      else

          current_saved_flag = false;

      end
  end

  jump_to_unsaved = false;
  jump_to_house = nan;

  included_mask = true(365,1);

  default_confidence = 90;
  confidence = default_confidence;
  az_locked_flag = false;
  ze_locked_flag = false;
  abort = true;
  undo_stack = {};
  redo_stack = {};
  curr_view = 'compare';
  default_tags = {'unspecified'};
  tags = default_tags;
  comments = '';

  % Load tag file if exists
%   tags = default_tags;

  tag_group_file = '../data/saved_tag_groups.mat';
  if(exist(tag_group_file,'file'))
      load(tag_group_file);
  end

  tag_value = 1;
  if(current_saved_flag)
      [ismember_flag,idx] = ismember(orientation_data_cell.orientation_data{data_idx}{2},tags);
      if(ismember_flag)
         tag_value = idx;
      else
          tags{end+1} = orientation_data_cell.orientation_data{data_idx}{2};
          tag_value = length(tags);
      end

      if(isempty(orientation_data_cell.orientation_data{data_idx}{3}))
          comments = '';
      else
          comments = orientation_data_cell.orientation_data{data_idx}{3};
      end

      az = orientation_data_cell.orientation_data{data_idx}{1}(1);
      ze = orientation_data_cell.orientation_data{data_idx}{1}(2);
      cap = orientation_data_cell.orientation_data{data_idx}{1}(3);
      confidence = orientation_data_cell.orientation_data{data_idx}{1}(4);
  end
  original_az = az;
  original_ze = ze;
  original_cap = cap;

  uicontrol('Style','text','String','House:',...
            'Units','normalized',...
            'HorizontalAlignment','left',...
            'Position',[0.8 0.94 0.05 0.04]);

  edit_go_text_val = uicontrol('Style','edit','String',num2str(data_idx),...
            'Units','normalized',...
            'Position',[0.85 0.95 0.05 0.04],...
            'Callback',{@edit_go});

  uicontrol('Style','pushbutton','String','Go',...
            'Units','normalized',...
            'Position',[0.9 0.95 0.05 0.04],...
            'Callback',{@go_to_house,edit_go_text_val});

  uicontrol('Style','pushbutton','String','Prev View',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.90, 0.075, 0.04], ...
             'Callback', {@prev_view});
  

  uicontrol ('Style', 'pushbutton', 'String', 'Next View', ...
             'Units', 'normalized', ...
             'Position', [0.875, 0.90, 0.075, 0.04], ...
             'Callback', {@next_view});

  uicontrol ('Style', 'text', 'String', 'azimuth', ...
             'Units', 'normalized', ...
             'Position', [0.8, 0.84, 0.15, 0.04]);

  az_edit = ...
  uicontrol ('Style', 'edit', 'String', num2str (az),...
             'Units', 'normalized', ...
             'Position', [0.8, 0.78 + 0.04, 0.15, 0.04], ...
             'Callback', {@edit_az});

  slider_range = 180;
  az_slider = ...
  uicontrol ('Style', 'slider', 'String', 'azimuth', ...
             'Units', 'normalized', ...
             'Position', [0.8, 0.60 + 0.13 + 0.04, 0.15, 0.04], ...
             'min', -180, 'max', 180, 'sliderstep', 2.5/slider_range * [1 1], 'value', az);
  addlistener (az_slider, 'Value', 'PostSet', @slide_az);

  uicontrol ('Style', 'text', 'String', 'zenith',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.55 + 0.13 + 0.04, 0.15, 0.04]);

  ze_edit = ...
  uicontrol ('Style', 'edit', 'String', num2str (ze),...
             'Units', 'normalized', ...
             'Position', [0.8, 0.53 + 0.13 + 0.04, 0.15, 0.04], ...
             'Callback', {@edit_ze});

  slider_range = 90;
  ze_slider = ...
  uicontrol ('Style', 'slider', 'String', 'zenith', ...
             'Units', 'normalized', ...
             'Position', [0.8, 0.48 + 0.13 + 0.04, 0.15, 0.04], ...
             'Min', 0, 'Max', 90, 'sliderstep', 1/slider_range * [1 1], 'value', ze);
  addlistener (ze_slider, 'Value', 'PostSet', @slide_ze);


  uicontrol ('Style', 'text', 'String', 'capacity',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.43 + 0.13 + 0.04, 0.15, 0.04]);

  cap_edit = ...
  uicontrol ('Style', 'edit', 'String', num2str (cap),...
             'Units', 'normalized', ...
             'Position', [0.8, 0.40 + 0.13 + 0.04, 0.15, 0.04], ...
             'Callback', {@edit_cap});

  slider_range = max (3, 2 * cap);
  cap_slider = ...
  uicontrol ('Style', 'slider', 'String', 'capacity', ...
             'Units', 'normalized', ...
             'Position', [0.8, 0.35 + 0.13 + 0.04, 0.15, 0.04], ...
             'Min', 0, 'Max', slider_range, 'sliderstep', 0.1/slider_range * [1 1], 'value', cap);
  addlistener (cap_slider, 'Value', 'PostSet', @slide_cap);

  uicontrol ('Style', 'pushbutton', 'String', 'Undo',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.25 + 0.17 + 0.04, 0.05, 0.04], ...
             'Callback', {@undo});
  uicontrol ('Style', 'pushbutton', 'String', 'Redo',...
             'Units', 'normalized', ...
             'Position', [0.85, 0.25 + 0.17 + 0.04, 0.05, 0.04], ...
             'Callback', {@redo});

  uicontrol ('Style', 'pushbutton', 'String', 'Reset',...
             'Units', 'normalized', ...
             'Position', [0.90, 0.25 + 0.17 + 0.04, 0.05, 0.04], ...
             'Callback', {@reset});

  negative_con_chkbox = uicontrol('Style','checkbox',...
            'Units','normalized',...
            'Enable','off',...
            'Position',[0.8 0.22 0.15 0.04],...
            'String','Show negative est. con.',...
            'Callback',{@show_negative_con_chkbox_cb});

  if(~isnan(data_idx))
      if(current_saved_flag)

          saved_status_string = sprintf('data_idx: %d, Status: Saved\naz: %.2f, ze: %.2f, cap: %.2f\nConfidence: %.1f%%',...
                                        data_idx,...
                                        orientation_data_cell.orientation_data{data_idx}{1}(1),...
                                        orientation_data_cell.orientation_data{data_idx}{1}(2),...
                                        orientation_data_cell.orientation_data{data_idx}{1}(3),...
                                        orientation_data_cell.orientation_data{data_idx}{1}(4));

      else

          saved_status_string = sprintf('data_idx: %d, Status: Unsaved\naz: N/A, ze: N/A, cap: N/A\nConfidence: N/A',...
                                        data_idx);

      end

      uicontrol ('Style', 'edit', 'String', saved_status_string,...
                 'Units', 'normalized', ...
                 'Enable','off',...
                 'Max',10,...
                 'Min',0,...
                 'HorizontalAlignment','left',...
                 'Position', [0.8, 0.14+0.03, 0.15, 0.05]);
  end

%   uicontrol ('Style', 'edit', 'String', '

  %done_button =
  uicontrol ('Style', 'pushbutton', 'String', 'Done',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.08+0.04, 0.075, 0.04], ...
             'Callback', {@done});

  uicontrol ('Style', 'pushbutton', 'String', 'Cancel',...
             'Units', 'normalized', ...
             'Position', [0.875, 0.08+0.04, 0.075, 0.04], ...
             'Callback', {@cancel});


  uicontrol ('Style', 'pushbutton', 'String', 'Next Unsaved Anomalous',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.03+0.04, 0.15, 0.04], ...
             'Callback', {@next_unsaved_anomalous});

  uicontrol ('Style', 'pushbutton', 'String', 'Next Unsaved',...
             'Units', 'normalized', ...
             'Position', [0.8, 0.02, 0.075, 0.04],...
             'Callback', {@next_unsaved});

  uicontrol ('Style', 'pushbutton', 'String', 'Next Anomalous',...
             'Units', 'normalized', ...
             'Position', [0.875, 0.02, 0.075, 0.04],...
             'Callback', {@next_anomalous});

  uicontrol('Style','text',...
            'Units','normalized',...
            'Position',[0.8 0.35 0.025 0.04],...
            'String','Lock:',...
            'HorizontalAlignment','left');

  az_lock_chkbox = uicontrol('Style','checkbox',...
            'Units','normalized',...
            'Position',[0.825 0.36 0.025 0.04],...
            'String','az',...
            'Callback',{@az_lock_chkbox_cb});

  ze_lock_chkbox = uicontrol('Style','checkbox',...
            'Units','normalized',...
            'Position',[0.85 0.36 0.025 0.04],...
            'String','ze',...
            'Callback',{@ze_lock_chkbox_cb});

  uicontrol('Style','pushbutton','String','Re-estimate',...
            'Units','normalized',...
            'Position', [0.875 0.36,0.075,0.04],...
            'Callback',{@reestimate});

  uicontrol('Style','text','String','Tag:',...
            'Units','normalized',...
            'HorizontalAlignment','left',...
            'Position',[0.8 0.265+0.04 0.05 0.04]);

  uicontrol('Style','pushbutton','String','Add comment',...
            'Units','normalized',...
            'Position',[0.80 0.25+0.04 0.05 0.03],...
            'Callback',{@add_comment});

  uicontrol('Style','pushbutton','String','Add new tag',...
            'Units','normalized',...
            'Position',[0.85 0.25+0.04 0.05 0.03],...
            'Callback',{@add_new_tag});

  uicontrol('Style','pushbutton','String','Remove tag',...
            'Units','normalized',...
            'Position',[0.9 0.25+0.04 0.05 0.03],...
            'Callback',{@remove_tag});

  tag_popupmenu = uicontrol('Style','popupmenu','String',tags,...
            'Units','normalized',...
            'Position',[0.85 0.27+0.04 0.10 0.04],...
            'Value',tag_value);

  uicontrol('Style','text','String','Confidence:',...
            'Units','normalized',...
            'Position',[0.8 0.22+0.04 0.075 0.02],...
            'HorizontalAlignment','left');

  confidence_text_val = uicontrol('Style','edit','String',num2str(confidence),...
                                  'Units','normalized',...
                                  'Position',[0.875 0.22+0.04 0.075 0.02],...
                                  'Callback',{@edit_conf_val});

  clr_line_selection = uicontrol('Style','pushbutton','String','Clear Line Selection',...
            'Units','normalized',...
            'Enable','off',...
            'Position', [0.8 0.40,0.15,0.04]);

  draw_traces;
  tic;

  set (f, 'WindowButtonDownFcn', @start_select_or_move);
  set (f, 'WindowButtonUpFcn',   @final_select_or_move);

  % initialize for use below
  downtime = 0;
  coordinates = 0;

  if (~exist ('OCTAVE_VERSION', 'builtin'))
    movegui (f, 'center');
  end
  set (f, 'Visible', 'on');

  uiwait (f);

  function modifiers = get_modifiers (h)
    sel = '';
    % Find actual object that contains selection type
    while isempty (sel)
      try
        sel = get (h, 'selectiontype');
        fprintf ('got selectiontype\n');
      catch
        try
          h = get (h, 'Parent');
          fprintf ('got parent\n');
        catch
          sel = 'oops';
        end
      end
    end

    if (strcmp (sel, 'extend'))
      modifiers = {'shift'};
      %display (modifiers);
    elseif (strcmp (sel, 'alt'))
      modifiers = {'control'};
      %display (modifiers);
    else
      modifiers = {};
    end
  end

  function start_select_or_move (~, ~)
    coordinates = get (ha, 'CurrentPoint');
    coordinates = coordinates (1, 1:2);
    downtime = toc;
  end

  function final_select_or_move (object_handle, ~)
    modifiers = get_modifiers (object_handle);

    coordinates_now = get (ha, 'CurrentPoint');
    coordinates_now = coordinates_now (1, 1:2);
    move = coordinates_now - coordinates;
  end

  function draw_traces
    switch curr_view
      case 'compare'
        set(clr_line_selection,'Enable','off');
        set(negative_con_chkbox,'Enable','off');
        gen_cap  = funcs.angle_coefficient(sunPos, az, ze);
    %     capFactor = max (0, cosd (ze) * sun_pos.s1 ...
    %                       + sind (ze) * cosd (sun_pos.pp - az) .* sun_pos.s2);
        gen = bsxfun (@times, cap/2, gen_cap);    % 48 samples

        solar_range_mod = mod (solar_range, 48);
        sr = min (solar_range_mod):max (solar_range_mod);
%{
        gg = reshape (gen, length (sr), size (gen, 1) / length (sr), size (gen, 2));
        ss = reshape (seen(:, :, solar_range), size (seen, 1), size (seen, 2), length (sr), length (solar_range) / length (sr));
        plot (ha, sr, max (gen(:,  hot_days), [], 2), 'k', ...
                  sr, max (gen(:, cold_days), [], 2), 'k--', ...
                  sr, squeeze (max (seen(:, hot_days,  sr), [], 3))', 'g', ...
                  sr, squeeze (max (seen(:, cold_days, sr), [], 3))', 'r')
        g1 = gg(:, :, 1)';
        g2 = gg(:, :, 2)';
        s1 = squeeze (ss(1, 1, :, :))';
        s2 = squeeze (ss(1, 2, :, :))';
        [~, mn] = min (g1 - s1);
        gg1 = g1(sub2ind (size (g1), mn, 1:size (g1, 2)));
        ss1 = s1(sub2ind (size (s1), mn, 1:size (s1, 2)));
        [~, mn] = min (g2 - s2);
        gg2 = g2(sub2ind (size (g2), mn, 1:size (g2, 2)));
        ss2 = s2(sub2ind (size (s2), mn, 1:size (s2, 2)));
%}
        [~, mn] = min (gen(hot_days, :) - squeeze (seen(1,hot_days,sr)));
        gg1 = gen(sub2ind (size (gen), hot_days(mn), 1:length (sr)));
        ss1 = squeeze (seen(:, sub2ind ([size(seen, 2), size(seen,3)], hot_days(mn), sr)));
        [~, mn] = min (gen(cold_days, :) - squeeze (seen(2,cold_days,sr)));
        gg2 = gen(sub2ind (size (gen), cold_days(mn), 1:length (sr)));
        ss2 = squeeze (seen(:, sub2ind ([size(seen, 2), size(seen,3)], cold_days(mn), sr)));
        plot (ha, sr, gg1, 'k', sr, gg2, 'k--', ...
                  sr, ss1, 'g', sr, ss2, 'r')

        h_lines = findobj(gca,'Type','line');
        xlabel('Time of Day (Half-an-hour)');
        ylabel('Generation (kWh)');
        %legend([h_lines(10);h_lines(9);h_lines(8);h_lines(4)],{'Estimated from model (Summer)','Estimated from model (Winter)','Data from summer days','Data from winter days'});
        colorbar ('off');
      case 'flat'
        set(clr_line_selection,'Enable','off');
        set(negative_con_chkbox,'Enable','off');
        imagesc (ha, data);
        xlabel('Day of Year');
        ylabel('Time of Day (Half-an-hour)');
        colorbar;
        title('Consumption (with Vampire) - Generation');
      case 'corrected_clouds'
        set(negative_con_chkbox,'Enable','on');
        set(clr_line_selection,'Enable','off');
        imagesc (ha, data + get_solar_cor .* irradiation);
        xlabel('Day of Year');
        ylabel('Time of Day (Half-an-hour)');
        colorbar;
        title('Consumption (with Vampire) - Generation + cloud-corrected Estimated Generation');
        negative_marker = find(data + get_solar_cor .* irradiation < 0);
        [negative_marker_i,negative_marker_j] = ind2sub(size(data),negative_marker);
        if(negative_con_chkbox.Value == negative_con_chkbox.Max)
            hold on;
            plot(ha, negative_marker_j,negative_marker_i,'r.');
            hold off;
        end
      case 'corrected'
        set(negative_con_chkbox,'Enable','on');
        set(clr_line_selection,'Enable','off');
        imagesc (ha, data + get_solar_cor);
        xlabel('Day of Year');
        ylabel('Time of Day (Half-an-hour)');
        colorbar;
        title('Consumption (with Vampire) - Generation + Estimated Generation');
        negative_marker = find(data + get_solar_cor < 0);
        [negative_marker_i,negative_marker_j] = ind2sub(size(data),negative_marker);
        if(negative_con_chkbox.Value == negative_con_chkbox.Max)
            hold on;
            plot(ha, negative_marker_j,negative_marker_i,'r.');
            hold off;
        end
      case 'corrected_plot'
        set(negative_con_chkbox,'Enable','off');
        tmp = data + get_solar_cor;
        h_lines = plot (ha, tmp(:, connected));
        set(h_lines(~included_mask(connected)),'LineWidth',2.5);
        set(h_lines(included_mask(connected)),'LineWidth',0.5);
        set(h_lines, 'ButtonDownFcn', {@line_selected, h_lines});
        set(clr_line_selection,'Enable','on');
        set(clr_line_selection,'Callback',{@clear_line_selection,h_lines});
        xlabel('Time of Day (Half-an-hour)');
        ylabel('Consumption (with Vampire) - Generation + Estimated Generation');
        title('All days in a year');
        colorbar ('off');
      case 'corrected_summer'
        set(negative_con_chkbox,'Enable','off');
        set(clr_line_selection,'Enable','off');
        g = get_solar_cor;
        d = data + g;
        tmp = false (size (disconnected));
        tmp(150:250) = true;
        plot (ha, d(:, ~disconnected & tmp));
        xlabel('Time of Day (Half-an-hour)');
        ylabel('Consumption (with Vampire) - Generation + Estimated Generation');
        title('Only days in winter');
      case 'corrected_winter'
        set(negative_con_chkbox,'Enable','off');
        set(clr_line_selection,'Enable','off');
        g = get_solar_cor;
        d = data + g;
        tmp = true (size (disconnected));
        tmp(80:280) = false;
        plot (ha, d(:, ~disconnected & tmp));
        xlabel('Time of Day (Half-an-hour)');
        ylabel('Consumption (with Vampire) - Generation + Estimated Generation');
        title('Only days in summer');
      case 'year_plot'
        set(negative_con_chkbox,'Enable','off');
        tmp = data + get_solar_cor;
        h_lines = plot (ha, tmp(round (10*size(data,1)/24):round(15*size(data,1)/24), connected)');
        xlabel('Day of year');
        ylabel('Consumption (with Vampire) - Generation + Estimated Generation');
        title('plot vs year');
        colorbar ('off');
    end
  end

  function solar_correction = get_solar_cor
    solar_correction = zeros (size (data));
    for i = 1:size (data, 2)
      for j = 0:(s.dark_start-s.dark_end)
        % Estimate generation at this time, assuming no cloud
        p1 = cosd(ze);
        p2 = sind(ze).*cosd(s.pp(i,j+1) - az);
        solar_correction(j+s.dark_end, i) ...
                  = max (0,   bsxfun (@times, s.s1(i,j+1), p1) ...
                            + bsxfun (@times, s.s2(i,j+1), p2))';
      end
    end
    solar_correction = solar_correction * (24 / size(data,1) * cap);
    solar_correction(:, disconnected) = 0;
  end

    function edit_conf_val(~,~)
        num_representation = str2num(get(confidence_text_val,'String'));
        if(isempty(num_representation) || (num_representation < 0) || (num_representation > 100))
            set(confidence_text_val,'String',num2str(default_confidence));
            confidence = default_confidence;
        else
            confidence = num_representation;
        end
    end

    function edit_go(object, event, h)
        num_representation = str2num(get(object,'String'));
        if(isempty(num_representation) || ~ismember(num_representation,meta.solar_users))
            set(object,'String',num2str(data_idx));
            return;
        end
    end

    function go_to_house(object, event, h)
        jump_to_house = str2num(h.String);
        delete(f);
    end

    function add_comment(object, event, h)
        comment_box = figure('Position',[0 0 300 200], 'MenuBar', 'none', ...
                            'ToolBar', 'none', 'Visible', 'off', ...
                            'Name', 'Add comment', 'Resize', 'off', ...
                            'WindowStyle','modal');

        movegui(comment_box,'center');

        comment_text = uicontrol('Parent',comment_box,'Style','edit','String', comments,...
                         'Units', 'normalized', ...
                         'Enable','on',...
                         'Max',10,...
                         'Min',0,...
                         'HorizontalAlignment','left',...
                         'Position', [0.1, 0.20, 0.8, 0.8]);


        uicontrol('Parent',comment_box,'Style','pushbutton','String','OK',...
                  'Units','normalized',...
                  'Position',[0.1, 0.05, 0.3, 0.15],...
                  'Callback',{@okay_comment_box});

        uicontrol('Parent',comment_box,'Style','pushbutton','String','Cancel',...
                  'Units','normalized',...
                  'Position',[0.6, 0.05, 0.3, 0.15],...
                  'Callback',{@cancel_comment_box});

        set(comment_box,'Visible','on');

        uiwait(comment_box);


        function cancel_comment_box(~,~)
            delete(comment_box);
        end

        function okay_comment_box(~,~)
            comments = comment_text.String;
            delete(comment_box);
        end
    end

    function add_new_tag(~,~)
        new_tag_name = inputdlg('Tag name','Add new tag');

        if(isempty(new_tag_name) || isempty(new_tag_name{1}))
            return;
        else
            % Append new tag to existing tags
            tags{end+1} = new_tag_name{1};
            tag_popupmenu.String = tags;
            save(tag_group_file,'tags');
        end
    end

    function remove_tag(~,~)
        % Remove currently selected
        if(tag_popupmenu.Value > 1)
            tags(tag_popupmenu.Value) = [];
            tag_popupmenu.Value = tag_popupmenu.Value-1;
            tag_popupmenu.String = tags;
            save(tag_group_file,'tags');
        end
    end

  function next_view (~, ~)
    switch curr_view
      case 'compare'
        curr_view = 'flat';
      case 'flat'
        curr_view = 'corrected';
      case 'corrected'
        curr_view = 'corrected_clouds';
      case 'corrected_clouds'
        curr_view = 'corrected_plot';
      case 'corrected_plot'
        curr_view = 'corrected_winter';
      case 'corrected_winter'
        curr_view = 'corrected_summer';
      case 'corrected_summer'
        curr_view = 'year_plot';
      case 'year_plot'
        curr_view = 'compare';
    end
    draw_traces;
  end

    function prev_view(~,~)
        switch curr_view
          case 'year_plot'
              curr_view = 'corrected_summer';
          case 'corrected_summer'
              curr_view = 'corrected_winter';
           case 'corrected_winter'
             curr_view = 'corrected_plot';
          case 'corrected_plot'
            curr_view = 'corrected_clouds';
          case 'corrected_clouds'
             curr_view = 'corrected';
           case 'corrected'
             curr_view = 'flat';
          case 'flat'
            curr_view = 'compare';
          case 'compare'
            curr_view = 'year_plot';
        end
        draw_traces;
    end

  %% Sliders
  function edit_az (~, ~)
    set_az (str2double (get (az_edit, 'String')));
  end

  function slide_az (~, ~)
    set_az (get (az_slider, 'value'));
  end

  function set_az (new_val)
    old_val = str2num(get (az_edit,'String'));
    add_undo ('az', old_val);
    az = new_val;
    %display (az);
    set (az_edit, 'String', num2str (new_val));
    set (az_slider, 'value', new_val);
    draw_traces;
    find_mismatch
  end

  function edit_ze (~, ~)
    set_ze (str2double (get (ze_edit, 'String')));
  end

  function slide_ze (~, ~)
    set_ze (get (ze_slider, 'value'));
  end

  function set_ze (new_val)
    old_val = str2num(get (ze_edit,'String'));
    add_undo ('ze', old_val);
    ze = new_val;
    %display (ze);
    set (ze_edit, 'String', num2str (new_val));
    set (ze_slider, 'value', new_val);
    draw_traces;
    find_mismatch
  end

  function edit_cap (~, ~)
    set_cap (str2double (get (cap_edit, 'String')));
  end

  function slide_cap (~, ~)
    set_cap (get (cap_slider, 'value'));
  end

  function set_cap (new_val)
    old_val = str2num(get (cap_edit,'String'));
    add_undo ('cap', old_val);
    cap = new_val;
    %display (cap);
    set (cap_edit, 'String', num2str (new_val));
    set (cap_slider, 'value', new_val);
    draw_traces;
    find_mismatch
  end

    %% Parameter lock checkbox callbacks
    function az_lock_chkbox_cb(~,~)

        % Only az or ze can be locked. Didn't use a radio button for this
        % since radio button must have one selection, whereas in reality,
        % az and ze can be both unlocked.
        if(ze_lock_chkbox.Value == ze_lock_chkbox.Max)
            ze_lock_chkbox.Value = ze_lock_chkbox.Min;
        end

        ze_locked_flag = false;
        az_locked_flag = logical(az_lock_chkbox.Value);
    end

    function ze_lock_chkbox_cb(~,~)
        if(az_lock_chkbox.Value == az_lock_chkbox.Max)
            az_lock_chkbox.Value = az_lock_chkbox.Min;
        end
        az_locked_flag = false;
        ze_locked_flag = logical(ze_lock_chkbox.Value);
    end

    %% Line selection for manual outlier removal.
    function line_selected(object, event, h)
        c = connected(object == h);
        included_mask(c) = ~included_mask(c);

        if(~included_mask(c))
            fprintf(1,'Day (line) selected: %d\n',c);
        end

        set(h(~included_mask(connected)),'LineWidth',2.5);
        set(h(included_mask(connected)),'LineWidth',0.5);
    end

    function clear_line_selection(object, event, h)
        included_mask = true(365,1);
        set(h(included_mask(connected)),'LineWidth',0.5);
    end

  %% Chekbox callback for showing red spots on the imagesc plot of the corrected data
  function show_negative_con_chkbox_cb(obj,evt,h)
    draw_traces;
  end

  function find_mismatch
    local_data = data';
    local_data(~included_mask,:) = inf;

    % New "seen".  with or without vampires, summer or winter, half hour      clear seen
    % Assign "most summery" half of non-missing days to "one",
    % and "most wintry" half to "two".
    g = find(~disconnected);
    max_len = 5;
    sp = s.SunPos(g(max_len:max_len:end), :);
    local_seen = seen(:, :, solar_range);
    local_sunPos = sunPos;

    % Recompute "seen" with outlying days removed.

    [mx, mx_pos] = max(squeeze (local_seen(1,:,:)), [], 2);
    mx(hot_days)  = max (mx(hot_days));
    mx(cold_days) = max (mx(cold_days));

    % Only consider data points at least 25% of peak
    % to exclude both times of high demand and ambient light near dawn/dusk
    big = bsxfun(@gt, squeeze (local_seen(1,:,:)), mx/4);
    % If there was a mismatch of azimuths,
    % only match to peak of less reliable one.

    % Find times where the generation is big enough to be trustworthy.
    % This avoids fitting morning/evening (subject to shadows and
    % reflections), and also times with regular daily use, such as the
    % "morning routine" or a pool pump.
    % This is used for finding a lower-bound on capacity,
    % and so it more conservative than the 25% used for the matching.
    big_all = bsxfun (@gt, -local_data(:, s.dark_end:s.dark_start), -0.75 * min(local_data, [], 2));
    sun_pos.pp = s.pp(big_all);
    sun_pos.s1 = s.s1(big_all);
    sun_pos.s2 = s.s2(big_all);

    % solar_mismatch() will ensure that  cap  is big enough to account
    % for all of the observations recorded in few_data.
    % We could use all tiemsteps, but for efficiency,
    % we only consider the most negative readings over
    % contiguous 5-day intervals.
    my_data = local_data(:, s.dark_end:s.dark_start);

    max_seen = double (max (local_seen(:)));
    cost = funcs.solar_mismatch (double ([az, ze]), local_sunPos, ...
                                 double (local_seen), big, max_seen, my_data(big_all), sun_pos);
    fprintf ('mismatch %g\n', cost);
  end

  %% Reestimation with the selected outlying days removed.
  % TODO: Refactor to reduce duplication
  function reestimate(~,~)
      % TODO: Reestimate orientation params using only days where bitmasks
      % are true.

      % Lines which are not included in the optimisation procedure are
      % assigned to infinity so that the minimisation in
      % compute_peak_sunpos() and compute_seen().
      local_data = data';
      local_data(~included_mask,:) = inf;
      local_data_no_vamp = data_no_vamp';
      local_data_no_vamp(~included_mask,:) = inf;

%       [spos, wpos, bpos] = compute_peak_sunpos(reshape(local_data,[1 size(local_data,1) size(local_data,2)]), reshape(local_data_no_vamp,[1 size(local_data_no_vamp,1) size(local_data_no_vamp,2)]), meta, s);
%
%       % ASK: sps, spw selection based on "disconnected"?
%       wpos = bpos;
%
%       sps = spos;
%       spw = wpos;
%       local_sunPos = [sps;spw];
%
%       [seen_summer, seen_winter] = compute_seen(reshape(local_data,[1 size(local_data,1) size(local_data,2)]), reshape(local_data_no_vamp,[1 size(local_data_no_vamp,1) size(local_data_no_vamp,2)]), smr, rr, meta);
%       seen = [seen_summer, seen_winter];
%       % seen = compute_seen(data_no_vamp,
      % New "seen".  with or without vampires, summer or winter, half hour      clear seen
      % Assign "most summery" half of non-missing days to "one",
      % and "most wintry" half to "two".
      g = find(~disconnected);
      max_len = 5;
      sp = s.SunPos(g(max_len:max_len:end), :);
      local_seen = seen(:, :, solar_range);
      local_sunPos = sunPos;

      % Recompute "seen" with outlying days removed.

      locked_az = az;
      locked_ze = ze;

      [mx, mx_pos] = max(squeeze (local_seen(1,:,:)), [], 2);
      mx(hot_days)  = max (mx(hot_days));
      mx(cold_days) = max (mx(cold_days));

      % Only consider data points at least 25% of peak
      % to exclude both times of high demand and ambient light near dawn/dusk
      big = bsxfun(@gt, squeeze (local_seen(1,:,:)), mx/4);

      % Find times where the generation is big enough to be trustworthy.
      % This avoids fitting morning/evening (subject to shadows and
      % reflections), and also times with regular daily use, such as the
      % "morning routine" or a pool pump.
      % This is used for finding a lower-bound on capacity,
      % and so it more conservative than the 25% used for the matching.
      big_all = bsxfun (@gt, -local_data(:, s.dark_end:s.dark_start), -0.75 * min(local_data, [], 2));
      sun_pos.pp = s.pp(big_all);
      sun_pos.s1 = s.s1(big_all);
      sun_pos.s2 = s.s2(big_all);

      % solar_mismatch() will ensure that  cap  is big enough to account
      % for all of the observations recorded in few_data.
      % We could use all tiemsteps, but for efficiency,
      % we only consider the most negative readings over
      % contiguous 5-day intervals.
      len = 1:5 * floor (size (local_data, 1) / 5);
      my_data = local_data(:, s.dark_end:s.dark_start);
      d = local_data(len, s.dark_end:s.dark_start);
      few_big = bsxfun (@rdivide, -d, min(local_data(len, :), [], 2));
      few_big = reshape (few_big, len(end)/5, 5 * (s.dark_start - s.dark_end + 1));
      d = reshape  (d,  len(end)/5, 5 * (s.dark_start - s.dark_end + 1));
      % pp = reshape (s.pp, len(end)/5, 5 * (s.dark_start - s.dark_end + 1));
      % s1 = reshape (s.s1, len(end)/5, 5 * (s.dark_start - s.dark_end + 1));
      % s2 = reshape (s.s2, len(end)/5, 5 * (s.dark_start - s.dark_end + 1));
      [~, few_idx] = max (few_big, [], 2);
      idx = sub2ind (size (few_big), 1:size (few_big, 1), few_idx');
      sun_pos_few.pp = s.pp(idx);
      sun_pos_few.s1 = s.s1(idx);
      sun_pos_few.s2 = s.s2(idx);
      few_data = d(idx);


      if length(local_data) > 1
          % Initial estimates, chosen heuristically:
          % Capacity is the maximum observed nett generation
          % Azimuth is roughly proportional to the square of
          % the time of the peak nett generation (with midday being 0),
          % clipped to be between due east and due west
          cap = max(mx);
          az = (((mean(mx_pos) + solar_range(1)-1) * -360/meta.SamPerDay) + 180);
          az = az * abs(az/10);
          az = max (-90, min (90, az));
          ze = 10;		% ~ sun's summer zenith angle in Melbourne
      else                % initialise from previous run
          cap = local_data;
          az = smr;
          ze = rr;
      end
      options = optimoptions('fmincon', 'Display', 'off', 'Algorithm','active-set');
      % DEBUG CODE
%       sseen = zeros (size (data, 1), (s.dark_start - s.dark_end)+1);
%       capFactor = sseen;
%       for i = 1:size(data, 2)
%         for j = 0:(s.dark_start-s.dark_end)
%           % Estimate generation at this time, assuming no cloud
%           p1 = cosd(ze);
%           p2 = sind(ze).*cosd(s.pp(i,j+1) - az);
%           capFactor(i,j+1) ...
%             = max (0,   bsxfun (@times, s.s1(i,j+1), p1) ...
%             + bsxfun (@times, s.s2(i,j+1), p2))';
%         end
%       end
% %       dat = squeeze (evalin ('caller', 'data(i,:,s.dark_end:s.dark_start)'));
%       dat = data(s.dark_end:s.dark_start,:)';
%       imagesc (figure1, (dat + capFactor * cap/2)');    % 48 samples
%       plot (figure2, (dat + capFactor * cap/2)');       % 48 samples
%       plot (figure3, dat');
%       % END DEBUG CODE

      max_seen = double (max (local_seen(:)));
      [cost, dy] = funcs.solar_mismatch (double ([az, ze]), local_sunPos, ...
                                   double (local_seen), big, max_seen, my_data(big_all), sun_pos);
      if ~isfinite (cost)
        [az, ze, big] = find_feasible (az, ze, local_sunPos, ...
                                             double (local_seen), big, ...
                                             max_seen, my_data(big_all), sun_pos, dy, solar_range);
      end

      if(az_locked_flag)
          initial_params = ze;
          lower_bound = 1;
          upper_bound = 45;

          f2 = @(X) funcs.solar_mismatch ([locked_ze, X], local_sunPos, double (local_seen),...
                                                   big, max_seen, my_data(big_all), sun_pos);
      elseif(ze_locked_flag)
          f2 = @(X) funcs.solar_mismatch ([X, locked_ze], local_sunPos, double (local_seen),...
                                                   big, max_seen, my_data(big_all), sun_pos);

          initial_params = az;
          lower_bound = -90;
          upper_bound = 90;
      else
          f2 = @(X) funcs.solar_mismatch (X, local_sunPos, double (local_seen),...
                                                   big, max_seen, my_data(big_all), sun_pos);

          initial_params = [az ze];
          lower_bound = [-90 1];
          upper_bound = [90 45];
      end

      old_cost = 1e90;
      for n = 1:2
        tic
        % Lower limit of 1 for ze, as ze=0 gives plateau for az.
        [X, cost] = fmincon (f2, ...
                             initial_params, ...
                             [], [], [], [], lower_bound, upper_bound, [], ...
                             options);

        if(az_locked_flag)
            az = locked_az;
            ze = X;
        elseif(ze_locked_flag)
            az = X;
            ze = locked_ze;
        else
            az  = X(1);
            ze  = X(2);
        end
        [~, ~, cap] = funcs.solar_mismatch ([az, ze], local_sunPos, double (local_seen), big, ...
                              max_seen, my_data(big_all), sun_pos);

        % fprintf ('pass 1: %g\n', toc);
        if ~isfinite (cap)
          tic
          [X, cost] = fmincon (f2, ...
                               initial_params, ...
                               [], [], [], [], lower_bound, upper_bound, [], ...
                               options);
            if(az_locked_flag)
                az = locked_az;
                ze = X;
            elseif(ze_locked_flag)
                az = X;
                ze = locked_ze;
            else
                az  = X(1);
                ze  = X(2);
            end
          [~, ~, cap] = funcs.solar_mismatch ([az, ze], local_sunPos, double (local_seen), big, ...
                                max_seen, my_data(big_all), sun_pos);
          %fprintf ('\t\t\tpass 2: %g\n', toc);
        end
        if cost > 0.99 * old_cost && cost <= old_cost
          break;
        end
        old_cost = cost;
      end
      disp ([az, ze, cap, cost])
      gen   = cap*funcs.angle_coefficient(local_sunPos, az, ze);
%       plot (figure11, capFactor');

      capacity = cap;
      mx_start = 5;
      dy = squeeze (seen(1, :, solar_range)) - gen;
      mismatch = max(max(abs(dy(:,mx_start:end-mx_start))));
      mismatch = mismatch / max(gen(:));
      draw_traces;

      action{1} = 'az';
      action{2} = az;
      apply_action(action);

      action{1} = 'ze';
      action{2} = ze;
      apply_action(action);

      action{1} = 'cap';
      action{2} = cap;
      apply_action(action);

  end

  %% Undo
  function add_undo (var, val)
    undo_stack{end+1} = {var, val};
    redo_stack = {}; % TODO: redo functionality is still broken.
  end

  function undo (~, ~)
    if ~isempty (undo_stack)
      action = undo_stack{end};
      undo_stack = undo_stack(1:end-1);
      redo_stack{end+1} = action;
      apply_action (action)
    end
  end

  function redo (~, ~)
    if ~isempty (redo_stack)
      action = redo_stack{end};
      redo_stack = redo_stack(1:end-1);
      undo_stack{end+1} = action;
      apply_action (action);
    end
  end

  function reset(~,~)

      % Reset to values obtained from auto-fitting.

      action{1} = 'az';
      action{2} = original_az;
      apply_action(action);

      action{1} = 'ze';
      action{2} = original_ze;
      apply_action(action);

      action{1} = 'cap';
      action{2} = original_cap;
      apply_action(action);

      redo_stack = {};
      undo_stack = {};

  end

  function apply_action (action)
    if isequal (action{1}, 'az')
      az = action{2};
      set (az_edit, 'String', num2str (az));
      set (az_slider, 'value', az);
      draw_traces;
    elseif isequal (action{1}, 'ze')
      ze = action{2};
      set (ze_edit, 'String', num2str (ze));
      set (ze_slider, 'value', ze);
      draw_traces;
    else
      cap = action{2};
      set (cap_edit, 'String', num2str (cap));
      set (cap_slider, 'value', cap);
      draw_traces;
    end
  end

  function done (~, ~)
    % TODO: Save results to file
    abort = false;
    if(~isnan(data_idx))
        orientation_data = orientation_data_cell.orientation_data;
        orientation_data{data_idx}{1} = [az, ze, cap, confidence];
        orientation_data{data_idx}{2} = tag_popupmenu.String{tag_popupmenu.Value};
        orientation_data{data_idx}{3} = comments;
    end
    delete (f);
  end

  % TODO: function to skip to the next house

  function cancel (~, ~)
    abort = true;
    % assign nan value to saved data
    if(~isnan(data_idx))
%         orientation_data = orientation_data_cell.orientation_data;
%         orientation_data(data_idx,:) = [az, ze, cap];% nan(1,3);
%         if(exist(save_file,'file'))
%
%             load(save_file);
%             orientation_data(data_idx,:) = nan(1,3);
%
%         else
%
%             orientation_data(1,:) = nan(1,3);
%
%             save(save_file,'orientation_data');
%
%         end
%         save(save_file,'orientation_data');
    end
%     done;
    delete(f);
  end

    function next_unsaved(~, ~)

        % if current_saved_flag is true, skip until we reach a data entry
        % we haven't saved yet.
        if(current_saved_flag)

            jump_to_unsaved = true;
            delete(f);

        end

    end

    function next_anomalous(~,~)
        jump_to_house = find_anomalous_id(data_idx);
        if(~isnan(jump_to_house))
            delete(f);
        end
    end

    function next_unsaved_anomalous(~,~)
        if(current_saved_flag)
            jump_to_house = find_anomalous_id(data_idx);
            if(~isnan(jump_to_house))
                jump_to_unsaved = true;
                delete(f);
            end
        end
    end

end

