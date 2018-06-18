function generate_aircon_training (cor_ac, data, gex, redo_old, compare, meta, todo)
  %redo_old = false;   % Re-edit cases we have already checked?
  %compare = true;
  if nargin < 7
    todo = true (1, size (data, 1));
  end

  % Errors:

  if exist ([meta.metaDataPath, '../data/aircon_labelled.mat'], 'file')
    load ([meta.metaDataPath, '../data/aircon_labelled'])
  else
    for r = 1:size (data, 1)
      aircon_labelled(r).done = false;
      aircon_labelled(r).cor = sparse (size (data, 3), size (data, 2));
    end
  end
  if ~isfield (aircon_labelled, 'manual_confidence')
    init = sparse (numel (aircon_labelled(1).cor), 1);
    init(1) = 1;
    for r = 1:length (aircon_labelled)
      aircon_labelled(r).manual_confidence = init;
    end
  end

  step = 1;   % go forwards.  ("prev" sets step to -1.)
  i = 1;

  % Allow "prev" to work with redo_old==false
  edited_this_session = false (size (aircon_labelled));

  while i <= length (aircon_labelled)
    if todo(i) && (redo_old || aircon_labelled(i).done ~= true ...
                            || edited_this_session(i))
      edited_this_session(i) = true;
      aircon = cor_ac(:, :, i);

      if ~compare && aircon_labelled(i).done ~= false        % Re-display edited version, if not comparing new alg
        aircon = reshape (aircon_labelled(i).cor, size (aircon));
      end

      cv = squeeze (data (i, :, :))';

      % Call gui
      confidence = reshape (cumsum (aircon_labelled(i).manual_confidence), ...
                            size (cv));
      [aircon, ok, confidence] = aircon_ground_truth (cv, aircon, i, gex, confidence, meta);
      step = 1;         % (leave "prev" search)
      if ok == -2       % Quit
        break;
      elseif ok == -1   % Prev
        step = -1;
        i = i + step;
      elseif ok         % Next
        aircon_labelled(i).cor = aircon;
        aircon_labelled(i).manual_confidence = [confidence(1)
                                                diff(confidence(:))];
        if ok == 1
          aircon_labelled(i).done = true;
        else
          aircon_labelled(i).done = -1;
        end
        save ../data/aircon_labelled aircon_labelled
        i = i + step;
      else              % Skip
        i = i + step;
      end
    else
       i = i + step;
    end
    if i < 1
      i = 1;
      step = 1;
    end
  end

end

function scratch
  % code to copy-and-paste into the command window

  figure(1); hold off;
  curves = zeros (length (meta.summer), size (cor_ac, 3));
  for i = 1:size (cor_ac, 3)
    c = cor_ac(:, meta.summer, i);
    s = sum (c ~= 0);
    curves(:, i) = sort(s)';
    %plot (sort (s), (1:length (s)) / length (s));
    %hold on;
    %keyboard
  end

  quantiles = [0.6 0.7, 0.8, 0.9];
  suspect = cell (size (quantiles));
  for i = 1:length(quantiles)
    [c, idx] = sort (curves(ceil (quantiles(i) * end), :));
    first = find (c, 1);
    last = first + round (0.8 * (length (c) - first));
    thresh = c(last) / 0.8;
    suspect{i} = idx (c >= thresh);
  end

  samPerDay = 48;

  % Reshape to an array indexed by (hrs, user)
  [d, u] = find (diff (curves) ~= 0);
  vals = curves (sub2ind (size (curves), d+1, u));
  times = zeros (samPerDay, size (curves, 2));
  deleted = [];
  for i = 1:samPerDay
    day = (vals == i);
    times (i, u(day)) = d(day);
    [smallest, idx] = min (times(i,:));
    next = min (times(i, 1: size (times, 2) ~= idx));
    if next - smallest > 10
      deleted = [deleted, idx];
      idx = (u == idx);
      d = d(~idx);
      u = u(~idx);
    else
      keyboard
    end
  end


  sus = suspect{1};
  for i = 2:length (suspect)
    sus = intersect (sus, suspect{i});
  end
  suspect = sus;

  for i = suspect(30:end)
    figure(3); imagesc (cor_ac(:, :, i));
    figure(4); imagesc (squeeze (corrected(i, :, :))');
    figure(5); imagesc (squeeze (corrected(i, :, :))' - cor_ac(:, :, i));
    keyboard
  end

  % Too much:
  % i == 30, 79, 111
end
