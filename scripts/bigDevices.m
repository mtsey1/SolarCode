function [cor, conf, thresh, jump_times, bins] = bigDevices (user_data, vamp, days_train, hrs_train, days_use, hrs_use)
  if nargin < 3
    hrs_train  = size (user_data, 1);
    days_train = size (user_data, 2);
  end
  if nargin < 5
    hrs_use  = hrs_train;
    days_use = days_train;
  end
  
  cv = user_data(hrs_use, days_use);
  train_data = user_data(hrs_train, days_train);
  train_data = train_data(:);
  if all (isnan (train_data))
    cor = zeros (size (user_data));
    conf = 0;
    thresh = 0;
    jump_times = [];
    bins = [];
    return;
  end
  %n = 2;		% time signal must be ~flat each side of a jump

  %    diffs = diff(train_data(:));
  %    diffs(~isfinite(diffs)) = 0;
  %    ab = cumsum(abs(diffs));
  %    % declare a "jump" if change this time step is as big as the sum of
  %    % the changes in the  n  previous and  n  subsequent time steps.
  %    diffs = abs(diffs);
  %    jump_times = find(diffs((n+1):(end-n)) > 0.5*(ab((2*n+1):end)-ab(1:end-(2*n)))) + n;
  %    jump_sizes = diffs(jump_times);
  sud = sort (train_data);
  [alt(1, 1:2), rectj, ij, jj, aj, bj, pj, ipj, binsj] = common_jumps (train_data, sud);

  % peak detection in histogram of absolute power levels
  off_peak = [];
  on_peak = [];
  lbins = 30; ml = max (train_data); mml = min (train_data);
  %while any (cl == 0)
  lbins = lbins+1;
  binsl = (0.5:1:lbins) * (ml-mml)/lbins - mml;
  [al, bl] = hist (train_data,binsl);
  %figure(5); hist(train_data, binsl);
  % break ties by minor smoothing
  [pl, cl, sepl, rectl] = secondPeak (al + 0.1*([al(2:end),0] + [0,al(1:end-1)]),5);
  %end
  pl(cl==0) = NaN;	%Ignore non-peaks
  % Find "big" peak that isn't from small switches etc.
  if abs (pl(2) - pl(1)) < max (pl(isfinite (pl)))/3
    big = pl(3);
    on_peak = 3;
    off_peak = 2;
  else
    big = pl(2);
    on_peak = 2;
    off_peak = 1;
  end
  % Look for two peaks that roughly add to a third: a sign of >= 2 A/Cs
  % This time, we need to consider "noise" introduced by small devices
  [pp, ipl] = sort (pl);
  if pp(2) < pp(3)/2 && cl(ipl(2)) > cl(ipl(1))
    ppp = [pp+pp(1), pl+pp(2)];	% + not -, as compared with pp(.)+pp(.)
  else
    ppp = pp+pl(1);	% + not -, as compared with pp(.)+pp(.)
  end
  ppp = [ppp-1, ppp, ppp+1];
  ppp = ppp(:);
  if ismember (pp(3) + pp(4), ppp) && min (abs ([pp(3), pp(4), pp(3)+pp(4)-pp(1)]-big)) < 2
    if abs(pp(3) - big) < 2
      alt(2,1) = bl(pp(3));
      alt(2,2) = bl(pp(4));
    else
      alt(2,1) = bl(pp(4));
      alt(2,2) = bl(pp(3));
    end
    alt(2,:) = alt(2,:) - min (bl(pp(isfinite(pp))));
    il=3; jl=4;
    on_peak = 3;
    %fprintf('Two ACs (level): %d & %d:  bl(%d)=%f  bl(%d)=%f\n', ii, jj, pp(ii), bl(pp(ii)), pp(jj), bl(pp(jj)));
  elseif ismember (pp(2) + pp(3), ppp) && min (abs ([pp(3), pp(2), pp(3)+pp(2)-pp(1)]-big)) < 2
    if abs(pp(2) - big) < 2
      alt(2,1) = bl(pp(2));
      alt(2,2) = bl(pp(3));
    else
      alt(2,1) = bl(pp(3));
      alt(2,2) = bl(pp(2));
    end
    alt(2,:) = alt(2,:) - min (bl(pp(isfinite(pp))));
    il=2; jl=3;
    on_peak = 2;
    %fprintf('Two ACs (level): %d & %d:  bl(%d)=%f  bl(%d)=%f\n', ii, jj, pp(ii), bl(pp(ii)), pp(jj), bl(pp(jj)));
  elseif ismember (pp(2) + pp(4), ppp) && min (abs ([pp(4), pp(2), pp(4)+pp(2)-pp(1)]-big)) < 2
    if abs(pp(2) - big) < 2
      alt(2,1) = bl(pp(2));
      alt(2,2) = bl(pp(4));
    else
      alt(2,1) = bl(pp(4));
      alt(2,2) = bl(pp(2));
    end
    alt(2,:) = alt(2,:) - min (bl(pp(isfinite(pp))));
    il=2; jl=4;
    on_peak = 2;
    %fprintf('Two ACs (level): %d & %d:  bl(%d)=%f  bl(%d)=%f\n', ii, jj, pp(ii), bl(pp(ii)), pp(jj), bl(pp(jj)));

  else
    % Just look for the largest peak of a "large" device
    % If the second largest peak is also small, skip it
    thresh = train_data(ceil (length (train_data)*0.97))/3;
    
    pl = pl(isfinite(pl));
    on_peak = min (on_peak, length (pl));
    
    while length (pl) > 2 && (abs (bl(pl(1))-bl(pl(2))) < thresh)
      pl = pl([1 3:length(pl)]);
    end
    if length (pl)>1
      alt(2,1) = bl(max (pl(1), pl(2))) - bl(min (pl(1),pl(2)));
      on_peak = 2;
      off_peak = 1;
      %fprintf('%f-%f = %f\n',  bl(max(pl(1), pl(2))),  bl(min(pl(1),pl(2))), alt(2,1) );
    else
      alt(2,1) = 0;
    end
    alt(2,2) = 0;
  end

  bins = {binsj; binsl};

  dev = [mean([alt(1,1), alt(2,1)]), mean([alt(2,1), alt(2,2)])];
%{
  fprintf ('sepl ');  fprintf (' %g', sepl(:,3));  fprintf ('\n');
  fprintf ('rectl '); fprintf (' %g', rectl(:,3)); fprintf ('\n');
  fprintf ('pl '); fprintf (' %g', binsl(pl(isfinite (pl)))); fprintf ('\n');
  if ~isempty (jump_times)
    fprintf ('sepj ');  fprintf (' %g', sepj(:,3));  fprintf ('\n');
    fprintf ('rectj '); fprintf (' %g', rectj(:,3)); fprintf ('\n');
    fprintf ('pj '); fprintf (' %g', binsj(pj(isfinite (pj)))); fprintf ('\n');
  end
%}
  if ~isempty (ij)
    pl = pl(isfinite (pl));
    pj = pj(isfinite (pj));
    % look for match between a pair of levels and a jump
    lev_diff = bsxfun (@minus, bl(pl), bl(pl)');
    lev_diff = triu (lev_diff) + tril (nan (size (lev_diff)));
    pj = pj(isfinite (pj));
    match = abs (bsxfun (@minus, bj(pj), lev_diff(:)));

    wgts = bsxfun (@times, sqrt (al(pl)), sqrt (al(pl))');
    wgts = bsxfun (@times, sqrt ([1, aj(pj(2:end))]), wgts(:));
    match = match ./ wgts;

    [val, row] = min (match);
    [~, col] = min (val);
    row = row(col);
    larger = ceil (row / length (pl));
    smaller = row - (larger-1) * length (pl);
    [gj, gap_j] = max (rectj(:,3));
    [gl, gap_l] = max (rectl(:,3));
    [ggj, next_j] = max (rectj(rectj(:,3) < gj, 3));
    [ggl, next_l] = max (rectl(rectl(:,3) < gl, 3));
    if gap_j > size (rectj, 1)/2
      big_j = gap_j - size (rectj, 1)/2;
    else
      big_j = ipj(gap_j);
    end
    if gap_l > size (rectl, 1)/2
      big_l   = gap_l - size (rectl, 1)/2;
      small_l = big_l - 1;
    elseif gap_l > 1
      big_l   = ipl(gap_l);
      small_l = ipl(gap_l - 1);
    else
      big_l = ipl(2);
      small_l = ipl(1);
    end
    ggj = gj / ggj;
    ggl = gl / ggl;
%{
    fprintf ('Is %g = %g - %g?  (%g%%) big_j %d %d big_l %d %d small_l %d %d gj/ggj %g gl/ggl %g\n', ...
      bj(pj(col)), bl(pl(larger)), bl(pl(smaller)), ...
      abs (bj(pj(col)) - (bl(pl(larger)) - bl(pl(smaller)))) / bj(pj(col)) * 100, ...
      big_j, col, ...
      big_l, larger, ...
      small_l, smaller, ...
      ggj, ggl);
%}
%    disp(alt);
    if abs (bj(pj(col)) - (bl(pl(larger)) - bl(pl(smaller)))) / bj(pj(col)) < 0.1
      if abs (alt(1,1) - alt(2,1)) / alt(1,1) < 0.1 ...
            && abs (alt(1,1) - bj(pj(col))) < 0.1
        pow = min (alt(1,1), alt(2,1));
      elseif smaller == 1
        two_jumps = (alt(1,2) ~= 0);
        two_level = (alt(2,2) ~= 0);
        if two_jumps && ~two_level
          pp = sort (pj);
          if aj(pp(ij)) + aj(pp(jj)) < aj(pj(col))
            alt = [bj(pj(col)), 0; (bl(pl(larger))-bl(pl(smaller))), 0];
            on_peak = larger;
            off_peak = smaller;
          end
        elseif two_level && ~two_jumps
          pp = sort (pl);
          if al(pp(il)) + al(pp(jl)) < (al(pl(larger)) + al(pl(smaller)))
            alt = [bj(pj(col)), 0; (bl(pl(larger))-bl(pl(smaller))), 0];
            on_peak = larger;
            off_peak = smaller;
          end
        elseif two_jumps && two_level
          alt = -sort (-alt, 2);
          if     alt(1,1) < 0.9 * alt(2,1) || alt(2,1) < 0.9 * alt(1,1) ...
              || alt(1,2) < 0.9 * alt(2,2) || alt(2,2) < 0.9 * alt(1,2)
            alt = [bj(pj(col)), 0; (bl(pl(larger))-bl(pl(smaller))), 0];
            on_peak = larger;
            off_peak = smaller;
          end
        else
          if bj(pj(col)) == alt(1,1) ...
             || bl(pl(larger)) - bl(pl(smaller)) == alt(2,1)
            alt = [bj(pj(col)), 0; (bl(pl(larger))-bl(pl(smaller))), 0];
            on_peak = larger;
            off_peak = smaller;
          end
        end
        pow = min (alt(1,1), alt(2,1));
      else
        pow = 100 + min (alt(1,1), alt(2,1));
      end
    elseif abs (bj(pj(col)) - (bl(pl(larger)) - bl(pl(smaller)))) / bj(pj(col)) < 0.2
      if bj(pj(col)) == alt(1,1) ...
          && bl(pl(larger)) - bl(pl(smaller)) == alt(2,1)
        pow = min (alt(1,1), alt(2,1));
      else
        pow = 200 + min (alt(1,1), alt(2,1));
      end
    else
      pow = 300 + min (alt(1,1), alt(2,1));
    end
  else
    [gl] = max (rectl(:,3));
    [ggl] = max (rectl(rectl(:,3) < gl, 3));
    ggl = gl / ggl;
    gj  = NaN;
    ggj = NaN;
    pow = NaN;
    col = 1;
    larger = on_peak;
    smaller = off_peak;
  end

  % Update gap(level) to match the actual chosen power
  if on_peak ~= 2 || off_peak ~= 1
    gl = find_rect (al(min (pl([smaller,larger])):max (pl([smaller,larger]))));
    gl = gl(end);
  end
  
  th = min (alt(:,1));

  % estimate which "vampires" are due to a "big device" being on overnight.
  %[av, bv] = hist (vamp(days_use));
  v = user_data(1:size(user_data, 1)/4.8, days_use);
  [av, bv] = hist (rolling_min (v(:), 5));
  [pv, cntv, sepv, rectv] = secondPeak (av + 0.1*([av(2:end),0] + [0,av(1:end-1)]),5);
%  disp ([bv(pv); cntv]);
%{
  if rectv > evalin('base', 'rv')
    keyboard
  end
%}
  df = bv(pv(2)) - bv(pv(1));
  if th > df && 0.7 * th < df
    th = df;
  end
  
  % 
  cv = bsxfun (@minus, cv, vamp(days_use));
  
  a = 0; ACs;   % load known air conditioners from file
  %thresh = (gl > 50) || (gl > 30 && abs (bj(pj(col)) - (bl(pl(larger)) - bl(pl(smaller)))) / bj(pj(col)) < 0.03);
  gap_thresh = length (train_data) * 0.04;
  thresh =((gl      - 100 * abs (bj(pj(col)) - (bl(pl(larger)) - bl(pl(smaller)))) / bj(pj(col))) / 30 > 1);
  conf   = (gj + gl - 100 * abs (bj(pj(col)) - (bl(pl(larger)) - bl(pl(smaller)))) / bj(pj(col))) / 30;
  idx = find (a(:,1) == evalin('caller', 'i'));
% TODO: consider rect when reconciling power (or weighting level pairs)
manual_check = [26, 43, 44, 55, 67, 89, 95, 104, 109, 110];
  if 0 && ( ~isempty (idx) || any (evalin('caller', 'i') == manual_check))
    fprintf ('power: %g %g %g    thresh: %d %g (%g) %g (%g) -> %d, %d\n', ...
      a(idx,2), pow, th, a(idx,end), gj, ggj, gl, ggl, thresh, conf);
  end
%  disp (alt);
  ccv = cv;
  ccv(:) = -rolling_min (-rolling_min (cv(:)));
  
  if th < 0.8 * max (alt(:,1))
    col = find (bj > th, 1);
    if isempty (col)
      col = 1;
    end
    gj = find_rect (aj(1:col));
    gj = gj(end);
  end

  % Find rough time-dependence of power consumption of A/C
  acv = cv;
  acv(cv < th) = NaN;
  acv(sum (isfinite (acv), 2) < 4, :) = NaN;
  offset = medfilt1 (median (acv, 2, 'omitnan'));
  % TODO: make offset concave, or at least omit spikes

  offset1 = cv;
  offset1(:) = th;
  % Check if the non-constancy of offset is worth modelling,
  % or if it is probably just noise.
%  if sum (abs (diff (offset(~isnan(offset))))) < 2.1 * (max (offset) - min (offset))
  if sum (abs (diff (offset(~isnan(offset))))) < 2.2 * max (offset) - offset(find(isfinite(offset), 1)) - offset(find(isfinite(offset), 1, 'last'))
    offset = offset - min (offset);
    offset(isnan (offset)) = 0;
    offset1 = bsxfun (@plus, offset1, offset);
  end
  cor = zeros (size (user_data));
  cor(hrs_use, days_use) = -min(cv, (ccv > th) .* offset1);

  thresh = th;
%{
  if size (cv, 2) > 10 % && ~thresh
    figure(2); imagesc (cv + cor(hrs_use, days_use));  figure(3); imagesc (-cor(hrs_use, days_use));  figure(4); plot (cv);  figure(5); plot (cv + cor(hrs_use, days_use));  figure(6); hist (train_data, binsl);  if ~isempty (jump_times), figure(7); hist (jump_sizes, binsj); end; figure(1); imagesc (cv);
    keyboard
    if gj < 8 && gl < 20
      cor = 0 * cor;
    end
  end
%}
end

function [pwrs, rectj, ij, jj, aj, bj, pj, ipj, binsj] = common_jumps (user_data, sud)
  [jump_sizes, jump_times] = find_jumps (user_data);

  aj = [];
  bj = [];
  pj = [];
  ipj = [];
  ij = [];
  jj = [];
  rectj = 0;
  jbins = 15; mj = max (jump_sizes);
  binsj = [];
  if ~isempty (jump_times)
    jbins = jbins + 1;
    binsj = (0.5:1:jbins) * mj / jbins;
    [aj, bj] = hist (jump_sizes, binsj);
    %figure(4); hist(jump_sizes, binsj);
    % break ties by minor smoothing
    [pj, cj, ~, rectj] = secondPeak (aj + 0.1*([aj(2:end),0] + [0,aj(1:end-1)]),5);
    pj(cj==0) = NaN;	%Ignore non-peaks
    % Find "big" peak that isn't from small switches etc
    if abs (pj(2) - pj(1)) < max (bj(pj(isfinite (pj))))/3
      big = pj(3);
    else
      big = pj(2);
    end

    % Look for two peaks that add to a third: a sign of >= 2 A/Cs
    % The lowest and highest power peaks won't be individual ACs.
    % Exclude if neither they nor their sum explain the "big" peak
    [pp, ipj] = sort(pj);
    if ismember (pp(2)+pp(3), [pp(4), pp(5)]) && ismember (big, [pp(2), pp(3), pp(2)+pp(3)])
      if pp(2) == big
        pwrs(1) = bj(pp(2));
        pwrs(2) = bj(pp(3));
      else
        pwrs(1) = bj(pp(3));
        pwrs(2) = bj(pp(2));
      end
      ij=2; jj=3;
      %fprintf('Two ACs: %d & %d:  bj(%d)=%f  bj(%d)=%f\n', ii, jj, pp(ii), bj(pp(ii)), pp(jj), bj(pp(jj)));
    elseif pp(3)+pp(4) == pp(5) && ismember (big, [pp(4), pp(3), pp(4)+pp(3)])
      if pp(3) == big
        pwrs(1) = bj(pp(3));
        pwrs(2) = bj(pp(4));
      else
        pwrs(1) = bj(pp(4));
        pwrs(2) = bj(pp(3));
      end
      ij=3; jj=4;
      %fprintf('Two ACs: %d & %d:  bj(%d)=%f  bj(%d)=%f\n', ii, jj, pp(ii), bj(pp(ii)), pp(jj), bj(pp(jj)));
    elseif pp(2)+pp(4) == pp(5) && ismember (big, [pp(2), pp(4), pp(2)+pp(4)])
      if pp(2) == big
        pwrs(1) = bj(pp(2));
        pwrs(2) = bj(pp(4));
      else
        pwrs(1) = bj(pp(4));
        pwrs(2) = bj(pp(2));
      end
      ij=2; jj=4;
      %fprintf('Two ACs: %d & %d:  bj(%d)=%f  bj(%d)=%f\n', ii, jj, pp(ii), bj(pp(ii)), pp(jj), bj(pp(jj)));
    else
      % Just look for the largest peak of a "large" device
      % If the second largest peak is also small, skip it
      thresh = sud(ceil (length (user_data)*0.97))/3;
      pj = pj(isfinite (pj));
      while length (pj) > 2 && (abs (bj(pj(1))-bj(pj(2))) < thresh)
        pj = pj([1 3:length(pj)]);
      end
      if length (pj)>1
        pwrs(1) = bj(max (pj(1), pj(2)));
      else
        pwrs(1) = bj(pj(1));
      end
      pwrs(2) = 0;
    end
  else
    pwrs(1) = NaN;
    bj = 0;
    pj = 1;
  end
end