function [power, cor, freq] = periodic (data)
  ac_len = 80;
  seem_periodic = zeros (ac_len, size (data,1));
  interp_periodic = zeros (ac_len, size (data,1));
  quality = zeros (1, size(data,1));
  for i = 1:size(data,1)
    dd = squeeze(data(i,:,:))';
    if any (isnan (dd(:)))
      q = 0;
    else
      diffs = diff (dd(1:end-1)) - diff (dd(2:end));        % remove trends
      %diffs = diffs(isfinite (diffs));
      diffs(abs (diffs) > 0.5 * skipNaN (@mean, abs (diffs))) = NaN;
      ac = autocorr (diffs, 1:ac_len);
      v = var (diffs(isfinite (diffs)));
      q = sumsq (ac(2:end)) / (v .^ 2);
      seem_periodic(:, i) = ac / v;
%      interp_periodic(:, i) = ac2 / (v.^2);
    end
    quality(i) = q;


  end

  [~, idx] = sort (-quality);

  for i = 1:idx
    figure(1)
    plot (seem_periodic(:, idx(i)));

    scale = 10;
    ac2 = autocorr (interp (seem_periodic(:, idx(i)), scale), 1:ac_len);
    first_neg = find (ac2 < 0, 1)
    next_pos = find (ac2(first_neg+1:end) > 0, 1) + first_neg
    next_neg = find (ac2(next_pos +1:end) < 0, 1) + next_pos
    [~, period] = max (interp (ac2(next_pos:next_neg), scale))
    period = (period-1)/scale + next_pos
    period = (period-1)/scale

    figure(2)
%    plot (interp_periodic(:, idx(i)));
    plot (ac2);

    figure(3)
    imagesc (squeeze (data(idx(i), :, :))');


    dd = squeeze(data(idx(i),:,:))';
    aa = dd(isfinite (dd(:)));
    bb = medfilt1 ((aa(2:2:end-2) - (aa(1:2:end-2) + aa(3:2:end))/2), 5);
    cc = fft (bb);
    cc(abs (cc) < 0.5 * max (abs (cc))) = 0;
    ee = real (ifft (cc));
    cor = zeros (size (dd));
    cor (end-2*length(ee)  :2:end-2) = -ee/2;
    cor (end-2*length(ee)+1:2:end-1) = +ee/2;
    ff = dd + cor;

    figure(4);
    imagesc (ff);

    z = fft (dd);
    y = z ./ abs (z);
    figure(5)
    plot (abs (diff (y (13,:))));       % for period-4 patterns

    figure(6)
    imagesc (abs (z));

    len_bb = length (bb(:));
    len_pad = 256 * (ceil (len_bb / 256) - len_bb / 256);
    cc = reshape ([bb(:); zeros(len_pad,1)], [256, (len_bb+len_pad)/256]);
keyboard
    cc = fft (cc);
    cc(abs (cc(:)) < 0.5 * max (abs (cc(:)))) = 0;
    ee = real (ifft (cc));
    ee = ee(1:len_bb);
    cor = zeros (size (dd));
    cor (end-2*length(ee)  :2:end-2) = -ee/2;
    cor (end-2*length(ee)+1:2:end-1) = +ee/2;
    ff = dd + cor;

    figure(8);
    imagesc (ff);

    figure(7);
    imagesc (abs (z));

idx(i)
keyboard
  end

%%  54

%  spec = fft (ac);
%  spec = abs (spec (ceil(end/2):end));
%
%  ac2 = autocorr (ac(6:end-6), 1:50);

%  figure (3);
%  plot (ac2);
%
%  figure(4);
%  plot (spec);
%
%  ac2 = autocorr (spec, 1:50);
%  figure(5);
%  plot (ac2);

end

function [power, cor, freq] = periodic_old (data)
%  dd = squeeze(data (1,:,1:8))';
  ac_len = 20;
  seem_periodic = zeros(ac_len, 10);
  candidates = [];
  idx = 0;
  for i = 1:size(data,1)
    dd = squeeze(data(i,:,:))';

    diffs = diff (dd(1:end-1)) - diff (dd(2:end));        % remove trends

    diffs = diffs(isfinite (diffs));
    ac = autocorr (diffs, 1:ac_len);
    ac2 = ac(2:end) / var (diffs);

    second = ac2;
    [~, mx] = max (abs (ac2));
    second(mx) = 0;
    second = max (second);

    keep = true;
    keep = keep && (second > 0.1);
    %keep = keep && (max(abs(ac2(4:8))) > 0.5);
    %keep = keep && (max(abs(ac2(9:15))) > 0.4);
    %keep = keep && (mean(abs(ac2)) > 0.3);
    if keep
      idx = idx + 1;
      if idx > size (seem_periodic, 2)          % resize by doubling
        seem_periodic(1, 2*idx) = 0;
      end
      seem_periodic(:, idx) = ac / var (diffs);
      candidates(idx) = i;
    end
  end

  figure (2);
  plot (seem_periodic(:,1:idx));
  idx

  for i = 1:idx
    figure(1)
    plot (seem_periodic(2:end,i));

    figure(3)
    imagesc (squeeze (data(candidates(i), :, :))');
candidates(i)
keyboard
  end

%%  54

%  spec = fft (ac);
%  spec = abs (spec (ceil(end/2):end));
%
%  ac2 = autocorr (ac(6:end-6), 1:50);

%  figure (3);
%  plot (ac2);
%
%  figure(4);
%  plot (spec);
%
%  ac2 = autocorr (spec, 1:50);
%  figure(5);
%  plot (ac2);

end

%%%%%%%%%%%%%%%%%%%%
function ac = autocorr (x, offsets)
  offsets = min (abs (offsets), length (x)-1);

  ac = zeros (size (offsets));
  for i = 1:length(offsets)
    x1 = x(1:end-offsets(i));
    x2 = x(offsets(i)+1:end);
    x1 = x1 - skipNaN (@mean, x1);
    x2 = x2 - skipNaN (@mean, x2);
    ac(i) = skipNaN (@mean, x1 .* x2);
  end
end

%%%%%%%%%%%%%%%%%%%%
function ac = autocorr_NaNs_mess_up (x, offsets)
  offsets = min (abs (offsets), length (x)-2);

  ac = zeros (size (offsets));
  for i = 1:length(offsets)
    ac(i) = cov (x(1:end-offsets(i)), x(offsets(i)+1:end));
  end
end
