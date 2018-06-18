
% To be called at the end of calc_cor_solar
cs = reshape (sum (cor_solar), [365, 48])';
light_hrs = s.dark_start - s.dark_end + 1;
cs = full (cs(s.dark_end:s.dark_start,:));
cf = full (squeeze (sum (s.capFactor))');
rr = cs ./ cf;

[b, dark_until] = max (rr(1:end/2, :) < 0.10);   % first index s.t. rr < 0.075
[~, dark_until(b == 0)] = min (rr(1:end/2, b == 0));
dark_until = round (medfilt1 (dark_until));
rr(mod (1:numel (rr), light_hrs) <= repelem (dark_until, light_hrs)) = 0;

r = rr(end:-1:end/2, :);
[b, dark_from]  = max (r < 0.20 & r > 0);
[~, dark_from(b == 0)] = min (rr(end:-1:end/2, b == 0));
dark_from = round (medfilt1 (dark_from));
rr(mod (numel (rr):-1:1, light_hrs) <= repelem (dark_from, light_hrs)) = 0;

rr(rr > 1) = 0;