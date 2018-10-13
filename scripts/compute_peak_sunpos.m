function [spos, wpos, bpos] = compute_peak_sunpos(data, data_no_vamp, meta, s)
  if meta.hemisphere=='south'
    smr = [meta.January, meta.February, meta.November, meta.December];
    wtr = [meta.June,meta.July];
    broad = [meta.May, meta.June, meta.July, meta.August, meta.September];
  else
    smr = [meta.May, meta.June, meta.July, meta.August, meta.September];
    wtr = [meta.January,meta.December];
    broad = [meta.January, meta.February, meta.November, meta.December];      
  end
  if size(size(data)) ==[1,2]
      data=reshape(data,[1,365,48]); 
      data_no_vamp=reshape(data_no_vamp,[1,365,48]); 
  end
  % Record maximum in summer, winter and a broadly-defined winter,
  % along with the day of each.
  sz = [size(data,1), size(s.SunPos, 2)];
  sz_year = [size(data,2), size(s.SunPos, 2)];
  hrs = repmat (1:size (s.SunPos, 2), size (data, 1), 1);
  solar_range = hrs(1,:) + s.solar_start - 1;

  % Sun positions on peak generation days in summer and winter
  [~, spos] = min (data_no_vamp(:, smr, :), [], 2);
  spos = s.SunPos (sub2ind (sz_year, smr(squeeze (spos(:, :, solar_range))), hrs));
  spos = reshape (spos, sz);

  [~, wpos] = min (data_no_vamp(:,wtr,:), [], 2);
  wpos = s.SunPos (sub2ind (sz_year, wtr(squeeze (wpos(:, :, solar_range))), hrs));
  wpos = reshape (wpos, sz);


  %broad_winter = reshape (max (-meta.SamPerDay/24*bsxfun (@minus,data(:,broad,:), vampires(:,broad)),[],2), [size(data,1), meta.SamPerDay]);
  [~, bpos] = min (data_no_vamp(:,broad,:), [], 2);
  bpos = s.SunPos (sub2ind (sz_year, broad(squeeze (bpos(:, :, solar_range))), hrs));
  bpos = reshape (bpos, sz);

end