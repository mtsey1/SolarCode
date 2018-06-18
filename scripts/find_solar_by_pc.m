function [s] = find_solar_by_pc (s, meta)
    
  fprintf('Calculating solar corrections\n');
  solar_by_pc(length(meta.pclist),meta.Days,s.dark_start-s.dark_end+1)=0;
  new_solar_by_pc = solar_by_pc;
  for i = 1:meta.Days
    for j = 0:(s.dark_start-s.dark_end)
      for k = 1:length(meta.pclist)
        % Index into solar_users  of those in post code k
        u = find(s.postcode(s.solar_users) == meta.pclist(k));
        % May have problems if some have very small cap factor
        % Only consider those that can generate >100W now
        if isempty(u)
            continue
        end
        u = u(s.capFactor(u,i,j+1) ~= 0);       % Ignore non-generating panels
        if isempty(u)
            continue
        end
        v = u(s.capFactor(u,i,j+1) .* s.solar_cap(u)' > 0.1);
        %if (i == 1 && j == 0) fprintf('Postcode %d has %d solar users, with %d > 100W at day %d hour %d\n', pclist(k), length(u), length(v), i, (s.dark_end+j)/2); end
        if ~isempty(v)
            u = v;
        end
        
        % solar_by_pc  is the number of kwh per 30-min slot that would
        % be generated by a 1kW solar system if the sun is perpendicular
        solar_by_pc(k,i,j+1) = min(max(s.seen(u,i,j+1)./s.capFactor(u,i,j+1)),0.5);
      end
    end
  end
  
  % Sanity checks.
  % If solar_by_pc varies greatly one day for some postcode, should for nearby
  % solar_by_pc should be similar for nearby postcodes
  
  % Make new_solar_by_pc the max of all neighbouring postcodes.
  for k = 1:length(meta.pclist)
    neighbours = meta.postcode_neighbours(meta.postcode_neighbours(:,1) ...
                                          == meta.pclist(k),:);
    [~, neighbours] = ismember(neighbours, meta.pclist);    % postcode->idx
    neighbours = neighbours(neighbours>0);
    if ~isempty(neighbours)
        new_solar_by_pc(k,:,:) = max(solar_by_pc(neighbours,:,:), [],1);
    end
  end
  s.solar_by_pc = new_solar_by_pc;
end
