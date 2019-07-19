% Find postcodes "near" each Victorian postcode
thresh = 0.02;	% How "near".  Latitude and longitude must both be within thresh

load postcodes_duplicates.txt
pd = postcodes_duplicates;
unique = setdiff(pd(:,1), []);
for i = 1:length(unique)
    neighbours(i,1) = unique(i);
    %if unique(i) == 3764	% strange postcode, with many geocodes
        %continue
    %end
    me = find(pd(:,1) == unique(i));
    if isempty(me)
        continue
    end
    lo_lat = min(pd(me, 2));
    hi_lat = max(pd(me, 2));
    lo_lon = min(pd(me, 3));
    hi_lon = max(pd(me, 3));
    nb = find((pd(:,2) > lo_lat-thresh) .* ...
              (pd(:,2) < hi_lat+thresh) .* ...
              (pd(:,3) > lo_lon-thresh) .* ...
              (pd(:,3) < hi_lon+thresh));

    % If too many neighbours, don't allow "tolerance" on lat and long
    if length(nb) > 10
	nb1 = find((pd(:,2) > lo_lat) .* ...
		   (pd(:,2) < hi_lat) .* ...
		   (pd(:,3) > lo_lon) .* ...
		   (pd(:,3) < hi_lon));
        if length(nb1 > 3)
	    nb = nb1;
	end
    end

    nb = setdiff(pd(nb,1), pd(me,1));	% get unique nearby postcodes

    % If no postcode nearby, find nearest
    if isempty(nb)
	not_me = setdiff(1:size(pd,1), me);
        [dummy nb] = min((pd(not_me,2)-lo_lat).^2 + (pd(not_me,3)-lo_lon).^2);
	nb = pd(not_me(nb),1);
    end
    neighbours(i,1+(1:length(nb))) = nb;
end
