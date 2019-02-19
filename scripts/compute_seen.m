function [seen_summer, seen_winter] = compute_seen(data, data_no_vamp, r1, r2, meta)

    % Within r1 and r2, find the three biggest samples,
    % and the biggest without excluding vampires.
    %
    % Find top three in r1 and r2
    d = sort (data_no_vamp(1, r1, :));
    d = d(1, 1:min(3, size (d,2)), :);
    d = d * -meta.SamPerDay/24;
    nv = min (data(1, r1, :), [], 2) * -meta.SamPerDay/24;
    a = [d, nv];
    seen_summer = reshape(a, [size(a,2), 1, size(d,3)]);
    
    %For visualising
%     figure;
%     plot(linspace(0,23,48),[squeeze(d)' squeeze(nv)]);
%     xlabel('Time of day');
%     ylabel('Generation - Consumption');
%     title(sprintf('House %d - Summer',i));
%     legend({'Largest gen (no vampire)','2nd largest gen (no vampire)','3rd largest gen (no vampire)','Largest gen (with vampire)'},'Location','best');

    d = sort (data_no_vamp(1, r2, :));
    d = d (1, 1:min (3, size (d,2)), :);
    d = d * -meta.SamPerDay/24;
    nv = min (data(1, r2, :), [], 2) * -meta.SamPerDay/24;
    a = [d, nv];
    seen_winter = reshape(a, [size(a,2), 1, size(d,3)]);

end