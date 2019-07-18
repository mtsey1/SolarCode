todo_list = [1,3];

todo = zeros(10);
todo(todo_list) = 1;

bulkdatapath = 'C:/Users/llandrew/git/rsrch/NILM/UnitedEnergy/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of increased generation due to diffuse light from clouds
if todo(1)
    figure(1)

    if ~exist ('datasyd', 'var')
        load 100dualdata.mat
    end
    aa = -squeeze (datasyd(4,155:180,:))';
    sun   = min (aa(18:23, :)) > 0.3;
    nosun = min (aa(18:23, :)) <= 0.3;

    hold off
    plot ((1:48)/2, 2*aa(:, sun), '-')
    hold on
    plot ((1:48)/2, 2*aa(:, nosun), '--')
    hold off

    xlabel ('Time (hours since midnight)')
    xticks ([0, 6, 12, 18, 24])
    ylabel ('Generation (kW)')

    set (gcf, 'PaperPosition', [0.1 0 9 7]);
    print -depsc cloud_diffuse.eps;
end

if any (todo(2:3)) && ~exist ('state', 'var')
    load ([bulkdatapath, 'phase1_start3000.mat'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cloud estimates for a week
if todo(2)
    figure(2)
    hold off

    utilization = zeros(48, 365);
    utilization(9:39, :) = squeeze (state.solar_by_pc(2,:,:))';
    week = 1: 2*24*7;
    plot (week/(2*24), 2*utilization(3600 + week))

    ylabel ('Est. utilization')
    xlabel ('Day')
    axis ([0, 7, 0, 1]);

    set (gcf, 'PaperPosition', [0.1 0 9 7]);
    print -depsc cloud_est.eps;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cloud estimates for a week
if todo(3)
    figure(3)
    hold off
    
    if ~exist ('irradiance2013', 'var')
        load ('Data sets/irradiance2013.txt')
    end
    
    
    %days = [1:90, 270:365];
    %days = [91:270];
    days = 1:365;
    utilization2 = squeeze (state.solar_by_pc(2, :, :))';
    utilization3 = squeeze (state.solar_by_pc(30, :, :))';
    
    plot (irradiance2013(days), sum (utilization2(:, days)), 'o', ...
          irradiance2013(days), sum (utilization3(:, days)), 'x')
    xlabel ('Measured MJ/m^2')
    ylabel ('Inferred')

    set (gcf, 'PaperPosition', [0.1 0 9 7]);
    print -depsc generation_correlation.eps;
end
