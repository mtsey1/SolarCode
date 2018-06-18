function vampires = find_vampires(c);
% TODO:  Vectorize.

    mc = min(c,[],2);

                % allow vampires to be off for 2 days in a row
    mcf = medfilt1(double(mc),5);       % (medfilt doesn't support single)
    vamps_off = (mcf > 2*mc);    % only consider  big  drops in power
    mcf(~vamps_off) = mc(~vamps_off);

    w_big   = 151;      % off-by-one error in location if window size is even
    w_small =  51;

    vampires    = -rolling_min(-rolling_min(mcf, w_small), w_small);

    guardlen = ceil(w_big/2);
    guard = median(vampires)*ones(guardlen,1);
    mce = [guard; mcf; guard];
    vampiresBigG = -rolling_min(-rolling_min(mce, w_big),   w_big);
    vampiresBigG = vampiresBigG(guardlen+1:end-guardlen);

    if skipNaN(@sum,vampiresBigG) > 0.9 * skipNaN(@sum,vampires)
%%figure(100); plot([mc, vampires, vampiresBigG]);
        vampires = vampiresBigG;
%%else
%%figure(100); plot([mc, vampiresBigG, vampires]);
    end
    % enforce vampires off on the days they are
    vampires = max (0, min(vampires, mc));
end
