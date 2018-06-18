function cancelled = plot_rectangles (rectangles, meta)
    cancelled = zeros (meta.Days, meta.SamPerDay);
    for i = 1:size (rectangles, 1)
        r = rectangles(i);
        st = r.on_off(1);
        en = r.on_off(3);
        b = r.burst(1:end-1);
        missed = r.missed;
        ht = r.power;
        cancelled(st+1:en,b) = cancelled(st+1:en,b)-ht;
        cancelled(missed, b) = cancelled(missed, b)+ht;
    end
end
