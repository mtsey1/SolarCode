old = load ('aircon_checkpoint0');
new = load ('aircon_checkpoint1');
load ../data/aircon_labelled.mat

labelled = zeros (size (new.cor_ac));
for n = 1:length (aircon_labelled)    % Can't vectorise as some spare etc.
  labelled(:, :, n) = aircon_labelled(n).cor;
end
old_err = (old.cor_ac ~= 0) ~= (labelled ~= 0);
new_err = (new.cor_ac ~= 0) ~= (labelled ~= 0);

known = ([aircon_labelled.done] == true);
known(2) = false;

figure(1); imagesc (sum(old_err(:, :, known), 3) -sum(new_err(:, :, known), 3))

only_old = old_err(:, :, known) & ~ new_err(:, :, known);
only_new = new_err(:, :, known) & ~ old_err(:, :, known);
[sum(only_old(:)), sum(only_new(:))]

extra_positive = (new.cor_ac(:, :, known) ~= 0) & ~(old.cor_ac(:, :, known) ~= 0);
extra_negative = (old.cor_ac(:, :, known) ~= 0) & ~(new.cor_ac(:, :, known) ~= 0);
[sum(extra_positive(:)), sum(extra_negative(:))]

for i = 1:length(aircon_labelled)
  a = (new.cor_ac(:, :, i) ~= 0);
  b = (old.cor_ac(:, :, i) ~= 0);
  if (aircon_labelled(i).done == true) && any (a(:) ~= b(:))
    me_year = squeeze (new.corrected(i, :, :))';
    c = (aircon_labelled(i).cor ~= 0);
    figure(2); show_aircon (a, me_year, [1, 1, 1]);
    figure(3); show_aircon (b, me_year, [1, 1, 1]);
    figure(4); show_aircon (c, me_year, [1, 1, 1]);
    %figure(4); show_aircon (a ~= b, squeeze (new.corrected(i, :, :))', [1, 1, 1]);
    figure(5); imagesc (a - c);
    figure(6); imagesc (b - c);
    figure(7); imagesc ((a ~= c) - (b ~= c));
    [x, y] = ind2sub (size (a), find(a(:) ~= b(:), 1));
    fprintf ('i = %d row = %d col = %d\n', i, x, y);
    keyboard;
  end
end

for n = find ([aircon_labelled.done] == true)
  nn = (new.cor_ac(:, :, n) ~= 0);
  oo = (old.cor_ac(:, :, n) ~= 0);
  rr = (aircon_labelled(n).cor ~= 0);
  
  if any (nn(:) ~= oo(:))
    figure(1); imagesc (nn - rr);
    figure(2); imagesc (oo - rr);
    figure(3); imagesc (oo - nn);
    figure(4); imagesc (aircon_labelled(n).cor);
    figure(5); imagesc (pr);
    figure(6); imagesc (squeeze (corrected(n, :, :))');
    
    keyboard
  end
end