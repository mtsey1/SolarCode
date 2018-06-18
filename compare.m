for i = 1:length (tmp)
  a = tmp{i};
  b = tmp1{i};
  if isempty (a) ~= isempty (b)
    fprintf ('%4d mismatch %d %d\n', i, isempty (a), isempty (b));
  elseif ~isempty (a)
    aa = [a.on_off]';
    bb = [b.on_off]';
    aaa = aa(:, [1,3])+1;
    bbb = bb(:, [1,3]);
    if any (any (aaa ~= bbb))
      fprintf ('%4d 1,3 error\n', i);
      c = find (any (aaa ~=bbb, 2));
      %[aaa(c,:), bbb(c,:)]
      %keyboard
    end
  end
end
