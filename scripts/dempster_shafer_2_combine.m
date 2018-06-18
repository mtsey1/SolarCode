function answer = dempster_shafer_2_combine (arg1, arg2)
 % Demster-Shafer of A, B, A-or-B
 % (1,:) = yes, (2,:) = no, (3,:) = either
  % yes
  answer(1,:) = arg1(1,:) .* arg2(1,:) +  arg1(1,:) .* arg2(3,:) ...
                           + arg1(3,:) .* arg2(1,:);
  % no
  answer(2,:) = arg1(2,:) .* arg2(2,:) +  arg1(2,:) .* arg2(3,:) ...
                           + arg1(3,:) .* arg2(2,:);
  % either
  answer(3,:) = arg1(3,:) .* arg2(3,:);

  one_minus_k = answer(1,:) + answer(2,:) + answer(3,:);
  answer = bsxfun (@rdivide, answer, one_minus_k);

%  yes    = arg1.yes * arg2.yes + arg1.yes * arg2.either ...
%	                       + arg1.either * arg2.yes;
%  either = arg1.either * arg2.either;
%  no     = arg1.no * arg2.no + arg1.no * arg2.either ...
%	                     + arg1.either * arg2.no;
%
%  one_minus_k = yes + either + no;
%  answer.yes    = yes    / one_minus_k;
%  answer.either = either / one_minus_k;
%  answer.no     = no     / one_minus_k;
end
