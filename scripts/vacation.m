
data = [];
for i = [1 4 7 10];
  filename = sprintf('month%d.txt', i);
  A = load(filename);
  data = [data; A(:, [1 2 51]) ];
end

% Rename users contiguously.
local_users = 0;
for UE_user = 1:size(data,1)
  matches = find (data(:,1) == UE_user);
%  if length(matches) == 365
  if length(matches) == 123
      local_users = local_users + 1;
      fullYear(local_users,:) = data(matches,3)';
      s =  sum(fullYear(local_users,:));
      if (s ~= 0)
	  features(local_users,:) = fullYear(local_users,:) / s;
      else
	  features(local_users,:) = fullYear(local_users,:);
      end
      local_to_UE_user(local_users) = UE_user;
      UE_to_local_user(UE_user) = local_users;
  elseif length(matches) ~= 0
      printf('User %d has %d entries\n', UE_user, length(matches));
  end
end
