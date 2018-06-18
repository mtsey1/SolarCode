function m = mod1 (value, modulus)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Modulo operator, returning values in [1, modulus], not [0, modulus-1].
  m = mod (value-1, modulus) + 1;
end
