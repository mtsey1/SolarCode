## Copyright (C) 2016 Lachlan Andrew
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} scripts/spikes (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Lachlan Andrew <lachlan@lachlan-VirtualBox>
## Created: 2016-12-20

function [cor] = spikes (user_data)
  ud = user_data(:);
  % Find spikes of duration 1
  sp = false (size (ud));
  heights = zeros (size (ud));
  mx = max (ud(1:end-2), ud(3:end));
  heights(2:end-1) = ud(2:end-1) - mx;
  sp(2:end-1) = (heights(2:end-1) > 3 * (mx - min (ud(1:end-2), ud(3:end))));
  heights(~sp) = 0;
  cor = heights;
  
  [a, b] = hist (heights(sp), 40);
  figure(1); imagesc (user_data);
  figure(2); hist (heights(sp), 40);
  
  % Find spikes of durations > 1
  
endfunction
