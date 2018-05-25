## Copyright (C) 2013 Carnë Draug <carandraug+dev@gmail.com>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>.

## This private function is to be common to all functions that have an option
## for the interpolation method to use. For matlab compatibility, bicubic and
## bilinear are accepted which are the same as cubic and linear. This does not
## actually check if the method is valid, we leave that to interp2. The reason
## is that if interp2 implements a new method, all this functions will
## automatically work with it.

function method = interp_method (method)
  method = tolower (method);
  switch method
    case "bicubic",   method = "cubic";
    case "bilinear",  method = "linear";
    case "triangle",  method = "linear";  # interpolation kernel
  endswitch
endfunction
