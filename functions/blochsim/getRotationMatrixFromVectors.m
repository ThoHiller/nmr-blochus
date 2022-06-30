function R = getRotationMatrixFromVectors(A,B)
%getRotationMatrixFromVectors calculates rotation matrix R to rotate A into B
%
% Syntax:
%       getRotationMatrixFromVectors(A,B)
%
% Inputs:
%       A - start vector
%       B - end vector
%
% Outputs:
%       R - 3x3 rotation matrix to rotate A into B
%
% Example:
%       R = getRotationMatrixFromVectors([1 0 0]',[0 0 1]')
%       yields R = 0  0 -1
%                  0  1  0
%                  1  0  0
%       so that R*[1 0 0]' = [0 0 1]'
%
% Other m-files required:
%       none
%
% Subfunctions:
%       none
%
% MAT-files required:
%       none
%
% See also: BLOCHUS
% Author: Thomas Hiller
% email: thomas.hiller[at]leibniz-liag.de
% License: GNU GPLv3 (at end)

%------------- BEGIN CODE --------------

% normalize both vectors
A = A(:)./norm(A);
B = B(:)./norm(B);
% cross product of both vectors
c = cross(A,B);
% angle between A and B
alpha = acos(dot(A,B));
% check if A and B are parallel / antiparallel
if abs(sum(c))<1e-128 && (alpha==0 || alpha==pi)
    % check if A == B (alpha=0 -> parallel)
    if alpha==0
        % in that case the rotation matrix is obviously identity
        R = eye(3);
    else % A == -B (antiparallel)
        R = -eye(3);
    end
else
    % skew-symmetric cross-product
    ssc = [ 0  -c(3) c(2);
           c(3)  0  -c(1);
          -c(2) c(1)  0 ];
    % rotation matrix R
    R = eye(3) + ssc + (ssc/norm(c))^2*(1-dot(A,B));
end

return

%------------- END OF CODE --------------

%% License:
% GNU GPLv3
%
% BLOCHUS
% Copyright (C) 2019 Thomas Hiller
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
