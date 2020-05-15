function onCheckRDP(src,~)
%onCheckRDP updates the checkboxes that activate relaxation during pulse
%(RDP) or relaxation during switch-off (RDS)
%
% Syntax:
%       onCheckRDP(src)
%
% Inputs:
%       src - handle of the calling object
%
% Outputs:
%       none
%
% Example:
%       onCheckRDP(src)
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

% get GUI handle
fig = ancestor(src,'figure','toplevel');

if ~isempty(fig) && strcmp(get(fig,'Tag'),'BLOCHUS')
    % get GUI data
    gui  = getappdata(fig,'gui');
    data = getappdata(fig,'data');
    
    switch get(src,'Tag')        
        case 'PrePolRDS'
            % update the RDS switch
            data.prepol.RDS = get(src,'Value');
        case 'PulseRDP'            
            % update the RDP switch
            data.pulse.RDP = get(src,'Value');
    end
    
    % update the data inside the GUI
    setappdata(fig,'data',data);
else
    warndlg({'onCheckRDP:','There is no figure with the BLOCHUS Tag open.'},...
        'BLOCHUS error');
end

end

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