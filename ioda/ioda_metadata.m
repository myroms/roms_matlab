function M = ioda_metadata (varargin)

%
% IODA_METADATA:  Sets metadata structure for IODA-type NetCDF4 files
%
% M = ioda_metadata (report)
%
% This function initializes the metadata structure, M, for IODA-enhanced
% NetCDF-4 files. It sets the names of the data assimilation control
% variables, the area-averaged, and the time-averaged parameters to
% empty values. Users will need to assign these values at a later stage.
%  
% The function is designed to be used before calling the "roms2ioda.m"
% converter, as the resulting structure M is passed as an argument to
% that function.
%  
% On Input:
%
%    report    Switch to report metadata (OPTIONAL)
%
% On Output:
%
%    M         IODA enhanced NetCDF-4 file metadata structure
%                (struct array)
%
%                M(:).name           variable short name
%                M(:).half_length    area-averaged half-length (km) scale
%                M(:).time_window    time-averaged window (hours)
%                M(:).ioda_vname     IODA NetCDF-4 variable name
%                M(:).standard_name  variable standard name
%
%                M(:).name = 'SSH', 'SST', 'SSS', 'uv_CODAR', 
%                            'ptemp', 'temp', 'salt'
%                 
% USAGE:
%
%       M = ioda_metadata;
%  or
%       M = ioda_metadata(true);
%
%  To set the area-averaged and time-averaged scales for SSH use:
%
%       M(strcmp({M.name}, 'SSH')).half_length = 30;       % km
%       M(strcmp({M.name}, 'SSH')).time_window = 36;       % hours
%
%  To set a time-averaged window of 24 hours for uv_CODAR velocities use:
%
%       M(strcmp({M.name}, 'uv_CODAR')).time_window = 24;  % hours
%
%  To display the updated values use:
%
%       disp( struct2table(M) );
%

% git %Id$
%=======================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                               %
%    Licensed under a MIT/X style license                               %
%    See License_ROMVariables.txt                   Hernan G. Arango    %
%=======================================================================%

% Initialize.

switch numel(varargin)
  case 0
    report = false;
  case 1
    report = varargin{1};
end

% Define structure array where is element is associated with a variable
% in ROMS data assimilation control vector.

M(1:7) = struct('name'          , [],                                 ...
                'half_length'   , [],                                 ...
                'time_window'   , [],                                 ...
                'ioda_vname'    , [],                                 ...
                'standard_name' , []);

% Initialize structure values.

M(1).name          = 'SSH';
M(1).ioda_vname    = 'absoluteDynamicTopography';
M(1).standard_name = 'absolute_dynamic_topography';

M(2).name          = 'SST';
M(2).ioda_vname    = 'seaSurfaceTemperature';
M(2).standard_name = 'sea_surface_temperature';

M(3).name          = 'SSS';
M(3).ioda_vname    = 'seaSurfaceSalinity';
M(3).standard_name = 'sea_surface_salinity';

M(4).name          = 'uv_CODAR';
M(4).ioda_vname    = {'waterZonalVelocity',                           ...
                      'waterMeridionalVelocity'};
M(4).standard_name = {'eastward_sea_water_velocity',                  ...
                      'meridional_sea_water_velocity'};

M(5).name          = 'ptemp';
M(5).ioda_vname    = 'waterPotentialTemperature';
M(5).standard_name = 'sea_water_potential_temperature';

M(6).name          = 'temp';
M(6).ioda_vname    = 'waterTemperature';
M(6).standard_name = 'sea_water_temperature';

M(7).name          = 'salt';
M(7).ioda_vname    = 'salinity';
M(7).standard_name = 'sea_water_salinity';

% Report.

if (report)
  nvars = length(M);

  disp(blanks(1));
  disp('IODA NetCDF-4 file Metadata Structure:');
  disp(blanks(1));
  
  for i = 1:nvars
    disp(M(i))
  end
end

return





