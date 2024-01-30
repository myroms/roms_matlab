function create_ioda_obs(S, file)

%
% CREATE_IODA_OBS:  Creates an IODA observation NetCDF4 file
%
% create_ioda_obs(S, file)
%
% This function creates an IODA observation file using specified structure,
% S, for data assimilation within JEDI.
%
% Currently, it creates a NetCDF-4 file with the following groups:
%
%   group: MetaData {
%     variables:
%       int64 dateTime(Location) ;
%       string date_time(Location) ;
%       float depth(Location) ;
%       float latitude(Location) ;
%       float longitude(Location) ;
%       int provenance(Location) ;  
%       int sequenceNumber(Location) ;
%       string variables_name(nvars) ;
%   }
%
%   group: ObsError {
%     variables:
%       float myObsVariableName(Location) ;
%   }
%
%   group: ObsValue {
%     variables:
%       float myObsVariableName(Location) ;
%   }
%
%   group: PreQC {
%     variables:
%       float myObsVariableName(Location) ;
%   }
%  
% On Input:
%
%    S        Observations file creation parameters (struct):
%
%               S.ncfile            NetCDF file name (string)
%               S.nlocs             number of observations
%               S.nvars             number of variables
%               S.ncvname(:)        IODA NetCDF4 variable name
%               S.variable_names(:) IODA variable standard name
%               S.units(:)          variables units
%               S.depth             depth of observation
%               S.provenance        observation provenance
%               S.TimeIODA          IODA reference time YYYYMMDDHH
%               S.DateIODA          IODA reference date string
%
%    file     Output observation file name (string, OPTIONAL). If provided,
%               it creates this file instead the one in S.ncfile.
%
% WARNING:
% =======
%
% This function works for Matlab version 2022b or higher because it uses
% variable type NC_STRING (NetCDF4 files only). For example:
%
% group: MetaData {
%
%   string date_time(Location) ;
%          date_time:long_name = "ISO 8601 UTC date and time string" ;
%   string variables_name(nvars) ;
%          variables_name:long_name = "observation UFO/IODA standard name" ;
%
% data:  
%
%   date_time = "2004-01-03T14:21:00Z", "2004-01-03T14:21:00Z", ... ;
%
%   variables_name = "absolute_dynamic_topography" ;
% }
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2024 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

if (nargin < 2)
  if (isfield(S,'ncfile'))
    ncfile=S.ncfile;
  else
    error(['CREATE_IODA_OBS - Cannot find file name field: ncfile, ',   ...
          'in structure array S']);
  end
else
  ncfile=file;
end

%--------------------------------------------------------------------------
% Initialize.
%--------------------------------------------------------------------------

% Get Matlab version.

Mversion = version('-release');
Vyear    = sscanf(Mversion, '%i');

% Set indices for local dimension structure vector element.

nlocs = 1;           % "Location" dimension index
nvars = 2;           % "nvars"    dimension index
nstr  = 5;           % "nstring"  dimension index

% Set indices for NetCDF/HDF5 Group local structure vector element.

Meta = 1;            % MetaData Group
ObsE = 2;            % ObsError Group
ObsV = 3;            % ObsValue Group
PreQ = 4;            % PreQC    Group

% Set variable representation kind.

nc_int    = netcdf.getConstant('nc_int');        % integer
nc_int64  = netcdf.getConstant('nc_int64');      % 64-bit integer
nc_real   = netcdf.getConstant('nc_float');      % floating-point
nc_string = netcdf.getConstant('nc_string');     % string type Matlab 2022b

% Control switches.

do_depth = false;              % surface observations
if (isfield(S, 'depth'))
  if (~isempty(S.depth))
    do_depth = true;           % sub-surface observations
  end
end  
  
do_provenance = false;         % do not include observation origin
if (isfield(S, 'provenance'))
  if (~isempty(S.provenance))
    do_provenance = true;      % include observation origin
  end 
end
  
%--------------------------------------------------------------------------
% Set NetCDF4 file dimensions.
%--------------------------------------------------------------------------

D(nlocs).name = 'Location';     D(nlocs).size = S.nlocs;
D(nvars).name = 'nvars';        D(nvars).size = S.nvars;

%--------------------------------------------------------------------------
% Create NetCDF4 observation file.
%--------------------------------------------------------------------------

disp(' ');
disp(['*** Creating observations file:  ', ncfile]);

mode = netcdf.getConstant('NETCDF4'); 
mode = bitor(mode, netcdf.getConstant('NC_CLOBBER'));

ncid = netcdf.create(ncfile, mode);

%--------------------------------------------------------------------------
% Define NetCDF4 file dimensions.
%--------------------------------------------------------------------------

for dim = 1:length(D)
  name = char(D(dim).name);
  size = D(dim).size;
  D(dim).did = netcdf.defDim(ncid, name, size);
end

%--------------------------------------------------------------------------
% Define NetCDF4 file global attributes.
%--------------------------------------------------------------------------

varid = netcdf.getConstant('nc_global');

netcdf.putAtt(ncid, varid, '_ioda_layout', 'ObsGroup');
netcdf.putAtt(ncid, varid, '_ioda_layout_version', int32(3));
netcdf.putAtt(ncid, varid, 'odb_version', int32(1));
netcdf.putAtt(ncid, varid, 'date_time', S.TimeIODA);
netcdf.putAtt(ncid, varid, 'datetimeReference', S.DateIODA);

string='Native ROMS 4D-Var observations file converted to IODA';
netcdf.putAtt(ncid, varid, 'description', string);

netcdf.putAtt(ncid, varid, 'sourceFiles',  S.source);

history=['Created from Matlab script: ', mfilename, ' on ', date_stamp];
netcdf.putAtt(ncid, varid, 'history', history);

%--------------------------------------------------------------------------
% Define NetCDF4 file variables: same as dimension names.
%--------------------------------------------------------------------------

for i = 1:length(D)
  Vname = char(D(i).name);
  switch Vname
    case 'Location'
      D(i).vid = netcdf.defVar(ncid, Vname, nc_int, D(i).did);
      netcdf.putAtt(ncid, D(i).vid, 'suggested_chunck_dim',             ...
                    int32(D(i).size));
    case 'nvars'
      D(i).vid = netcdf.defVar(ncid, Vname, nc_int, D(i).did);
      netcdf.putAtt(ncid, D(i).vid, 'suggested_chunck_dim',             ...
                    int32(100));
  end
end

%--------------------------------------------------------------------------
% Define NetCDF4 files Groups: 'MetaData', 'ObsError', 'ObsValue', and
%                              'PreQc',
%--------------------------------------------------------------------------

G(Meta).gid = netcdf.defGrp(ncid, 'MetaData');
G(ObsE).gid = netcdf.defGrp(ncid, 'ObsError');
G(ObsV).gid = netcdf.defGrp(ncid, 'ObsValue');
G(PreQ).gid = netcdf.defGrp(ncid, 'PreQC');

%--------------------------------------------------------------------------
% Define variables in the 'MetaData' Group.
%--------------------------------------------------------------------------

MetaVars = {'dateTime', 'latitude', 'longitude', 'sequenceNumber',      ...
            'variables_name'};
MetaVars = [MetaVars, 'date_time'];
if (do_depth)
  MetaVars = [MetaVars, 'depth'];
end
if (do_provenance)
  MetaVars = [MetaVars, 'provenance'];
end
G(Meta).vars = sort(MetaVars);

for i = 1:length(G(Meta).vars)
  Vname = char(G(Meta).vars(i));
  switch Vname
    case 'dateTime'
      G(Meta).vid(i) = netcdf.defVar(G(Meta).gid, Vname, nc_int64,      ...
                                     D(nlocs).did);
      string = 'elapsed observation time since reference';
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'long_name', string);
      epoch  = datenum(num2str(S.datetime_ref),'yyyymmddHH');
      string = ['seconds since ' datestr(epoch, 'yyyy-mm-ddTHH:MM:SSZ')];
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'units', string);
    case 'date_time'
      G(Meta).vid(i) = netcdf.defVar(G(Meta).gid, Vname, nc_string,     ...
                                     D(nlocs).did);
      string = 'ISO 8601 UTC date and time string';
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'long_name', string);
      string = blanks(0);
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), '_FillValue', string,  ...
                    nc_string);
    case 'depth'
      G(Meta).vid(i) = netcdf.defVar(G(Meta).gid, Vname, nc_real,       ...
                                     D(nlocs).did);
      string = 'observation depth below sea level';
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'long_name', string);
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'units', 'meter');
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'negative', 'downwards');
    case 'latitude'
      G(Meta).vid(i) = netcdf.defVar(G(Meta).gid, Vname, nc_real,       ...
                                     D(nlocs).did);
      string = 'observation latitude';
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'long_name', string);
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'units', 'degrees_north');
    case 'longitude'
      G(Meta).vid(i) = netcdf.defVar(G(Meta).gid, Vname, nc_real,       ...
                                     D(nlocs).did);
      string = 'observation longitude';
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'long_name', string);
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'units', 'degrees_east');
    case 'provenance'
      G(Meta).vid(i) = netcdf.defVar(G(Meta).gid, Vname, nc_int,        ...
                                     D(nlocs).did);
      string = 'observation origin identifier';
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'long_name', string);
      if (isfield(S, 'flag_values'))
        netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'flag_values',       ...
	              S.flag_values);
        netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'flag_meanings',     ...
	              S.flag_meanings);
      end	
    case 'sequenceNumber'
      G(Meta).vid(i) = netcdf.defVar(G(Meta).gid, Vname, nc_int,        ...
                                     D(nlocs).did);
      string = 'observations sequence number';
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'long_name', string);
    case 'variables_name'
      G(Meta).vid(i) = netcdf.defVar(G(Meta).gid, Vname, nc_string,     ...
                                     D(nvars).did);
      string = 'observation UFO/IODA standard name';
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'long_name', string);
      string = blanks(0);
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), '_FillValue', string,  ...
                    nc_string);
  end
end

%--------------------------------------------------------------------------
% Define variables in the 'ObsError' Group.
%--------------------------------------------------------------------------

for i = 1:S.nvars
  Vname = char(S.ncvname(i));
  G(ObsE).vid(i) = netcdf.defVar(G(ObsE).gid, Vname, nc_real,           ...
                                 D(nlocs).did);
  string = 'observation error standard deviation';
  netcdf.putAtt(G(ObsE).gid, G(ObsE).vid(i), 'long_name', string);
  string = char(S.units(i));
  if (length(string) > 0)
    netcdf.putAtt(G(ObsE).gid, G(ObsE).vid(i), 'units', string);
  end
  if (do_depth)
    string = 'longitude, latitude, depth, dateTime';
  else
    string = 'longitude, latitude, dateTime';
  end
  netcdf.putAtt(G(ObsE).gid, G(ObsE).vid(i), 'coordinates', string);
end

%--------------------------------------------------------------------------
% Define variables in the 'ObsValue' Group.
%--------------------------------------------------------------------------

for i = 1:S.nvars
  Vname = char(S.ncvname(i));
  G(ObsV).vid(i) = netcdf.defVar(G(ObsV).gid, Vname, nc_real,           ...
                                 D(nlocs).did);
  switch Vname
    case 'absolute_dynamic_topography'
      string = 'observation value, ADT = MDT + SLA';
    otherwise
      string = 'observation value';
  end
  netcdf.putAtt(G(ObsV).gid, G(ObsV).vid(i), 'long_name', string);
  string = char(S.units(i));
  if (length(string) > 0)
    netcdf.putAtt(G(ObsV).gid, G(ObsV).vid(i), 'units', string);
  end
  if (do_depth)
    string = 'longitude, latitude, depth, dateTime';
  else
    string = 'longitude, latitude, dateTime';
  end
  netcdf.putAtt(G(ObsV).gid, G(ObsV).vid(i), 'coordinates', string);
end

%--------------------------------------------------------------------------
% Define variables in the 'PreQC' Group.
%--------------------------------------------------------------------------

for i = 1:S.nvars
  Vname = char(S.ncvname(i));
  G(PreQ).vid(i) = netcdf.defVar(G(PreQ).gid, Vname, nc_int,            ...
                                 D(nlocs).did);
  string = 'observation preset quality control filter identifier';
  netcdf.putAtt(G(PreQ).gid, G(PreQ).vid(i), 'long_name', string);
  if (do_depth)
    string = 'longitude, latitude, depth, dateTime';
  else
    string = 'longitude, latitude, dateTime';
  end
  netcdf.putAtt(G(PreQ).gid, G(PreQ).vid(i), 'coordinates', string);
end

%--------------------------------------------------------------------------
% End NetCDF4 file definition.
%--------------------------------------------------------------------------

netcdf.endDef(ncid);

%--------------------------------------------------------------------------
% Write some variables.
%--------------------------------------------------------------------------

% Ungroupled variables.

for i = 1:length(D)
  size = D(i).size;
  values = zeros([1 size]);
  netcdf.putVar(ncid, D(i).vid, values);
end

%--------------------------------------------------------------------------
% Close NetCDF4 file.
%--------------------------------------------------------------------------

netcdf.close(ncid);

return