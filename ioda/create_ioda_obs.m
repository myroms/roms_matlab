function create_ioda_obs(S, file)

%
% CREATE_IODA_OBS:   Created an IODA observation NetCDF4 file
%
% create_ioda_obs(S, file)
%
% This function creates an IODA observation file using specified structure,
% S, for data assimilation within JEDI.
%
% On Input:
%
%    S        Observations file creation parameters (struct):
%
%               S.ncfile            NetCDF file name (string)
%               S.ndatetime         number of datetime
%               S.nlocs             number of observations
%               S.nrecs             number of records
%               S.nstring           number of string
%               S.nvars             number of variables
%               S.variable_names(:) IODA variable standard name
%               S.units(:)          variables units
%               S.datetime_ref      IODA reference time YYYYMMDDHH
%               S.depth             depth of observation
%               S.provenance        observation provenance
%
%    file     Output observation file name (string, OPTIONAL). If provided,
%               it creates this file instead the one in S.ncfile.

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2022 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
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

% Set indices for local dimension structure vector element.

ndate = 1;           % ndatetime index
nlocs = 2;           % nlocs     index
nobs  = 3;           % nobs      index 
nrecs = 4;           % nrecs     index
nstr  = 5;           % nstring   index
nvars = 6;           % nvars     index

% Set indices for local group structure vector element.

Meta = 1;            % MetaData
ObsE = 2;            % ObsError
ObsV = 3;            % ObsValue
PreQ = 4;            % PreQC
VarM = 5;            % VarMetaData

% Set variable representation kind.

nc_char  = netcdf.getConstant('nc_char');       % string element
nc_int   = netcdf.getConstant('nc_int');        % integer
nc_int64 = netcdf.getConstant('nc_int64');      % 64-bit integer
nc_real  = netcdf.getConstant('nc_float');      % floating-point

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

D(ndate).name = 'ndatetime';    D(ndate).size = S.ndatetime;
D(nlocs).name = 'nlocs';        D(nlocs).size = S.nlocs;
D(nobs ).name = 'nobs';         D(nobs ).size = S.nobs;
D(nrecs).name = 'nrecs';        D(nrecs).size = S.nrecs;
D(nstr ).name = 'nstring';      D(nstr ).size = S.nstring;
D(nvars).name = 'nvars';        D(nvars).size = S.nvars;

%--------------------------------------------------------------------------
% Create NetCDF4 observation file.
%--------------------------------------------------------------------------

disp(' ');
disp(['*** Creating observations file:  ', ncfile]);

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode, netcdf.getConstant('NETCDF4'));

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
netcdf.putAtt(ncid, varid, '_ioda_layout_version', int32(2));
netcdf.putAtt(ncid, varid, 'nrecs', int32(S.nrecs));
netcdf.putAtt(ncid, varid, 'nvars', int32(S.nvars));
netcdf.putAtt(ncid, varid, 'nlocs', int32(S.nlocs));
netcdf.putAtt(ncid, varid, 'nobs',  int32(S.nobs));

netcdf.putAtt(ncid, varid, 'source_file',  S.source);

history=['Created from Matlab script: ', mfilename, ' on ', date_stamp];
netcdf.putAtt(ncid, varid, 'history', history);

%--------------------------------------------------------------------------
% Define NetCDF4 file variables: same as dimension names.
%--------------------------------------------------------------------------

for i = 1:length(D)
  name = char(D(i).name);
  D(i).vid = netcdf.defVar(ncid, name, nc_real, D(i).did);
  netcdf.putAtt(ncid, D(i).vid, 'suggested_chunck_dim', '100LL')
end

%--------------------------------------------------------------------------
% Define NetCDF4 files Groups: 'MetaData', 'ObsError', 'ObsValue', 'PreQc',
%                              and 'VarMetaData'.
%--------------------------------------------------------------------------

G(Meta).gid = netcdf.defGrp(ncid, 'MetaData');
G(ObsE).gid = netcdf.defGrp(ncid, 'ObsError');
G(ObsV).gid = netcdf.defGrp(ncid, 'ObsValue');
G(PreQ).gid = netcdf.defGrp(ncid, 'PreQC');
G(VarM).gid = netcdf.defGrp(ncid, 'VarMetaData');

%--------------------------------------------------------------------------
% Define variables in the 'MetaData' Group.
%--------------------------------------------------------------------------

MetaVars = {'dateTime', 'date_time', 'latitude', 'longitude',           ...
            'record_number'};
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
      G(Meta).vid(i) = netcdf.defVar(G(Meta).gid, Vname, nc_char,       ...
                                     [D(ndate).did D(nlocs).did]);
      string = 'ISO 8601 UTC date and time string';
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'long_name', string);
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
    case 'record_number'
      G(Meta).vid(i) = netcdf.defVar(G(Meta).gid, Vname, nc_int,        ...
                                     D(nlocs).did);
      string = 'group variable record number for multi-value observations';
      netcdf.putAtt(G(Meta).gid, G(Meta).vid(i), 'long_name', string);
  end
end

%--------------------------------------------------------------------------
% Define variables in the 'ObsError' Group.
%--------------------------------------------------------------------------

for i = 1:S.nvars
  Vname = char(S.variable_names(i));
  G(ObsE).vid(i) = netcdf.defVar(G(ObsE).gid, Vname, nc_real,           ...
                                 D(nlocs).did);
  string = 'observation error standard deviation';
  netcdf.putAtt(G(ObsE).gid, G(ObsE).vid(i), 'description', string);
  string = char(S.units(i));
  if (length(string) > 0)
    netcdf.putAtt(G(ObsE).gid, G(ObsE).vid(i), 'units', string);
  end
end

%--------------------------------------------------------------------------
% Define variables in the 'ObsValue' Group.
%--------------------------------------------------------------------------

for i = 1:S.nvars
  Vname = char(S.variable_names(i));
  G(ObsV).vid(i) = netcdf.defVar(G(ObsV).gid, Vname, nc_real,           ...
                                 D(nlocs).did);
  switch Vname
    case 'absolute_dynamic_topography'
      string = 'observation value, ADT = MDT + SLA';
    otherwise
      string = 'observation value';
  end
  netcdf.putAtt(G(ObsV).gid, G(ObsV).vid(i), 'description', string);
  string = char(S.units(i));
  if (length(string) > 0)
    netcdf.putAtt(G(ObsV).gid, G(ObsV).vid(i), 'units', string);
  end
end

%--------------------------------------------------------------------------
% Define variables in the 'PreQC' Group.
%--------------------------------------------------------------------------

for i = 1:S.nvars
  Vname = char(S.variable_names(i));
  G(PreQ).vid(i) = netcdf.defVar(G(PreQ).gid, Vname, nc_int,            ...
                                 D(nlocs).did);
  string = 'observation preset quality control filter identifier';
  netcdf.putAtt(G(PreQ).gid, G(PreQ).vid(i), 'description', string);
end

%--------------------------------------------------------------------------
% Define variables in the 'VarMetaData' Group.
%--------------------------------------------------------------------------

G(VarM).vid = netcdf.defVar(G(VarM).gid, 'variable_names', nc_char,     ...
                            [D(nstr).did D(nvars).did]);
string = 'observation UFO/IODA standard name';
netcdf.putAtt(G(VarM).gid, G(VarM).vid(i), 'description', string);

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

% VarMetaData values.

for i = 1:S.nvars
  Vname = char(S.variable_names(i));
  start = [0 i-1];
  count = [length(Vname) 1];
  netcdf.putVar(G(VarM).gid, G(VarM).vid, start, count, Vname);
end

%--------------------------------------------------------------------------
% Close NetCDF4 file.
%--------------------------------------------------------------------------

netcdf.close(ncid);

return