function [S]=ioda_read(ncfile)

%
% IODA_READ:  Reads IODA observation NetCDF4 file
%
% [S]=ioda_read(ncfile)
%
% This function reads IODA observation NetCDF4 file and stores all
% the variables in structure array, S.
%
% On Input:
%
%    ncfile  IODA Observations NetCDF4 file name (string)
%
% On Output:
%
%    S       IODA Observations data (structure array):
%
%              S.ncfile           NetCDF file name (string)
%              S.souce            Native 4D-Var source file (string)
%              S.nlocs            number of observations
%              S.nvars            number of variables
%              S.epoch            IODA time reference YYYYMMDDHH
%              S.datenum          epoch date number
%              S.dateTimeRef      IODA time reference string
%              S.variable_names   UFO/IODA standard name
%              S.dateTime         seconds since yyyy-mm-ddTHH:MM:SSZ
%              S.date_time        date and time ISO 8601 UTC string
%              S.depth            depth of observations
%              S.longitude        longitude of observations
%              S.latitude         latitude  of observations
%              S.provenance       observation origin identifier
%              S.values           observation values
%              S.units            observation units
%              S.errors           observation error
%              S.PreQC            observation quality control
%

% git $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

% Initialize.

S = struct('ncfile'           , [],                                     ...
           'source'           , [],                                     ...
           'nlocs'            , [],                                     ...
           'nvars'            , [],                                     ...
           'epoch'            , [],                                     ...
           'datenum'          , [],                                     ...
           'datetimeRef'      , [],                                     ...
           'iodaVarName'      , [],                                     ...
           'variables_name'   , [],                                     ...
           'shortname'        , [],                                     ...
           'units'            , [],                                     ...
           'dateTime'         , [],                                     ...
           'date_time'        , [],                                     ...
           'latitude'         , [],                                     ...
           'longitude'        , [],                                     ...
           'provenance'       , [],                                     ...
           'sequenceNumber'   , []);

% Inquire NetCDF4 file.

I = ncinfo(ncfile);

S.ncfile = ncfile;
S.nlocs  = I.Dimensions(strcmp({I.Dimensions.Name}, 'Location' )).Length;
S.nvars  = I.Dimensions(strcmp({I.Dimensions.Name}, 'nvars'    )).Length;

% Get IODA reference time YYYYMMDDHH global attribute.

S.epoch       = I.Attributes(strcmp({I.Attributes.Name}, 'date_time'  )).Value;
S.datetimeRef = I.Attributes(strcmp({I.Attributes.Name}, 'datetimeReference'  )).Value;
S.source      = I.Attributes(strcmp({I.Attributes.Name}, 'sourceFiles')).Value;

S.datenum     = datenum(num2str(S.epoch), 'yyyymmddHH');

% Read in 'MetaData' Group.

G = ncinfo(ncfile, 'MetaData');

S.longitude = double(ncread(ncfile, '/MetaData/longitude'));
S.latitude  = double(ncread(ncfile, '/MetaData/latitude'));

if (any(strcmp({G.Variables.Name}, 'depth')))
  S.depth = double(ncread(ncfile, '/MetaData/depth'));
end

if (any(strcmp({G.Variables.Name}, 'provenance')))
  S.provenance = double(ncread(ncfile, '/MetaData/provenance'));
else
  S = rmfield(S, 'provenance');
end

if (any(strcmp({G.Variables.Name}, 'sequenceNumber')))
  S.sequenceNumber = double(ncread(ncfile, '/MetaData/sequenceNumber'));
else
  S = rmfield(S, 'sequenceNumber');
end

S.variables_name = cellstr(ncread(ncfile, '/MetaData/variables_name'))';

S.dateTime = double(ncread(ncfile, '/MetaData/dateTime'));

if (any(strcmp({G.Variables.Name}, 'date_time')))
  S.date_time = ncread(ncfile, '/MetaData/date_time');
end

% Set IODA NetCDF variables.

for i = 1:S.nvars
  string = I.Groups(2).Variables(i).Name;
  S.iodaVarName{i} = string;
end

% Read in 'EffectiveError' Group.

if (any(strcmp({I.Groups.Name}, 'EffectiveError')))
  for i = 1:S.nvars
    Vname = strcat('/EffectiveError/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.EffectiveError{i} = field;
  end
end

% Read in 'EffectiveQC' Group.

if (any(strcmp({I.Groups.Name}, 'EffectiveQC')))
  for i = 1:S.nvars
    Vname = strcat('/EffectiveQC/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.EffectiveQC{i} = field;
  end
end

% Read in 'ObsBias' Group.

if (any(strcmp({I.Groups.Name}, 'ObsBias')))
  for i = 1:S.nvars
    Vname = strcat('/ObsError/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.ObsBias{i} = field;
  end
end

% Read in 'ObsError' Group.

if (any(strcmp({I.Groups.Name}, 'ObsError')))
  for i = 1:S.nvars
    Vname = strcat('/ObsError/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.ObsError{i} = field;
  end
end

% Read in 'ObsValue' Group.

if (any(strcmp({I.Groups.Name}, 'ObsValue')))
  for i = 1:S.nvars
    Vname = strcat('/ObsValue/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.ObsValue{i} = field;
    S.units{i}  = nc_getatt(ncfile, 'units', Vname);
  end
end

% Read in 'PreQC' Group.

if (any(strcmp({I.Groups.Name}, 'PreQC')))
  for i = 1:S.nvars
    Vname = strcat('/PreQC/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.PreQC{i} = field;
  end
end

% Read in 'hofx' Group: Model at observation locations, H(x).

if (any(strcmp({I.Groups.Name}, 'hofx')))
  for i = 1:S.nvars
    Vname = strcat('/hofx/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.hofx{i} = field;
  end
end

% Read in 'hofx0' Group: Initial H(x).

if (any(strcmp({I.Groups.Name}, 'hofx0')))
  for i = 1:S.nvars
    Vname = strcat('/hofx0/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.hofx0{i} = field;
  end
end

if (any(strcmp({I.Groups.Name}, 'hofx0_1')))
  for i = 1:S.nvars
    Vname = strcat('/hofx0_1/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.hofx0_1{i} = field;
  end
end

if (any(strcmp({I.Groups.Name}, 'hofx0_2')))
  for i = 1:S.nvars
    Vname = strcat('/hofx0_2/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.hofx0_2{i} = field;
  end
end

if (any(strcmp({I.Groups.Name}, 'hofx0_3')))
  for i = 1:S.nvars
    Vname = strcat('/hofx0_3/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.hofx0_3{i} = field;
  end
end

% Read in 'hofx1' Group: Final H(x).

if (any(strcmp({I.Groups.Name}, 'hofx1')))
  for i = 1:S.nvars
    Vname = strcat('/hofx1/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.hofx1{i} = field;
  end
end

if (any(strcmp({I.Groups.Name}, 'hofx1_1')))
  for i = 1:S.nvars
    Vname = strcat('/hofx0_1/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.hofx1_1{i} = field;
  end
end

if (any(strcmp({I.Groups.Name}, 'hofx1_2')))
  for i = 1:S.nvars
    Vname = strcat('/hofx0_2/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.hofx1_2{i} = field;
  end
end

if (any(strcmp({I.Groups.Name}, 'hofx1_3')))
  for i = 1:S.nvars
    Vname = strcat('/hofx0_3/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.hofx1_3{i} = field;
  end
end

% Read in 'oman' Group: Observation minus analysis.

if (any(strcmp({I.Groups.Name}, 'oman')))
  for i = 1:S.nvars
    Vname = strcat('/oman/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.oman{i} = field;
  end
end

% Read in 'ombg' Group: Observation minus background.

if (any(strcmp({I.Groups.Name}, 'ombg')))
  for i = 1:S.nvars
    Vname = strcat('/ombg/', S.iodaVarName{i});
    field = double(ncread(ncfile, Vname));
    S.ombg{i} = field;
  end
end

return
