function [S]=ioda_read(ncfile)

%
% IODA_READ:  Reads IODA observation NetCDF-4 file
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
%    S       IODA Observations data (struct and cell arrays):
%
%              S.ncfile           NetCDF file name (string)
%              S.souce            Native 4D-Var source file (string)
%              S.nlocs            number of observations
%              S.nvars            number of variables
%              S.epoch            IODA time reference YYYYMMDDHH
%              S.datenum          epoch date number
%              S.dateTimeRef      IODA time reference string
%
%              Group MetaData:
%
%              S.dateTime         seconds since yyyy-mm-ddTHH:MM:SSZ
%              S.depth            depth of observations
%              S.longitude        longitude of observations
%              S.latitude         latitude  of observations
%              
%              S.provenance       observation origin identifier
%              S.sequenceNumber   observation sequence number
%              S.stateID          ROMS state variable index
%              S.surveyIndex      observation survey time indices
%              S.surveyTime       observation survey time
%              S.variable_names   UFO/IODA standard name
%              S.x_grid           observation fractional x-grid location
%              S.y_grid           observation fractional y-grid location
%              S.z_grid           observation fractional z-grid location
%  
%              Other Groups, cell arrays(nvars):
%
%              S.ObsValue         observation values
%              S.units            observation units
%              S.ObsError         observation error
%              S.PreQC            observation quality control
%

% git $Id$
%=========================================================================%
%  Copyright (c) 2002-2026 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

% Initialize.

S = struct('ncfile'           , [],                                     ...
           'roms_grid'        , [],                                     ...
           'source'           , [],                                     ...
           'nlocs'            , [],                                     ...
           'nvars'            , [],                                     ...
           'nsurvey'          , [],                                     ...
           'TimeIODA'         , [],                                     ...
           'DateIODA'         , [],                                     ...
           'datenum'          , [],                                     ...
           'datetime_ref'     , [],                                     ...
           'ncvname'          , [],                                     ...
           'variables_name'   , [],                                     ...
           'units'            , [],                                     ...
           'dateTime'         , [],                                     ...
           'timeAvgBegin'     , [],                                     ...
           'timeAvgEnd'       , [],                                     ...
           'depth'            , [],                                     ...
           'latitude'         , [],                                     ...
           'longitude'        , [],                                     ...
           'provenance'       , [],                                     ...
           'sequenceNumber'   , [],                                     ...
           'spatialAverage'   , [],                                     ...
           'stateID'          , [],                                     ...
           'surveyIndex'      , [],                                     ...
	   'surveyTime'       , [],                                     ...
	   'x_grid'           , [],                                     ...
	   'y_grid'           , [],                                     ...
	   'z_grid'           , []);

% Inquire NetCDF4 file.

I = ncinfo(ncfile);

S.ncfile  = ncfile;
S.nlocs   = I.Dimensions(strcmp({I.Dimensions.Name}, 'Location' )).Length;
S.nvars   = I.Dimensions(strcmp({I.Dimensions.Name}, 'nvars'    )).Length;
S.nsurvey = I.Dimensions(strcmp({I.Dimensions.Name}, 'survey'    )).Length;

% Get IODA reference time YYYYMMDDHH global attribute.

S.TimeIODA = I.Attributes(strcmp({I.Attributes.Name}, 'date_time')).Value;
S.DateIODA = I.Attributes(strcmp({I.Attributes.Name}, 'datetimeReference')).Value;
S.source   = I.Attributes(strcmp({I.Attributes.Name}, 'sourceFiles')).Value;

S.datenum  = datenum(num2str(S.TimeIODA), 'yyyymmddHH');
S.datetime_ref = S.TimeIODA;

% Get ROMS application grid dimensions.

if (any(strcmp({I.Attributes.Name}, 'roms_grid')))
  S.roms_grid = I.Attributes(strcmp({I.Attributes.Name}, 'roms_grid')).Value;
end

% Read in 'MetaData' Group.

G = ncinfo(ncfile, 'MetaData');

if (any(strcmp({G.Variables.Name}, 'dateTime')))
  S.dateTime = double(ncread(ncfile, '/MetaData/dateTime'));
else
  S = rmfield(S, 'dateTime');
end

if (any(strcmp({G.Variables.Name}, 'dateTimeAverageBegin')))
  S.timeAvgBegin = double(ncread(ncfile, '/MetaData/dateTimeAverageBegin'));
else
  S = rmfield(S, 'timeAvgBegin');
end

if (any(strcmp({G.Variables.Name}, 'dateTimeAverageEnd')))
  S.timeAvgEnd = double(ncread(ncfile, '/MetaData/dateTimeAverageEnd'));
else
  S = rmfield(S, 'timeAvgEnd');
end

if (any(strcmp({G.Variables.Name}, 'depth')))
  S.depth = double(ncread(ncfile, '/MetaData/depth'));
else
  S = rmfield(S, 'depth');
end

if (any(strcmp({G.Variables.Name}, 'longitude')))
  S.longitude = double(ncread(ncfile, '/MetaData/longitude'));
else
  S = rmfield(S, 'longitude');
end

if (any(strcmp({G.Variables.Name}, 'latitude')))
  S.latitude  = double(ncread(ncfile, '/MetaData/latitude'));
else
  S = rmfield(S, 'latitude');
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

if (any(strcmp({G.Variables.Name}, 'stateID')))
  S.stateID = double(ncread(ncfile, '/MetaData/stateID'));
else
  S = rmfield(S, 'stateID');
end

if (any(strcmp({G.Variables.Name}, 'spatialAverage')))
  S.areaAvgRadius = double(ncread(ncfile, '/MetaData/stateID'));
else
  S = rmfield(S, 'spatialAverage');
end

if (any(strcmp({G.Variables.Name}, 'surveyIndex')))
  S.surveyIndex = double(ncread(ncfile, '/MetaData/surveyIndex'));
else
  S = rmfield(S, 'surveyIndex');
end

if (any(strcmp({G.Variables.Name}, 'surveyTime')))
  S.surveyTime = double(ncread(ncfile, '/MetaData/surveyTime'));
else
  S = rmfield(S, 'surveyTime');
end

if (any(strcmp({G.Variables.Name}, 'x_grid')))
  S.x_grid = double(ncread(ncfile, '/MetaData/x_grid'));
else
  S = rmfield(S, 'x_grid');
end

if (any(strcmp({G.Variables.Name}, 'y_grid')))
  S.y_grid = double(ncread(ncfile, '/MetaData/y_grid'));
else
  S = rmfield(S, 'y_grid');
end

if (any(strcmp({G.Variables.Name}, 'z_grid')))
  S.z_grid = double(ncread(ncfile, '/MetaData/z_grid'));
else
  S = rmfield(S, 'z_grid');
end

S.variables_name = cellstr(ncread(ncfile, '/MetaData/variables_name'))';


% Set IODA NetCDF variables.

for i = 1:S.nvars
  string = I.Groups(2).Variables(i).Name;
  S.ncvname{i} = string;
end

% Read in 'EffectiveError' Group.

if (any(strcmp({I.Groups.Name}, 'EffectiveError')))
  for i = 1:S.nvars
    Vname = strcat('/EffectiveError/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.EffectiveError{i} = field;
  end
end

% Read in 'EffectiveQC' Group.

if (any(strcmp({I.Groups.Name}, 'EffectiveQC')))
  for i = 1:S.nvars
    Vname = strcat('/EffectiveQC/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.EffectiveQC{i} = field;
  end
end

% Read in 'ObsBias' Group.

if (any(strcmp({I.Groups.Name}, 'ObsBias')))
  for i = 1:S.nvars
    Vname = strcat('/ObsError/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.ObsBias{i} = field;
  end
end

% Read in 'ObsError' Group.

if (any(strcmp({I.Groups.Name}, 'ObsError')))
  for i = 1:S.nvars
    Vname = strcat('/ObsError/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.ObsError{i} = field;
  end
end

% Read in 'ObsValue' Group.

if (any(strcmp({I.Groups.Name}, 'ObsValue')))
  for i = 1:S.nvars
    Vname = strcat('/ObsValue/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.ObsValue{i} = field;
    S.units{i}  = nc_getatt(ncfile, 'units', Vname);
  end
end

% Read in 'PreQC' Group.

if (any(strcmp({I.Groups.Name}, 'PreQC')))
  for i = 1:S.nvars
    Vname = strcat('/PreQC/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.PreQC{i} = field;
  end
end

% Read in 'hofx' Group: Model at observation locations, H(x).

if (any(strcmp({I.Groups.Name}, 'hofx')))
  for i = 1:S.nvars
    Vname = strcat('/hofx/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.hofx{i} = field;
  end
end

% Read in 'hofxInitial' Group: Initial H(x), native ROMS.

if (any(strcmp({I.Groups.Name}, 'hofxInitial')))
  for i = 1:S.nvars
    Vname = strcat('/hofxInitial/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.hofxInitial{i} = field;
  end
end

% Read in 'hofxFinal' Group: Final H(x), native ROMS.

if (any(strcmp({I.Groups.Name}, 'hofxFinal')))
  for i = 1:S.nvars
    Vname = strcat('/hofxFinal/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.hofxFinal{i} = field;
  end
end

% Read in 'hofx0' Group: Initial H(x).

if (any(strcmp({I.Groups.Name}, 'hofx0')))
  for i = 1:S.nvars
    Vname = strcat('/hofx0/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.hofx0{i} = field;
  end
end

if (any(strcmp({I.Groups.Name}, 'hofx0_1')))
  for i = 1:S.nvars
    Vname = strcat('/hofx0_1/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.hofx0_1{i} = field;
  end
end

if (any(strcmp({I.Groups.Name}, 'hofx0_2')))
  for i = 1:S.nvars
    Vname = strcat('/hofx0_2/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.hofx0_2{i} = field;
  end
end

if (any(strcmp({I.Groups.Name}, 'hofx0_3')))
  for i = 1:S.nvars
    Vname = strcat('/hofx0_3/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.hofx0_3{i} = field;
  end
end

% Read in 'hofx1' Group: Final H(x).

if (any(strcmp({I.Groups.Name}, 'hofx1')))
  for i = 1:S.nvars
    Vname = strcat('/hofx1/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.hofx1{i} = field;
  end
end

if (any(strcmp({I.Groups.Name}, 'hofx1_1')))
  for i = 1:S.nvars
    Vname = strcat('/hofx0_1/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.hofx1_1{i} = field;
  end
end

if (any(strcmp({I.Groups.Name}, 'hofx1_2')))
  for i = 1:S.nvars
    Vname = strcat('/hofx0_2/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.hofx1_2{i} = field;
  end
end

if (any(strcmp({I.Groups.Name}, 'hofx1_3')))
  for i = 1:S.nvars
    Vname = strcat('/hofx0_3/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.hofx1_3{i} = field;
  end
end

% Read in 'oman' or 'Residual' Group: Observation minus analysis.

if (any(strcmp({I.Groups.Name}, 'oman')))
  for i = 1:S.nvars
    Vname = strcat('/oman/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.oman{i} = field;
  end
end

if (any(strcmp({I.Groups.Name}, 'Residual')))
  for i = 1:S.nvars
    Vname = strcat('/Residual/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.Residual{i} = field;
  end
end

% Read in 'ombg' or 'Innovation' Group: Observation minus background.

if (any(strcmp({I.Groups.Name}, 'ombg')))
  for i = 1:S.nvars
    Vname = strcat('/ombg/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.ombg{i} = field;
  end
end

if (any(strcmp({I.Groups.Name}, 'Innovation')))
  for i = 1:S.nvars
    Vname = strcat('/Innovation/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.Innovation{i} = field;
  end
end

% Read in 'Increment' Group: Analysis minus background.

if (any(strcmp({I.Groups.Name}, 'Increment')))
  for i = 1:S.nvars
    Vname = strcat('/Increment/', S.ncvname{i});
    field = double(ncread(ncfile, Vname));
    S.Increment{i} = field;
  end
end

return
