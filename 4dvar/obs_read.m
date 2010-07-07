function [S]=obs_read(ncfile);

%
% OBS_READ:  Reads ROMS 4D-Var observation NetCDF file
%
% [S]=obs_read(ncfile)
%
% This function reads ROMS 4D-Var observation NetCDF file and stores all
% the variables in structure array, S.
%
% On Input:
%
%    ncfile  4D-Var Observations NetCDF file name (string)
%
% On Output:
%
%    S       Observations data (structure array):
%
%              S.ncfile         NetCDF file name (string)
%              S.Nsurvey        number of observations surveys
%              S.Nstate         number of state variables
%              S.Ndatum         total number of observations
%              S.spherical      spherical grid switch
%              S.Nobs           number of observations per survey
%              S.survey_time    time for each survey time
%              S.variance       global variance per state variable
%              S.type           state variable associated with observation
%              S.time           time for each observation
%              S.depth          depth of observation
%              S.Xgrid          observation fractional x-grid location
%              S.Ygrid          observation fractional y-grid location
%              S.Zgrid          observation fractional z-grid location
%              S.error          observation error
%              S.value          observation value
%
%            The following optional variables will be read if available:
%
%              S.provenance     observation origin
%              S.lon            observation longitude
%              S.lat            observation latitude
% 
%            The following variable attributes will be read if available:
%
%              S.state_flag_values     obs_type 'flag_values' attribute
%              S.state_flag_meanings   obs_type 'flag_meanings attribute
%              S.origin_flag_values    obs_provenance 'flag_values' attribute
%              S.origin_flag_meanings  obs_provenance 'flag_meaning' attribute
%
%            The following global attributes will be read if available:
%
%              S.global_variables      'state_variables' global attribute
%              S.global_provenance     'obs_provenance' global attribute
%              S.global_sources        'obs_sources' global attribute
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Initialize output structure.

S=[];

%----------------------------------------------------------------------------
%  Read in all available variables.
%----------------------------------------------------------------------------

[Vnames,nvars]=nc_vname(ncfile);

S.ncfile=ncfile;
S.Nsurvey=0;
S.Nstate=0;
S.Ndatum=0;

got.provenance=0;

for n=1:nvars,
  name=deblank(Vnames(n,:));
  switch name
    case 'spherical'
      S.spherical=nc_read(ncfile,'spherical');
    case 'Nobs'
      S.Nobs=nc_read(ncfile,'Nobs');
    case 'survey_time'
      S.survey_time=nc_read(ncfile,'survey_time');
      S.Nsurvey=length(S.survey_time);
    case 'obs_variance'
      S.variance=nc_read(ncfile,'obs_variance');
      S.Nstate=length(S.variance);
    case 'obs_type'
      S.type=nc_read(ncfile,'obs_type');
    case 'obs_provenance'
      S.provenance=nc_read(ncfile,'obs_provenance');
      got.provenance=1;
    case 'obs_time'
      S.time=nc_read(ncfile,'obs_time');
    case 'obs_lon'
      S.lon=nc_read(ncfile,'obs_lon');
    case 'obs_lat'
      S.lat=nc_read(ncfile,'obs_lat');
    case 'obs_depth'
      S.depth=nc_read(ncfile,'obs_depth');
    case 'obs_Xgrid'
      S.Xgrid=nc_read(ncfile,'obs_Xgrid');
    case 'obs_Ygrid'
      S.Ygrid=nc_read(ncfile,'obs_Ygrid');
    case 'obs_Zgrid'
      S.Zgrid=nc_read(ncfile,'obs_Zgrid');
    case 'obs_error'
      S.error=nc_read(ncfile,'obs_error');
    case 'obs_value'
      S.value=nc_read(ncfile,'obs_value');
      S.Ndatum=length(S.value);
  end,
end,

%----------------------------------------------------------------------------
%  Read in relevant variable and global attributes.
%----------------------------------------------------------------------------

%  Read in 'flag_values' and 'flag_meanings attributes for variable
%  'obs_type'.

Avalue=nc_getatt(ncfile,'flag_values','obs_type');
if (~isempty(Avalue)),
  S.state_flag_values=Avalue;
end,

Avalue=nc_getatt(ncfile,'flag_meanings','obs_type');
if (~isempty(Avalue)),
  S.state_flag_meanings=Avalue;
end,

%  Read in 'flag_values' and 'flag_meanings attributes for variable
%  'obs_provenace'.

if (got.provenance), 
  Avalue=nc_getatt(ncfile,'flag_values','obs_provenance');
  if (~isempty(Avalue)),
    S.origin_flag_values=Avalue;
  end,

  Avalue=nc_getatt(ncfile,'flag_meanings','obs_provenance');
  if (~isempty(Avalue)),
    S.origin_flag_meanings=Avalue;
  end,
end,

% Read in 'title' global attribute.

Avalue=nc_getatt(ncfile,'title');
if (~isempty(Avalue)),
  S.title=Avalue;
end,

% Read in 'state_variables' global attribute.

Avalue=nc_getatt(ncfile,'state_variables');
if (~isempty(Avalue)),
  S.global_variables=Avalue;
end,

% Read in 'obs_provenance' global attribute.

Avalue=nc_getatt(ncfile,'obs_provenance');
if (~isempty(Avalue)),
  S.global_provenance=Avalue;
end,

% Read in 'obs_sources' global attribute.

Avalue=nc_getatt(ncfile,'obs_sources');
if (~isempty(Avalue)),
  S.global_sources=Avalue;
end,


return
