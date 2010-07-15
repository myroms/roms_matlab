function [S]=obs_depth(ncfile,S,zflag);

%
% OBS_DEPTH:  Computes the observations fractional z-grid locations
%
% [S]=obs_read(ncfile)
%
% This function computes the observation fractional z-grid locations. The
% application vertical coordinates parameters are extracted from provided
% ROMS history file.
%
% On Input:
%
%    ncfile  ROMS history file name (string)
%    S       Observations strcuture data or NetCDF file name:
%    zflag   Flag to compute ROMS depths:
%
%              zflag = 0        use zero free-surface
%              zflag > 1        use zflag record for free-surface
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

%  Read in observation if S is a 4D-Var observation NetCDF file.

if (ischar(S)),
  S=obs_read(S);
end,

%  Check if 'provenace', 'lon', and 'lat' fields are available.

has.lonlat=false;
if (isfield(S,'lon') & isfield(S,'lat')),
  has.lonlat=true;
end,

has.provenance=false;
if (isfield(S,'provenance')),
  has.provenance=true;
end,

%  Get number of data surveys.

Nsurvey=length(S.survey_time);

%----------------------------------------------------------------------------
%  Compute depths of ROMS application.
%----------------------------------------------------------------------------

igrid = 5;                                  % depth of W-points
tindex = zflag;                             % use zflag record for zeta

switch ( zflag ),
  case 0
    Zw=depths(ncfile, ncfile, igrid, 0, 0);
  case 1
    Zw=depths(ncfile, ncfile, igrid, 0, tindex);
end,

[Lr Mr Nw]=size(Zw);
Nr = Nw -1;

%----------------------------------------------------------------------------
%  Compute observations fractional z-grid location. Use W-points depths
%  to include observations in the bottom half and top half of the grid
%  cells.
%----------------------------------------------------------------------------

ind = find(S.depth < 0);

if (~isempty(ind));
  for n=1:length(ind),
    iobs = ind(n);
    I = 1.0 + floor(S.Xgrid(iobs));
    J = 1.0 + floor(S.Ygrid(iobs));
    z = reshape(Zw(I,J,:),1,Nw);
    S.Zgrid = interp1(z, [0:1:Nr], S.depth(iobs));
  end,
end,

%  Since the state variables are at vertical RHO-points, assign bottom half
%  locations to k=1.

ind = find((0.0 < S.Zgrid) & (S.Zgrid <= 0.5));
if (~isempty(ind));
  S.Zgrid(ind) = 1.0;
end,

%  Since the state variables are at vertical RHO-points, assign top half
%  locations to k=Km-1.

ind = find((Nr < S.Zgrid) & (S.Zgrid <= Nr+0.5));
if (~isempty(ind));
  S.Zgrid(ind) = 1.0;
end,

%  Remove NaNs from the vertical interpolation.

ind = find(isnan(S.Zgrid));

if (~isempty(ind)),
  S.type(ind)  = [];
  S.time(ind)  = []; 
  S.depth(ind) = [];
  S.Xgrid(ind) = [];
  S.Ygrid(ind) = [];
  S.Zgrid(ind) = [];
  S.error(ind) = [];
  S.value(ind) = [];

  disp(' ');
  disp([' Number of removed outlier obserbations = ', ...
	num2str(length(ind))]);
  disp(' ');
  
  if (has.lonlat),
    S.lat(ind) = [];
    S.lon(ind) = [];
  end,

  if (has.provenance),
    S.provenance(ind) = [];
  end,

  for i=1:S.Nsurvey,
    ind=find(S.time == S.survey_time(i));
    S.Nobs(i)=length(ind);
  end,
end,

S.Ndatum = length(S.time);

return
