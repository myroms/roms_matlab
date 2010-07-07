function [S]=obs_merge(files,dt);

%
% OBS_MERGE:  Merges data from several 4D-Var observation NetCDF files
%
% [S]=obs_merge(obs_file,dt)
%
% Given a cell array of 4D-Var NetCDF observation files, this function
% combines them into a single observation data structure. The DT argument
% is the minimum amount of time that surveys must be separated by (MUST
% BE GREATER THAN DT OF YOUR ROMS CONFIGURATION). Surveys that are closer
% than DT together will be combined.
%
% On Input:
%
%    files   Observations NetCDF file names (string cell array)
%    dt  
%
% On Output:
%
%    S       Observations data (structure array):
%
%              S.ncfile         NetCDF file name (string)
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
% NOTE: On exit, a new structure, S, is returned which can be written to
%       disk using 'c_observations' and 'obs_write'.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license           Brian Powell            %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Initialize internal parameters.

check_time=true;
if (nargin < 2),
  check_time=false;
end,

Nfiles=length(files);

%----------------------------------------------------------------------------
%  Process each observation file and build observation structure.
%----------------------------------------------------------------------------

has.lonlat     = false;
has.provenance = false;

%  Initialize structure with first file.

[S]=obs_read(char(files(1)));

if (isfield(S,'lon') & isfield(S,'lat')),
  has.lonlat=true;
end,

if (isfield(S,'provenance')),
  has.provenance=true;
end,

%  Proccess remaining files.

for m=2:Nfiles,

  [N]=obs_read(char(files(m)));

  if (has.provenance & isfield(N,'provenance')),
    has.provenance=true;
  else,
    has.provenance=false;
  end,

  if (has.lonlat & (isfield(N,'lon') & isfield(N,'lat'))),
    has.lonlat=true;
  else,
    has.lonlat=false;
  end,

%  Add new data to S structure.

  S.spherical    = ( S.spherical +  N.spherical   ) / 2;
  S.Nobs         = [ S.Nobs;        N.Nobs        ];
  S.survey_time  = [ S.survey_time; N.survey_time ];
  S.variance     = ( S.variance  +  N.variance    ) / 2.0;
  S.type         = [ S.type;        N.type        ];
  S.time         = [ S.time;        N.time        ];
  S.depth        = [ S.depth;       N.depth       ];
  S.Xgrid        = [ S.Xgrid;       N.Xgrid       ];
  S.Ygrid        = [ S.Ygrid;       N.Ygrid       ];
  S.Zgrid        = [ S.Zgrid;       N.Zgrid       ];
  S.error        = [ S.error;       N.error       ];
  S.value        = [ S.value;       N.value       ];

  if (has.provenance),
    S.provenance = [ S.provenance;  N.provenance  ];
  end,

  if (has.lonlat),
    S.lon        = [ S.lon;         N.lon         ];
    S.lat        = [ S.lat;         N.lat         ];
  end,

%  Order according to ascending 'survey_time'.

  [l,ind]=sort(S.survey_time);

  S.survey_time  = S.survey_time(ind);
  S.Nobs         = S.Nobs(ind);

  S.Nsurvey      = length(S.Nobs);

%  Order according to ascending 'obs_time'.

  [l,ind]=sort(S.time);

  S.type         = S.type(ind);
  S.time         = S.time(ind);
  S.depth        = S.depth(ind);
  S.Xgrid        = S.Xgrid(ind);
  S.Ygrid        = S.Ygrid(ind);
  S.Zgrid        = S.Zgrid(ind);
  S.error        = S.error(ind);
  S.value        = S.value(ind);

  if (has.provenance),
    S.provenance = S.provenance(ind);
  end,

  if (has.lonlat),
    S.lon        = S.lon(ind);
    S.lat        = S.lat(ind);
  end,

end,

%  Loop over 'survey_times' until they are well-spaced.

if (check_time),
  
  while (true),
    ntime = S.survey_time;
    sdt = diff(S.survey_time);
    l = find(sdt ~= 0 & sdt < dt);
    if (isempty( l )),
      break;
    end,
    otime = S.time;
    for i=1:length(l),
      lo = find(S.time == S.survey_time(l(i)   ) | ...
		S.time == S.survey_time(l(i)+1));
      ntime(l(i):l(i)+1) = nanmedian(S.time(lo));
      otime(lo) = ntime(l(i));
    end,
    S.time = otime;
    S.survey_time = ntime;
  end,

end,

%  Combine all surveys with similar times.

[l,i] = unique(S.survey_time);

if (length(l) ~= length(S.survey_time)),
  S.spherical = S.spherical(i);
  survey_time = S.survey_time;
  S.survey_time = l;
  Nobs = S.Nobs;
  S.Nobs = [];
  for i=1:length(l)
    m = find(survey_time == l(i));
    S.Nobs(i,1) = sum(Nobs(m));
  end,
end,

S.Nsurvey = length(S.survey_time);

%  Sort the observations by the 'obs_type' for each survey.

for i=1:S.Nsurvey,

  l = find(S.time == S.survey_time(i));
  [m,i]=sort(S.type(l));

  S.type(l)         = S.type(l(i));
  S.time(l)         = S.time(l(i)); 
  S.depth(l)        = S.depth(l(i));
  S.Xgrid(l)        = S.Xgrid(l(i));
  S.Ygrid(l)        = S.Ygrid(l(i));
  S.Zgrid(l)        = S.Zgrid(l(i));
  S.error(l)        = S.error(l(i));
  S.value(l)        = S.value(l(i));

  if (has.lonlat),
    S.lat(l)        = S.lat(l(i));
    S.lon(l)        = S.lon(l(i));
  end,

  if (has.provenance),
    S.provenance(l) = S.provenance(l(i));
  end,

end,

% Lastly, remove any NaN's that may have slipped through.

l=find(isnan(S.value) | isnan(S.time)   | ...
       isnan(S.Xgrid) | isnan(S.Ygrid)  | ...
       isnan(S.Zgrid) | isnan(S.depth));

if (~isempty(l)),
  S.type(l)         = [];
  S.time(l)         = []; 
  S.depth(l)        = [];
  S.Xgrid(l)        = [];
  S.Ygrid(l)        = [];
  S.Zgrid(l)        = [];
  S.error(l)        = [];
  S.value(l)        = [];

  if (has.lonlat),
    S.lat(l)        = [];
    S.lon(l)        = [];
  end,

  if (has.provenance),
    S.provenance(l) = [];
  end,

  for i=1:S.Nsurvey,
    ind=find(S.time == S.survey_time(i));
    S.Nobs(i)=length(ind);
  end,
end,

S.Ndatum = length(S.time);
S.files  = files;

return
