function [Sout]=super_obs(Sinp);

%
% SUPER_OBS:  Creates super observations when necessary
%
% [Sout]=super_obs(Sinp)
%
% This function checks the provided observation data and creates
% super observations if there is more than one datum of the same
% observation type per grid cell. At input, Sinp is either a
% 4D-Var observation NetCDF file or data structure.
%
% On Input:
%
%    Sinp    Observations data structure or NetCDF file name
%
%
% On Output:
%
%    Sout    Binned observations data (structure array):
%
%              Sout.ncfile       NetCDF file name (string)
%              Sout.Ndatum       total number of observations
%              Sout.spherical    spherical grid switch
%              Sout.Nobs         number of observations per survey
%              Sout.survey_time  time for each survey time
%              Sout.variance     global variance per state variable
%              Sout.type         state variable associated with observation
%              Sout.time         time for each observation
%              Sout.depth        depth of observation
%              Sout.Xgrid        observation fractional x-grid location
%              Sout.Ygrid        observation fractional y-grid location
%              Sout.Zgrid        observation fractional z-grid location
%              Sout.error        observation error
%              Sout.value        observation value
%
%            The following optional variables will be process if available:
%
%              Sout.provenance   observation origin
%              Sout.lon          observation longitude
%              Sout.lat          observation latitude
% 
% The nice and efficient sparse matrix approach is based on Bartolome
% Garau gliders data binning scripts.
%
  
  
% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Read observations if 'Sinp' is a 4D-Var observation NetCDF file.

if (ischar(Sinp)),
  Sinp=obs_read(Sinp);
end,

Lm = Sinp.grid_Lm_Mm_N(1);
Mm = Sinp.grid_Lm_Mm_N(2);
N  = Sinp.grid_Lm_Mm_N(3);

%  Check if 'provenace', 'lon', and 'lat' fields are available.

has.lonlat=false;
if (isfield(Sinp,'lon') & isfield(Sinp,'lat')),
  has.lonlat=true;
end,

has.provenance=false;
if (isfield(Sinp,'provenance')),
  has.provenance=true;
end,

%  Get number of data surveys.

Nsurvey=length(Sinp.survey_time);

%  Find observations associated with the same state variable.

state_vars=unique(Sinp.type);
Nstate=length(state_vars);

%  Initialize output structure.

Sout.ncfile       = Sinp.ncfile;
Sout.Nsurvey      = Nsurvey;
Sout.Nstate       = length(Sinp.variance);
Sout.Ndatum       = [];
Sout.spherical    = Sinp.spherical;
Sout.Nobs         = zeros(size(Sinp.Nobs));
Sout.survey_time  = Sinp.survey_time;
Sout.variance     = Sinp.variance;
Sout.type         = [];
Sout.time         = [];
Sout.depth        = [];
Sout.Xgrid        = [];
Sout.Ygrid        = [];
Sout.Zgrid        = [];
Sout.error        = [];
Sout.value        = [];

if (has.provenance),
  Sout.provenance = [];
end,

if (has.lonlat),
  Sout.lon        = [];
  Sout.lat        = [];
end,

%----------------------------------------------------------------------------
%  Compute super observations when needed.
%----------------------------------------------------------------------------

for m=1:Nsurvey,

%  Extract observations with the same survey time.

  ind_t=find(Sinp.time == Sinp.survey_time(m));

  T.type  = Sinp.type (ind_t);
  T.time  = Sinp.time (ind_t);
  T.depth = Sinp.depth(ind_t);
  T.Xgrid = Sinp.Xgrid(ind_t);
  T.Ygrid = Sinp.Ygrid(ind_t);
  T.Zgrid = Sinp.Zgrid(ind_t);
  T.error = Sinp.error(ind_t);
  T.value = Sinp.value(ind_t);

  if (has.provenance),
    T.provenance = Sinp.provenance(ind_t);
  end,

  if (has.lonlat);
    T.lon = Sinp.lon(ind_t);
    T.lat = Sinp.lat(ind_t);
  end,

  for n=1:Nstate,

%   Now extract observations associated with the same state variable.

    ind_v=find(T.type == state_vars(n));

    if isempty(ind_v),
      continue
    end,

    V.type  = T.type (ind_v);
    V.time  = T.time (ind_v);
    V.depth = T.depth(ind_v);
    V.Xgrid = T.Xgrid(ind_v);
    V.Ygrid = T.Ygrid(ind_v);
    V.Zgrid = T.Zgrid(ind_v);
    V.error = T.error(ind_v);
    V.value = T.value(ind_v);

    if (has.provenance),
      V.provenance = T.provenance(ind_v);
    end,

    if (has.lonlat);
      V.lon = T.lon(ind_v);
      V.lat = T.lat(ind_v);
    end,

%  Set binning parameters.

    Xmin = min(V.Xgrid);
    Xmax = max(V.Xgrid);

    Ymin = min(V.Ygrid);
    Ymax = max(V.Ygrid);

    Zmin = min(V.Zgrid);
    Zmax = max(V.Zgrid);

    dx = 1.0;
    dy = 1.0;
    dz = 1.0;
    
    minObs= 1;

%  Compute the index in each dimension of the grid cell in
%  which the observation is located.

    Xbin = 1.0 + floor((V.Xgrid - Xmin) ./ dx);
    Ybin = 1.0 + floor((V.Ygrid - Ymin) ./ dy);
    Zbin = 1.0 + floor((V.Zgrid - Zmin) ./ dz);

%  Similarly, compute the maximum averaging grid size.

    Xsize = 1.0 + floor((Xmax - Xmin) ./ dx);
    Ysize = 1.0 + floor((Ymax - Ymin) ./ dy);
    Zsize = 1.0 + floor((Zmax - Zmin) ./ dz);
    
    matsize = [Ysize, Xsize, Zsize];

%  Combine the indexes in each dimension into one index. It
%  is like stacking all the matrix in one column vector.

    varInd    = sub2ind(matsize, Ybin, Xbin, Zbin);
    onesCol   = ones(size(varInd));
    maxVarInd = max(varInd);

%  Process the grid information using sparse matrix:
%
%  - A sparse matrix store any non-zero values with their index in
%    the matrix
%  - It sum up all the elements that fall in the same matrix position
%
%  This is an efficient technique because:
%
%  - Many cells of the matrix will not have observations in it
%    (many zeros)
%  - Many observations will fall in the same grid cell
%    (they will be automatically summed up)
%
%  The sparse matrix will represent a column vector
%    
%  How many observations fall into each grid cell, summing all the
%  values that have the same index in the matrix

    binCounter = sparse(varInd, onesCol, onesCol, maxVarInd, 1);

%  Accumulate variable values at each cell.

    Xsum = sparse(varInd, onesCol, V.Xgrid, maxVarInd, 1);
    Ysum = sparse(varInd, onesCol, V.Ygrid, maxVarInd, 1);
    Zsum = sparse(varInd, onesCol, V.Zgrid, maxVarInd, 1);

%  Compute the mean of the observations:
%
%    Mean:     Xsum / numObs
%
%    (Be careful with variance if numObs == 1)

    [nonzeroIdx, dummyCol, dummyVals] = find(binCounter >= minObs);

    YmeanMat = Ysum(nonzeroIdx) ./ binCounter(nonzeroIdx);
    XmeanMat = Xsum(nonzeroIdx) ./ binCounter(nonzeroIdx);
    ZmeanMat = Zsum(nonzeroIdx) ./ binCounter(nonzeroIdx);

    clear dummy*;

%  Recover the full matrix after averaging.

    Sout.Xgrid = [Sout.Xgrid, transpose(full(XmeanMat))];
    Sout.Ygrid = [Sout.Ygrid, transpose(full(YmeanMat))];
    Sout.Zgrid = [Sout.Zgrid, transpose(full(YmeanMat))];

    Nsuper=length(full(XmeanMat));
    
%  Process the rest of the observation variables in similar fashion.

    if (has.lonlat);
      Fsum       = sparse(varInd, onesCol, V.lon, maxVarInd, 1);
      FmeanMat   = Fsum(nonzeroIdx) ./ binCounter(nonzeroIdx);
      Sout.lon   = [Sout.lon, transpose(full(FmeanMat))];

      Fsum       = sparse(varInd, onesCol, V.lat, maxVarInd, 1);
      FmeanMat   = Fsum(nonzeroIdx) ./ binCounter(nonzeroIdx);
      Sout.lat   = [Sout.lat, transpose(full(FmeanMat))];
    end,
    
    if (has.provenance),
      Fsum            = sparse(varInd, onesCol, V.provenance, maxVarInd, 1);
      FmeanMat        = Fsum(nonzeroIdx) ./ binCounter(nonzeroIdx);
      Sout.provenance = [Sout.provenance, transpose(full(FmeanMat))];
    end      

    Fsum       = sparse(varInd, onesCol, V.depth, maxVarInd, 1);
    FmeanMat   = Fsum(nonzeroIdx) ./ binCounter(nonzeroIdx);
    Sout.depth = [Sout.depth, transpose(full(FmeanMat))];

    Fsum       = sparse(varInd, onesCol, V.error, maxVarInd, 1);
    FmeanMat   = Fsum(nonzeroIdx) ./ binCounter(nonzeroIdx);
    Sout.error = [Sout.error, transpose(full(FmeanMat))];

    Fsum       = sparse(varInd, onesCol, V.value, maxVarInd, 1);
    FmeanMat   = Fsum(nonzeroIdx) ./ binCounter(nonzeroIdx);
    Sout.value = [Sout.value, transpose(full(FmeanMat))];

    Sout.type  = [Sout.type, ones([1 Nsuper]).*state_vars(n)];
   
    Sout.time  = [Sout.time, ones([1 Nsuper]).*Sinp.survey_time(m)];

    Sout.Nobs(m) = Sout.Nobs(m) + Nsuper;
    
  end,
end,

Sout.Ndatum = sum(Sout.Nobs);

%  Add other information fields from input structure. This
%  will facilitates writing of the new structure.

if (has.lonlat),
  Sout.do_longitude=true;
  Sout.do_latitude =true;
end,

if (isfield(Sinp,'title'));
  Sout.title = Sinp.title;
end,

if (isfield(Sinp,'state_flag_values'));
  Sout.state_flag_values = Sinp.state_flag_values;
end,

if (isfield(Sinp,'state_flag_meanings'));
  Sout.state_flag_meanings = Sinp.state_flag_meanings;
end,

if (isfield(Sinp,'origin_flag_values'));
  Sout.origin_flag_values = Sinp.state_flag_values;
end,

if (isfield(Sinp,'origin_flag_meanings'));
  Sout.origin_flag_meanings = Sinp.origin_flag_meanings;
end,

if (isfield(Sinp,'grd_file'));
  Sout.grd_file = Sinp.grd_file;
end,

if (isfield(Sinp,'grid_Lm_Mm_N'));
  Sout.grid_Lm_Mm_N = Sinp.grid_Lm_Mm_N;
end,

if (isfield(Sinp,'global_variables'));
  Sout.global_variables = Sinp.global_variables;
end,

if (isfield(Sinp,'global_provenance'));
  Sout.global_provenance = Sinp.global_provenance;
end,

if (isfield(Sinp,'global_sources'));
  Sout.global_sources = Sinp.global_sources;
end,

return
