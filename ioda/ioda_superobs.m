function [Sout]=ioda_superobs(InpName, OutName)

%
% IODA_SUPEROBS:  Creates super observations when necessary
%
% [Sout]=ioda_superobs(InpName, OutName)
%
% This function checks the provided observation data and creates
% super observations when there are more than one meassurement of
% the same state variable per grid cell. At input, Sinp is either a
% 4D-Var observation NetCDF file or data structure.
%
% An additional field (Sout.std) is added to the output observation
% structure containing the standard deviation of the binning which
% can be used as observation error. You may choose to assign this
% value as observation error before writing to NetCDF file:
%
%    Sout.error = max(Sout.error, Sout.std)
%
% That is, the larger observation variace is chosen. A zero value
% for "Sout.std" indicates that the observation did not required
% binning.  Use "c_observations.m" to create NetCDF file and
% "obs_write" to write out data.
%
% On Input:
%
%    InpName  Input IODA NetCDF-4 filename (string)
%    OutName    Output thinned IODA NetCDF-4 filename (string)
%
% On Output:
%
%    S        Binned observations data (structure array):
%
  
% git $Id$
%=======================================================================%
%  Copyright (c) 2002-2026 The ROMS Group                               %
%    Licensed under a MIT/X style license                               %
%    See License_ROMS.md                            Hernan G. Arango    %
%=======================================================================%

%  Read observations if 'Sinp' is a 4D-Var observation NetCDF file.

Sinp=ioda_read(InpName);

%  Check if 'depth' and 'provenace' are available.

has.depth = false;
if (isfield(Sinp, 'depth'))
  has.depth = true;
end

has.provenance = false;
if (isfield(Sinp, 'provenance'))
  has.provenance = true;
end


%  Set observations dynamical fields (cell array) in structure, S.

field_list = {'x_grid', 'y_grid', 'longitude', 'latitude',            ...
              'ObsError', 'ObsValue'};

if (has.depth)
  field_list = [field_list, 'depth', 'z_grid'];
end

if (has.provenance)
  field_list = [field_list, 'provenance'];
end

%  Insure that the vector fields in the input structure have the
%  singleton as the first dimension to allow vector concatenation.

for value = field_list
  field = char(value);
  if (size(Sinp.(field),1) > 1)
    Sinp.(field) = transpose(Sinp.(field));
  end
end

if (size(Sinp.surveyTime,1) > 1)
  Sinp.surveyTime = transpose(Sinp.surveyTime);
end

%  Add binning standard deviation field to input structure. Initialize
%  to input observation error.

Sinp.std = zeros(size(Sinp.ObsError));

%----------------------------------------------------------------------------
%  Find observations associated with the same state variable.
%----------------------------------------------------------------------------

%  Initialize output structure.

Sout = Sinp;

Sout.nlocs          = 0;
Sout.dateTime       = [];
Sout.sequenceNumber = [];
Sout.surveyIndex    = [];

for value = field_list
  field = char(value);             % convert from cell to string
  Sout.(field) = [];               % initilize to empty
end

Sout.std = [];                     % binning standard deviation

if (isfield(Sinp, 'PreQC'))
  Sout.PreQC = [];
end

%----------------------------------------------------------------------------
%  Compute super observations when needed.
%----------------------------------------------------------------------------

disp(blanks(1));
disp(['*** Thinning input observations file:  ', InpName]);

for m=1:Sinp.nsurvey  %%% SURVEY TIME LOOP %%%

%  Extract locations with the same survey time.

  ind_t=find(Sinp.dateTime == Sinp.surveyTime(m));

                                              % not included in dynamic
  T.dateTime = Sinp.dateTime(ind_t);          % fields for efficiency
  
  for value = field_list
    field = char(value);                      % convert from cell to string
    T.(field) = Sinp.(field);                 % initilize to Sinp structure
  end
  T.std = Sinp.std;                           % binning standard deviation
  
%  Set binning parameters. The processing is done in fractional (x,y,z) grid
%  locations

  Xmin = min(T.x_grid);
  Xmax = max(T.x_grid);

  Ymin = min(T.y_grid);
  Ymax = max(T.y_grid);

  Zmin = min(T.z_grid);
  Zmax = max(T.z_grid);

  dx = 1.0;
  dy = 1.0;
  dz = 1.0;
    
%  Compute the index in each dimension of the grid cell in which the
%  observation is located.

  Xbin = 1.0 + floor((T.x_grid - Xmin) ./ dx);
  Ybin = 1.0 + floor((T.y_grid - Ymin) ./ dy);
  Zbin = 1.0 + floor((T.z_grid - Zmin) ./ dz);

%  Similarly, compute the maximum averaging grid size.

  Xsize = 1.0 + floor((Xmax - Xmin) ./ dx);
  Ysize = 1.0 + floor((Ymax - Ymin) ./ dy);
  Zsize = 1.0 + floor((Zmax - Zmin) ./ dz);
    
  matsize = [Ysize, Xsize, Zsize];

%  Combine the indices in each dimension into one index. It is like stacking
%  all the matrix in one column vector.

  varInd    = transpose(sub2ind(matsize, Ybin, Xbin, Zbin));
  onesCol   = transpose(ones(size(varInd)));

%  Accumulate values in bins using "accumarray" function. Count how many
%  observations fall in each bin.

  count = accumarray(varInd, onesCol, [], @sum, [], true);
    
%  Bins with no observations are not keept. Index vector "isdata" will be
%  used in the conversion from the sparse output from "accumarray" to the
%  corresponding full vector assigned to Sout.

  isdata = find(count ~= 0);
  Nsuper = length(isdata);
    
%  Loop through list of fields that require binning. The "accumarray"
%  function take the sum of all values having the same bin index. The last
%  argument activates sparse. It turns out that method below is (somewhat)
%  faster than asking directly for the mean, thus:
%
%     binned = accumarray(varInd, V.(field), [], @mean, [], true);

  for fval = field_list
    field = char(fval);
      
    if (strcmp(field, 'ObsValue'))
      for n=1:Sout.nvars
        Value = T.ObsValue{n};
        binned = accumarray(varInd, Value, [], @sum, [], true);
        binned = full(binned(isdata) ./ count(isdata));
        Sout.ObsValue{n} = binned;     

        Vmean = binned;      
        binned = accumarray(varInd, Value.^2, [], @sum, [], true);
        binned = full(binned(isdata) ./ count(isdata)) - Vmean .^ 2;

        Sout.std{n} = sqrt(binned);
      end      
    elseif (strcmp(field, 'ObsError'))
      for n=1:Sout.nvars
        Error  = T.ObsError{n};
        binned = accumarray(varInd, Error, [], @sum, [], true);
        binned = full(binned(isdata) ./ count(isdata));
        Sout.ObsError{n} = binned;     
      end      
    else
      binned = accumarray(varInd, T.(field), [], @sum, [], true);
      binned = full(binned(isdata) ./ count(isdata));
      Sout.(field) = binned;
    end	
  end

%  Set "time" fields for binned observations.

  Sout.dateTime  = [Sout.dateTime, ones([1 Nsuper]).*Sinp.surveyTime(m)];

  Sout.nlocs = Sout.nlocs + Nsuper;

%  Make sure that "provenance", if available, is not a fractional number.
%  Binning provenance is weird here. Hopefully, this is always a full
%  number since the extrategy is to create one single observation file
%  per dataset or instrument.
    
  if (has.provenance)
    Sout.provenance = floor(Sout.provenance);
  end
  Sout.sequenceNumber = 1:1:Sout.nlocs;
  
end  %%% end of SURVEY TIME LOOP    %%%

%  Determine 'surveyIndex' defined ad observation survey time indices
%  as they appear in 'dateTime'.

[survey,IA,IC] = unique(Sout.dateTime, 'stable');
Sout.surveyTime  = int64(survey);
Sout.surveyIndex = int32(IA);

if (isfield(Sout, 'PreQC'))
  QC=zeros(size(Sout.longitude));
  for n=1:Sout.nvars
    Sout.PreQC{n} = QC;
  end
end

%  Assign binning standard deviation error as observation error.

for n=1:Sout.nvars
  Sout.ObsError{n} = max(Sout.ObsError{n}, Sout.std{n});
end

%------------------------------------------------------------------------
%  Create thinned (binned) output IODA NetCDF-4 file.
%------------------------------------------------------------------------

create_ioda_obs(Sout, OutName);

ioda_write(Sout, OutName);

return
