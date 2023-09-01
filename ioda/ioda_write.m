function ioda_write(S, ncfile)

%
% IODA_WRITE:  Writes IODA observation NetCDF4 file
%
% ioda_write(S, ncfile)
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2023 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%
  
disp(['*** Writing IODA observations file:  ', ncfile]);

ncwrite(ncfile, 'Location', int32(zeros(size(S.longitude))));
ncwrite(ncfile, 'nvars', int32(zeros(size(S.iodaVarName))));

% Write out 'MetaData' Group variables.

ncwrite(ncfile, 'MetaData/dateTime', int64(S.dateTime));

ncwrite(ncfile, 'MetaData/date_time', S.date_time);

if (~isempty(S.depth))
  ncwrite(ncfile, 'MetaData/depth', S.depth);
end

ncwrite(ncfile, 'MetaData/latitude', S.latitude);
ncwrite(ncfile, 'MetaData/longitude', S.longitude);

if (~isempty(S.provenance))
  ncwrite(ncfile, 'MetaData/provenance', S.provenance);
end

record_number  = int32(ones(size(S.longitude)));
ncwrite(ncfile, 'MetaData/sequenceNumber', record_number);

% Read in 'EffectiveError' Group.

if (isfield(S, 'EffectiveError'))
  for i = 1:S.nvars
    Vname = ['EffectiveError/' S.iodaVarName{i}];
    ncwrite(ncfile, Vname, S.EffectiveError{i});
  end
end

% Read in 'EffectiveQC' Group.

if (isfield(S, 'EffectiveQC'))
  for i = 1:S.nvars
    Vname = ['EffectiveQC/' S.iodaVarName{i}];
    ncwrite(ncfile, Vname, int32(S.EffectiveError{i}));
  end
end

% Read in 'ObsBias' Group.

if (isfield(S, 'ObsBias'))
  for i = 1:S.nvars
    Vname = ['ObsBias/' S.iodaVarName{i}];
    ncwrite(ncfile, Vname, S.ObsBias{i});
  end
end

% Write out 'ObsError' Group variables.

if (isfield(S, 'ObsError'))
  for i = 1:S.nvars
    Vname = ['ObsError/' S.iodaVarName{i}];
    ncwrite(ncfile, Vname, S.ObsError{i});
  end
end

% Write out 'ObsValue' Group variables.

if (isfield(S, 'ObsValue'))
  for i = 1:S.nvars
    Vname = ['ObsValue/' S.iodaVarName{i}];
    ncwrite(ncfile, Vname, S.ObsValue{i});
  end
end

% Write out 'PreQC' Group variables.

if (isfield(S, 'PreQC'))
  for i = 1:S.nvars
    Vname = ['PreQC/' S.iodaVarName{i}];
    ncwrite(ncfile, Vname, int32(S.PreQC{i}));
  end
end

% Write out 'hofx' Group variables.

if (isfield(S, 'hofx'))
  for i = 1:S.nvars
    Vname = ['hofx/' S.iodaVarName{i}];
    ncwrite(ncfile, Vname, S.hofx{i});
  end
end

return
