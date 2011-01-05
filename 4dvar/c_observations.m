function [status]=c_observations(S,file);

%
% C_OBSERVATIONS:  Creates 4D-Var observations NetCDF file
%
% [status]=c_observations(S,file)
%
% This function creates ROMS 4D-Var observation NetCDF file using specified
% in structure array, S.
%
% On Input:
%
%    S     Observations file creation parameters (structure array):
%
%            S.ncfile                NetCDF file name (string)
%            S.spherical             spherical grid switch
%            S.Nstate                number of state variables
%            S.Nsurvey               number of data surveys 
%            S.Ndatum                number of observations
%            S.state_flag_values     obs_type 'flag_values' attribute
%            S.state_flag_meanings   obs_type 'flag_meanings attribute
%            S.origin_flag_values    obs_provenance 'flag_values' attribute
%            S.origin_flag_meanings  obs_provenance 'flag_meaning' attribute
%
%          Optional fields:
%
%            S.title                 application title
%            S.grd_file              'grd_file' global attribute
%            S.global_variables      'state_variables' global attribute
%            S.global_provenance     'obs_provenance' global attribute
%            S.global_sources        'obs_sources' global attribute
%
%            S.do_longitude          switch to define longitude locations
%            S.do_latitude           switch to define latitude  locations
%
%    file  Output observation file name (string, OPTIONAL). If specified,
%            it creates this file instead the one in S.ncfile.
%
% On Output:
%
%    status  Error flag.
%
% Examples showing how to assign some of the structure fields:
%
%    newline=sprintf('\n');
%
%    S.state_flag_values=[1:1:7];
%
%    S.state_flag_meanings={'zeta', ...
%                           'ubar', ...
%                           'vbar', ...
%                           'u', ...
%                           'v', ...
%                           'temperature', ...
%                           'salinity'};
%
%    S.origin_flag_values=[1:1:11];
%
%    S.origin_flag_meanings={'gridded_AVISO_SLA', ...
%                            'blended_SST', ...
%                            'XBT_Met_Office', ...
%                            'CTD_temperature_Met_Office', ...
%                            'CTD_salinity_Met_Office', ...
%                            'ARGO_temperature_Met_Office', ...
%                            'ARGO_salinity_Met_Office', ...
%                            'CTD_temperature_CalCOFI', ...
%                            'CTD_salinity_CalCOFI', ...
%                            'CTD_temperature_GLOBEC', ...
%                            'CTD_salinity_GLOBEC'};
%
%    S.state_variables=[newline, ...
%      '1: free-surface (m) ', newline, ...
%      '2: vertically integrated u-momentum component (m/s) ', newline, ...
%      '3: vertically integrated v-momentum component (m/s) ', newline, ...
%      '4: u-momentum component (m/s) ', newline, ...
%      '5: v-momentum component (m/s) ', newline, ...
%      '6: potential temperature (Celsius) ', newline, ...
%      '7: salinity (nondimensional)'];
%
%    S.obs_provenance=[newline, ...
%       ' 1: gridded AVISO sea level anomaly ', newline, ...
%       ' 2: blended satellite SST ', newline, ...
%       ' 3: XBT temperature from Met Office ', newline, ...
%       ' 4: CTD temperature from Met Office ', newline, ...
%       ' 5: CTD salinity from Met Office ', newline, ...
%       ' 6: ARGO floats temperature from Met Office ', newline, ...
%       ' 7: ARGO floats salinity from Met Office ', newline, ...
%       ' 8: CTD temperature from CalCOFI ', newline, ...
%       ' 9: CTD salinity from CalCOFI ', newline, ...
%       '10: CTD temperature from GLOBEC ', newline, ...
%       '11: CTD salinity from GLOBEC'];
%
  
% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2011 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%----------------------------------------------------------------------------
%  Set some NetCDF parameters.
%----------------------------------------------------------------------------

[ncglobal]=mexnc('parameter','nc_global');
[ncbyte  ]=mexnc('parameter','nc_byte');
[ncchar  ]=mexnc('parameter','nc_char');
[ncshort ]=mexnc('parameter','nc_short');
[ncint   ]=mexnc('parameter','nc_int');
[ncfloat ]=mexnc('parameter','nc_float');
[ncdouble]=mexnc('parameter','nc_double');

%----------------------------------------------------------------------------
%  Get error covariance standard deviation creation parameters.
%----------------------------------------------------------------------------

if (nargin < 2),
  if (isfield(S,'ncfile')),
    ncfile=S.ncfile;
  else,
    error([ 'C_OBSERVATIONS - Cannot find file name field: ncname, ', ...
	    'in structure array S']);
  end,
else,
  ncfile=file;    
end,

if (isfield(S,'spherical')),
  spherical=S.spherical;
else,
  spherical=0;
end,

%----------------------------------------------------------------------------
%  Set dimensions.
%----------------------------------------------------------------------------

Dname.survey ='survey';          Dsize.survey = S.Nsurvey;
Dname.state  ='state_variable';  Dsize.state  = S.Nstate;
Dname.datum  ='datum';           Dsize.datum  = S.Ndatum;

%----------------------------------------------------------------------------
%  Set all posible variables names.
%----------------------------------------------------------------------------

Vname.spherical  = 'spherical';
Vname.Nobs       = 'Nobs';
Vname.survey     = 'survey_time';
Vname.variance   = 'obs_variance';
Vname.type       = 'obs_type';
Vname.provenance = 'obs_provenance';
Vname.time       = 'obs_time';
Vname.depth      = 'obs_depth';
Vname.lon        = 'obs_lon';
Vname.lat        = 'obs_lat';
Vname.Xgrid      = 'obs_Xgrid';
Vname.Ygrid      = 'obs_Ygrid';
Vname.Zgrid      = 'obs_Zgrid';
Vname.error      = 'obs_error';
Vname.value      = 'obs_value';

%----------------------------------------------------------------------------
%  Create 4D-Var observation NetCDF file.
%----------------------------------------------------------------------------

disp(' ');
disp(['*** Creating observations file:  ', ncfile]);

[ncid,status]=mexnc('create',ncfile,'nc_clobber');
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_OBSERVATIONS: CREATE - unable to create file: ', ncname]);
  return
end,

%----------------------------------------------------------------------------
%  Define dimensions.
%----------------------------------------------------------------------------

[did.survey,status]=mexnc('def_dim',ncid,Dname.survey,Dsize.survey); 
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_OBSERVATIONS: DEF_DIM - unable to define dimension: ', ...
	  Dname.survey]);
  return
end,

[did.state,status]=mexnc('def_dim',ncid,Dname.state,Dsize.state); 
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_OBSERVATIONS: DEF_DIM - unable to define dimension: ', ...
	  Dname.state]);
  return
end,

[did.datum,status]=mexnc('def_dim',ncid,Dname.datum,Dsize.datum); 
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_OBSERVATIONS: DEF_DIM - unable to define dimension: ', ...
	  Dname.datum]);
  return
end,

%----------------------------------------------------------------------------
%  Create global attributes.
%----------------------------------------------------------------------------

str='ROMS observations';
lstr=length(str);
[status]=mexnc('put_att_text',ncid,ncglobal,'type',ncchar,lstr,str);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:', ...
	  ' type.']);
  return
end,

if (isfield(S,'title')),
  lstr=length(S.title);
  [status]=mexnc('put_att_text',ncid,ncglobal,'title',ncchar,lstr,S.title);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error([ 'C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:' ...
	    ' title.']);
    return
  end,
end,

str='CF-1.4';
lstr=length(str);
[status]=mexnc('put_att_text',ncid,ncglobal,'Conventions',ncchar,lstr,str);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:', ...
	  ' Conventions.']);
  return
end,

if (isfield(S,'grd_file')),
  lstr=length(S.grd_file);
  [status]=mexnc('put_att_text',ncid,ncglobal,'grd_file',ncchar,lstr, ...
		 S.grd_file);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error([ 'C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:', ...
	   ' grd_file.']);
    return
  end,
end,

if (isfield(S,'grid_Lm_Mm_N')),
  nval=length(S.grid_Lm_Mm_N);
  [status]=mexnc('put_att_int',ncid,ncglobal,'grid_Lm_Mm_N',ncint,nval, ...
		 int32(S.grid_Lm_Mm_N));
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error([ 'C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:', ...
	   ' grd_Lm_Mm_N.']);
    return
  end,
end,

if (isfield(S,'global_variables')),
  lstr=length(S.global_variables);
  [status]=mexnc('put_att_text',ncid,ncglobal,'state_variables',ncchar, ...
                 lstr,S.global_variables);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error([ 'C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:' ...
	    ' state_variables.']);
    return
  end,
end,

if (isfield(S,'global_provenance')),
  lstr=length(S.global_provenance);
  [status]=mexnc('put_att_text',ncid,ncglobal,'obs_provenance',ncchar, ...
                 lstr,S.global_provenance);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error([ 'C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:' ...
	    ' obs_provenance.']);
    return
  end,
end,

str='squared state variable units';
lstr=length(str);
[status]=mexnc('put_att_text',ncid,ncglobal,'variance_units',ncchar, ...
	       lstr,str);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribure:', ...
	  ' variance_units.']);
  return
end,

if (isfield(S,'global_sources')),
  lstr=length(S.global_sources);
  [status]=mexnc('put_att_text',ncid,ncglobal,'obs_sources',ncchar, ...
                 lstr,S.global_sources);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error([ 'C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:' ...
	    'obs_sources.']);
    return
  end,
end,

str=['4D-Var observations, ',date_stamp];
lstr=length(str);
[status]=mexnc('put_att_text',ncid,ncglobal,'history',ncchar,lstr,str);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:', ...
	  ' history.']);
  return
end,

%----------------------------------------------------------------------------
%  Define configuration variables.
%----------------------------------------------------------------------------

% Define spherical switch.

Var.name          = Vname.spherical;
Var.type          = ncint;
Var.dimid         = [];
Var.long_name     = 'grid type logical switch';
Var.flag_values   = [0 1];
Var.flag_meanings = ['Cartesian', blanks(1), ...
                     'spherical'];
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

% Define observation variables.

Var.name          = Vname.Nobs;
Var.type          = ncint;
Var.dimid         = [did.survey];
Var.long_name     = 'number of observations with the same survey time';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.survey;
Var.type          = ncdouble;
Var.dimid         = [did.survey];
Var.long_name     = 'survey time';
Var.units         = 'days';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.variance;
Var.type          = ncdouble;
Var.dimid         = [did.state];
Var.long_name     = 'global temporal and spatial observation variance';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.type;
Var.type          = ncint;
Var.dimid         = [did.datum];
Var.long          = 'model state variable associated with observations';
Var.flag_values   = S.state_flag_values;
Var.flag_meanings = S.state_flag_meanings;
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.provenance;
Var.type          = ncint;
Var.dimid         = [did.datum];
Var.long_name     = 'observation origin';
Var.flag_values   = S.origin_flag_values;
Var.flag_meanings = S.origin_flag_meanings;
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.time;
Var.type          = ncdouble;
Var.dimid         = [did.datum];
Var.long_name     = 'time of observation';
Var.units         = 'days';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

if (isfield(S,'do_longitude')),
  if (S.do_longitude),
    Var.name          = Vname.lon;
    Var.type          = ncdouble;
    Var.dimid         = [did.datum];
    Var.long_name     = 'observation longitude';
    Var.units         = 'degrees_east';
    Var.standard_name = 'longitude';
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
  end,
end,

if (isfield(S,'do_latitude')),
  if (S.do_latitude),
    Var.name          = Vname.lat;
    Var.type          = ncdouble;
    Var.dimid         = [did.datum];
    Var.long_name     = 'observation latitude';
    Var.units         = 'degrees_north';
    Var.standard_name = 'latitude';
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
  end,
end,

Var.name          = Vname.depth;
Var.type          = ncdouble;
Var.dimid         = [did.datum];
Var.long_name     = 'depth of observation';
Var.units         = 'meters';
Var.negative      = 'downwards';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.Xgrid;
Var.type          = ncdouble;
Var.dimid         = [did.datum];
Var.long_name     = 'observation fractional x-grid location';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.Ygrid;
Var.type          = ncdouble;
Var.dimid         = [did.datum];
Var.long_name     = 'observation fractional y-grid location';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.Zgrid;
Var.type          = ncdouble;
Var.dimid         = [did.datum];
Var.long_name     = 'observation fractional z-grid location';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.error;
Var.type          = ncdouble;
Var.dimid         = [did.datum];
Var.long_name     = 'observation error covariance';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.value;
Var.type          = ncdouble;
Var.dimid         = [did.datum];
Var.long_name     = 'observation value';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

%----------------------------------------------------------------------------
%  Leave definition mode and close NetCDF file.
%----------------------------------------------------------------------------

[status]=mexnc('enddef',ncid);
if (status == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_OBSERVATIONS: ENDDEF - unable to leave definition mode.']);
  return
end,

[status]=mexnc('close',ncid);
if (status == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_OBSERVATIONS: CLOSE - unable to close NetCDF file: ', ncname]);
  return
end,

return

