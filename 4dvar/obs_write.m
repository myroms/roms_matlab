function [status]=obs_write(ncfile,S);

%
% OBS_WRITE:  Writes 4D-Var observation data into existing NetCDF file
%
% [status]=obs_write(S)
%
% This function writes all observation data in structure array S into
% a existing ROMS 4D-Var observation NetCDF file.
%
% On Input:
%
%    ncfile  4D-Var Observations NetCDF file name (string)
%
%    S       Observations data (structure array):
%
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
%            The following optional variables will be written if available:
%
%              S.provenance     observation origin
%              S.lon            observation longitude
%              S.lat            observation latitude
% 
% On Output:
%
%    status  Error flag
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

status=[];

%----------------------------------------------------------------------------
%  Check dimensions of input NetCDF file.
%----------------------------------------------------------------------------

[Dnames,Dsizes,unlimit]=nc_dim(ncfile);

for n=1:length(Dsizes),
  name=deblank(Dnames(n,:));
  switch name
    case 'survey'
      Nsurvey=Dsizes(n);
    case 'state_variable'
      Nstate=Dsizes(n);
    case 'datum'
      Ndatum=Dsizes(n);
  end,
end,

%  Check 'survey' and 'state_variable' dimensions.

if (Nsurvey ~= length(S.survey_time)),
  error(['OBS_WRITE: dimension mismatch between data and file, ', ...
	 'survey = ',num2str(Nsurvey),' ',num2str(length(S.survey_time))]);
end,

if (Nstate ~= length(S.variance)),
  error(['OBS_WRITE: dimension mismatch between data and file, ', ...
	 'state_variable = ',num2str(Nstate),' ', ...
	 num2str(length(S.variance))]);
end,
	 
%  Check size of the datum dimension in input NetCDF file. This dimension
%  can be unlimited or not. If not unlimited, issue an error if the size
%  of the observation array is not the same as the size of the dimension
%  in input NetCDF file. If the datum dimension is unlimited, it will
%  handle any size of the observation vector.

if (unlimit == -1),
  if (Ndatum ~= length(S.value)),
    error(['OBS_WRITE: dimension mismatch between data and file, ', ...
	   'datum = ',num2str(Ndatum),' ',num2str(length(S.value))]);
  end,
end,

disp(' ');

%----------------------------------------------------------------------------
%  Write all available variables in the structure.
%----------------------------------------------------------------------------

[Vnames,nvars]=nc_vname(ncfile);

notwritten=[];

for n=1:nvars,
  name=deblank(Vnames(n,:));
  switch name
    case 'spherical'
      if (isfield(S,'spherical')),
        status=nc_write(ncfile,'spherical',S.spherical);
      else,
	disp(['   field ''spherical''   not found in input structure']);
        notwritten=[notwritten 'spherical '];
      end,
    case 'Nobs'
      if (isfield(S,'Nobs')),
        status=nc_write(ncfile,'Nobs',S.Nobs);
      else,
	disp(['   field ''Nobs''        not found in input structure']);
        notwritten=[notwritten 'Nobs '];
      end,
    case 'survey_time'
      if (isfield(S,'survey_time')),
        status=nc_write(ncfile,'survey_time',S.survey_time);
      else,
	disp(['   field ''survey_time'' not found in input structure']);
        notwritten=[notwritten 'survey_time '];
      end,
    case 'obs_variance'
       if (isfield(S,'variance')),
        status=nc_write(ncfile,'obs_variance',S.variance);
      else,
	disp(['   field ''variance''    not found in input structure']);
        notwritten=[notwritten 'obs_variance '];
      end,
    case 'obs_type'
      if (isfield(S,'type')),
        status=nc_write(ncfile,'obs_type',S.type);
      else,
	disp(['   field ''type''        not found in input structure']);
        notwritten=[notwritten 'obs_type '];
      end,
    case 'obs_provenance'
      if (isfield(S,'provenance')),
        status=nc_write(ncfile,'obs_provenance',S.provenance);
      else,
	disp(['   field ''provenance''  not found in input structure']);
        notwritten=[notwritten 'obs_provenance '];
      end,
    case 'obs_time'
      if (isfield(S,'time')),
        status=nc_write(ncfile,'obs_time',S.time);
      else,
	disp(['   field ''time''        not found in input structure']);
        notwritten=[notwritten 'obs_time '];
      end,
    case 'obs_lon'
      if (isfield(S,'lon')),
        status=nc_write(ncfile,'obs_lon',S.lon);
      else,
	disp(['   field ''lon''         not found in input structure']);
        notwritten=[notwritten 'obs_lon '];
      end,
    case 'obs_lat'
      if (isfield(S,'lat')),
        status=nc_write(ncfile,'obs_lat',S.lat);
      else,
	disp(['   field ''lat''         not found in input structure']);
        notwritten=[notwritten 'obs_lat '];
      end,
    case 'obs_depth'
      if (isfield(S,'depth')),
        status=nc_write(ncfile,'obs_depth',S.depth);
      else,
	disp(['   field ''depth''       not found in input structure']);
        notwritten=[notwritten 'obs_depth '];
      end,
    case 'obs_Xgrid'
      if (isfield(S,'Xgrid')),
        status=nc_write(ncfile,'obs_Xgrid',S.Xgrid);
      else,
	disp(['   field ''Xgrid''       not found in input structure']);
        notwritten=[notwritten 'obs_Xgrid '];
      end,
    case 'obs_Ygrid'
      if (isfield(S,'Ygrid')),
        status=nc_write(ncfile,'obs_Ygrid',S.Ygrid);
      else,
	disp(['   field ''Ygrid''       not found in input structure']);
        notwritten=[notwritten 'obs_Ygrid '];
      end,
    case 'obs_Zgrid'
      if (isfield(S,'Zgrid')),
        status=nc_write(ncfile,'obs_Zgrid',S.Zgrid);
      else,
	disp(['   field ''Zgrid''       not found in input structure']);
        notwritten=[notwritten 'obs_Zgrid '];
      end,
    case 'obs_error'
      if (isfield(S,'error')),
        status=nc_write(ncfile,'obs_error',S.error);
      else,
	disp(['   field ''error''       not found in input structure']);
        notwritten=[notwritten 'obs_error '];
      end,
    case 'obs_value'
      if (isfield(S,'value')),
        status=nc_write(ncfile,'obs_value',S.value);
      else,
	disp(['   field ''value''       not found in input structure']);
        notwritten=[notwritten 'obs_value '];
      end,
  end,
end,

if (~isempty(notwritten)),
  disp(' ');
  disp(['   Following variables were not written:  ', notwritten]);
end,

return