function [f]=nc_read(fname,vname,tindex);

%
% NC_READ:  Read requested NetCDF variable
%
% [f]=nc_read(fname,vname,tindex)
%
% This function reads in a generic multi-dimensional field from a NetCDF
% file.  If only water points are available, this function fill the land
% areas with zero and returns full fields.
%
% On Input:
%
%    fname       NetCDF file name (character string)
%    vname       NetCDF variable name to read (character string)
%    tindex      Optional, time index to read (integer):
%                  *  If argument "tindex" is provided, only the requested
%                     time record is read when the variable has unlimitted
%                     dimension or the word "time" in any of its dimension
%                     names.
%
% On Output:
%
%    f           Field (scalar or array)
%
% calls:         nc_dim, nc_vinfo, nc_vname, ncread
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2009 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%----------------------------------------------------------------------------
% Inquire information from NetCDF file.
%----------------------------------------------------------------------------

% Inquire about file dimensions.

[dnames,dsizes]=nc_dim(fname);

for n=1:length(dsizes),
  name=deblank(dnames(n,:));
  switch name
    case 'xi_rho',
      Lr=dsizes(n);
    case 'xi_u',
      Lu=dsizes(n);
    case 'xi_v',
      Lv=dsizes(n);
    case 'eta_rho',
      Mr=dsizes(n);
    case 'eta_u',
      Mu=dsizes(n);
    case 'eta_v',
      Mv=dsizes(n);
    case 's_rho',
      Nr=dsizes(n);
    case 's_w',
      Nw=dsizes(n);
  end,
end,  
 
% Inquire about requested variable.

[vdnames,vdsizes,igrid]=nc_vinfo(fname,vname);

% Check if data is only available at water points.

is2d=0;
is3d=0;
water=0;
if (~isempty(vdsizes)),
  for n=1:length(vdsizes),
    name=deblank(vdnames(n,:));
    switch name
      case 'xy_rho',
        msknam='mask_rho';
        is2d=1; Im=Lr; Jm=Mr;
      case 'xy_u',
        msknam='mask_u';
        is2d=1; Im=Lu; Jm=Mu;
      case 'xy_v',
        msknam='mask_v';
        is2d=1; Im=Lv; Jm=Mv;
      case 'xyz_rho',
        msknam='mask_rho';
        is3d=1; Im=Lr; Jm=Mr; Km=Nr;
      case 'xyz_u',
        msknam='mask_u';
        is3d=1; Im=Lu; Jm=Mu; Km=Nr;
      case 'xyz_v',    
        msknam='mask_v';
        is3d=1; Im=Lv; Jm=Mv; Km=Nr;
      case 'xyz_w',    
        msknam='mask_rho';
        is3d=1; Im=Lr; Jm=Mr; Km=Nw;
    end,
  end,
end,

water=is2d | is3d;

% If water data only, read in Land/Sea mask.

got_mask=0;

if (water),
  [Vnames,nvars]=nc_vname(fname);
  got_mask=strmatch(msknam,Vnames);
  if (got_mask),
    mask=ncread(fname,msknam);
  else,
%   [fn,pth]=uigetfile(grid_file,'Enter grid NetCDF file...');
%   gname=[pth,fn];
    gname=input('Enter grid NetCDF file: ');
    mask=ncread(gname,msknam);
  end,
end,

%----------------------------------------------------------------------------
% Read in requested variable.
%----------------------------------------------------------------------------

if (water),

  if (nargin < 3),
    v=ncread(fname,vname);
  else
    v=ncread(fname,vname,tindex);
  end,
  [Npts,Ntime]=size(v);

  if (is2d),
    f=squeeze(ones([Im,Jm,Ntime])).*NaN;
    MASK=squeeze(repmat(mask,[1,1,Ntime]));
    ind=find(MASK > 0);
    f(ind)=v;
  elseif (is3d),    
    f=squeeze(ones([Im,Jm,Km,Ntime])).*NaN;
    MASK=squeeze(repmat(mask,[1,1,Km,Ntime]));
    ind=find(MASK > 0);
    f(ind)=v;
  end,

else,

  if (nargin < 3),
    f=ncread(fname,vname);
  else
    f=ncread(fname,vname,tindex);
  end,
  
end,

return

function [f]=ncread(fname,vname,tindex);

%
% NCREAD:  Internal routine to read requested NetCDF variable
%
% [f]=ncread(fname,vname,tindex)
%
% This function reads in a generic multi-dimensional field from a NetCDF
% file.
%
% On Input:
%
%    fname       NetCDF file name (character string)
%    vname       NetCDF variable name to read (character string)
%    tindex      Optional, time index to read (integer):
%                  *  If argument "tindex" is provided, only the requested
%                     time record is read if the variable has unlimitted
%                     dimension or the word "time" in any of its dimension
%                     names.
%
% On Output:
%
%    f           Field (scalar or array)
%

%  Set-up printing information switch.

global IPRINT

if (isempty(IPRINT)),
  IPRINT=1;
end,

%  Activate switch for reading specific record.

time_rec=0;
if (nargin > 2),
  time_rec=1;
end,

% Open NetCDF file.

[ncid]=mexnc('ncopen',fname,'nc_nowrite');
if (ncid == -1),
  error(['NCREAD: ncopen - unable to open file: ' fname])
  return
end,

% Supress all error messages from NetCDF.

[status]=mexnc('setopts',0);

% Define NetCDF parameters.

[NC_BYTE  ]=mexnc('parameter','nc_byte');
[NC_CHAR  ]=mexnc('parameter','nc_char');
[NC_DOUBLE]=mexnc('parameter','nc_double');
[NC_FLOAT ]=mexnc('parameter','nc_float');
[NC_INT   ]=mexnc('parameter','nc_int');
[NC_SHORT ]=mexnc('parameter','nc_short');

%----------------------------------------------------------------------------
% Inquire about requested variable.
%----------------------------------------------------------------------------

% Get variable ID.

[varid]=mexnc('ncvarid',ncid,vname);
if (varid < 0),
  [status]=mexnc('ncclose',ncid);
  nc_inq(fname);
  disp('  ');
  error(['NCREAD: ncvarid - cannot find variable: ',vname])
  return
end,

% Inquire about unlimmited dimension.

[ndims,nvars,natts,recdim,status]=mexnc('ncinquire',ncid);
if (status == -1),
  error(['NCREAD: ncinquire - cannot inquire file: ',fname])
end,

% Get information about requested variable.

[vname,nctype,nvdims,dimids,nvatts,status]=mexnc('ncvarinq',ncid,varid);
if (status == -1),
  error(['NCREAD: ncvarinq - unable to inquire about variable: ',vname])
end,

% Inquire about the _FillValue attribute.

got_FillValue=0;

for i = 0:nvatts-1
  [attnam,status]=mexnc('ncattname',ncid,varid,i);
  if (status == -1)
    error(['NCREAD: ncattname: error while inquiring attribute ' num2str(i)])
  end,
  lstr=length(attnam);
  if (strncmp(attnam(1:lstr),'_FillValue',10)),
    if (nctype == NC_DOUBLE),
      [spval,status]=mexnc('get_att_double',ncid,varid,'_FillValue'); 
    elseif (nctype == NC_FLOAT),
      [spval,status]=mexnc('get_att_float' ,ncid,varid,'_FillValue');
    elseif (nctype == NC_INT),
      [spval,status]=mexnc('get_att_int'   ,ncid,varid,'_FillValue');
    else
      [spval,status]=mexnc('ncattget'      ,ncid,varid,'_FillValue');
    end,
    if (status == -1),
      error(['NCREAD: ncattget error while reading _FillValue attribute'])
    end,
    got_FillValue=1;
  end,
end,

% Inquire about dimensions.

index=0;
for n=1:nvdims
  [name,dsize,status]=mexnc('ncdiminq',ncid,dimids(n));
  if (status == -1),
    error(['NCREAD: ncdiminq - unable to inquire about dimension ID: ',...
          num2str(dimids(n))])
  else
    lstr=length(name);
    dimnam(n,1:lstr)=name(1:lstr);
    dimsiz(n)=dsize;
    start(n)=0;
    count(n)=dsize;
    if ((dimids(n) == recdim) | ~isempty(findstr(name,'time'))),
      index=n;
    end,
  end,
end,

%  It reading specific time record, reset variable bounds.

nvdim=nvdims;
if (time_rec & (index > 0)),
  start(index)=tindex-1;
  count(index)=1;
  nvdims=nvdims-1;
else
  index=0;
end,

%----------------------------------------------------------------------------
% Read in requested variable.
%----------------------------------------------------------------------------

%  Read in scalar.

if (nvdim == 0),

  if (nctype == NC_DOUBLE),
    [f,status]=mexnc('get_var_double',ncid,varid);
  elseif (nctype == NC_FLOAT),
    [f,status]=mexnc('get_var_float' ,ncid,varid);
  elseif (nctype == NC_INT),
    [f,status]=mexnc('get_var_int'   ,ncid,varid);
  elseif (nctype == NC_CHAR),
    [f,status]=mexnc('get_var_text'  ,ncid,varid);
  else,
    [f,status]=mexnc('ncvarget1'     ,ncid,varid,[0]);
  end,

  if (status == -1),
    error(['NCREAD: ncvarget1 - error while reading: ',vname])
  end,

%  Read in a multidemensional array.

else,

  if (nctype == NC_DOUBLE),
    [f,status]=mexnc('get_vara_double',ncid,varid,start,count);
  elseif (nctype == NC_FLOAT),
    [f,status]=mexnc('get_vara_float' ,ncid,varid,start,count);
  elseif (nctype == NC_INT),
    [f,status]=mexnc('get_vara_int'   ,ncid,varid,start,count);
  else,
    [f,status]=mexnc('ncvarget'       ,ncid,varid,start,count);
  end,

  if (status == -1),
    error(['NCREAD: ncvarget - error while reading: ',vname])
  end,

  if (nvdims == 3),
    if (length(start) == 3),
      f=reshape(f,[count(3),count(2),count(1)]);
    elseif (length(start) == 4),
      f=reshape(f,[count(4),count(3),count(2)]);
    end,
  end,
  
  if (nvdims == 4),
    if (length(start) == 4),
      f=reshape(f,[count(4),count(3),count(2),count(1)]);
    elseif (length(start) == 5),
      f=reshape(f,[count(5),count(4),count(3),count(2)]);
    end,  
  end,
  
end,

% Print information about variable.

if (IPRINT),
  if (nvdims > 0),
    disp('  ')
    disp([vname ' has the following dimensions (input order):']);
    disp('  ')
    for n=1:nvdim,
      s=['           '  int2str(n) ') ' dimnam(n,:) ' = ' int2str(dimsiz(n))];
      disp(s);
    end,
    disp('  ')
    disp([vname ' loaded into an array of size:  [' int2str(size(f)) ']']);
    disp('  ')
  else
    disp('  ')
    disp([vname ' is a scalar and has a value of ',num2str(f)]);
    disp('  ')
  end,
end,

%----------------------------------------------------------------------------
%  Replace _FillValue with NaNs.
%----------------------------------------------------------------------------

if (got_FillValue),
  ind=find(f >= spval);
  if (~isempty(ind)),
    f(ind)=NaN;
  end,
end,

%----------------------------------------------------------------------------
%  Convert data to double precision.
%----------------------------------------------------------------------------

if (nctype == NC_FLOAT | nctype == NC_INT),
  f=double(f);
end,

%----------------------------------------------------------------------------
% Close NetCDF file.
%----------------------------------------------------------------------------

[status]=mexnc('ncclose',ncid);
if (status == -1),
  error(['NCREAD: ncclose - unable to close NetCDF file.'])
end

return
