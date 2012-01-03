function f=nc_read(Fname,Vname,Tindex,ReplaceValue,PreserveType)

%
% NC_READ:  Read requested NetCDF variable
%
% f=nc_read(Fname,Vname,Tindex,ReplaceValue,PreserveType)
%
% This function reads in a generic multi-dimensional field from a NetCDF
% file, URL file, or URL aggregated files.  If only the water points are
% available, this function fill the land areas with zero and returns full
% fields.
%
% On Input:
%
%    Fname         NetCDF file name or URL name (character string)
%    Vname         NetCDF variable name to read (character string)
%    Tindex        Optional, time record index to read (integer):
%                    If Tindex is provided, only the requested time record
%                    is read when the variable has unlimitted dimension or
%                    the word "time" in any of its dimension names.
%                    Otherwise, provide an empty [] argument.
%    ReplaceValue  Optional, value to use when "_FillValue" or "missing_value"
%                    attribute is found in variable. If not provided, a zero
%                    value will used.
%                    In some instances, like plotting, it is advantageous
%                    to set ReplaceValue = NaN to visualize better the land
%                    masking or the missing data.
%    PreserveType  Switch to preserve numerical data type. If false,
%                    convert numerical data to double precision. It
%                    has not effect data is already in double
%                    precision.
%
% On Output:
%
%    f             Field (scalar or array)
%
% calls:           nc_dim, nc_url, nc_vinfo, nc_Vname
%                  nc_info, nc_vargetr (SNCTOOLS java interface to OpenDAT)
%                  read_nc, read_nc_url (private functions below)
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% If the PreserveType switch is not provide, set default to convert data
% to double precision.

if (nargin < 5),
  PreserveType=false;
end

% If ReplaceValue is not provided, use zero as a fill value.

if (nargin < 4),
  ReplaceValue=0;
end

% If Tindex is not provided, use empty value.

if (nargin < 3)
  Tindex=[];
end

% Check if input file is URL from an OpenDAP server.

url=nc_url(Fname);

%----------------------------------------------------------------------------
% Inquire information from NetCDF file.
%----------------------------------------------------------------------------

% Inquire about file dimensions.

[dnames,dsizes]=nc_dim(Fname);

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
  end
end  
 
% Inquire about requested variable.

[vdnames,vdsizes,igrid]=nc_vinfo(Fname,Vname);

% Check if data is only available at water points.

is2d=false;
is3d=false;
water=false;

if (~isempty(vdsizes)),
  for n=1:length(vdsizes),
    name=deblank(vdnames(n,:));
    switch name
      case 'xy_rho',
        msknam='mask_rho';
        is2d=true; Im=Lr; Jm=Mr;
      case 'xy_u',
        msknam='mask_u';
        is2d=true; Im=Lu; Jm=Mu;
      case 'xy_v',
        msknam='mask_v';
        is2d=true; Im=Lv; Jm=Mv;
      case 'xyz_rho',
        msknam='mask_rho';
        is3d=true; Im=Lr; Jm=Mr; Km=Nr;
      case 'xyz_u',
        msknam='mask_u';
        is3d=true; Im=Lu; Jm=Mu; Km=Nr;
      case 'xyz_v',    
        msknam='mask_v';
        is3d=true; Im=Lv; Jm=Mv; Km=Nr;
      case 'xyz_w',    
        msknam='mask_rho';
        is3d=true; Im=Lr; Jm=Mr; Km=Nw;
    end
  end
end

water=is2d | is3d;

% If water data only, read in Land/Sea mask.

got_mask=false;

if (water),
  [Vnames,nvars]=nc_Vname(Fname);
  got_mask=strmatch(msknam,Vnames,'exact');
  if isempty (got_mask),
    got_mask=true;
  end,
  if (got_mask),
    if (url),
      mask=read_nc_url(Fname,msknam,Tindex,ReplaceValue,PreserveType);
    else
      mask=read_nc(Fname,msknam,Tindex,ReplaceValue,PreserveType);
    end
  else
%   [fn,pth]=uigetfile(grid_file,'Enter grid NetCDF file...');
%   gname=[pth,fn];
    gname=input('Enter grid NetCDF file: ');
    if (url),
      mask=read_nc_url(gname,msknam,Tindex,ReplaceValue,PreserveType);
    else
      mask=read_nc(gname,msknam,Tindex,ReplaceValue,PreserveType);
    end
  end
end

%----------------------------------------------------------------------------
% Read in requested variable.
%----------------------------------------------------------------------------

if (water),

  if (nargin < 3),
    if (url),
      v=read_nc_url(Fname,Vname,Tindex,ReplaceValue,PreserveType);
    else
      v=read_nc(Fname,Vname);
    end
  else
    if (url),
      v=read_nc_url(Fname,Vname,Tindex,ReplaceValue,PreserveType);
    else
      v=read_nc(Fname,Vname,Tindex,ReplaceValue,PreserveType);
    end
  end
  [Npts,Ntime]=size(v);

  if (is2d),
    f=squeeze(ones([Im,Jm,Ntime])).*ReplaceValue;
    MASK=squeeze(repmat(mask,[1,1,Ntime]));
    ind=MASK > 0;
    f(ind)=v;
  elseif (is3d),    
    f=squeeze(ones([Im,Jm,Km,Ntime])).*ReplaceValue;
    MASK=squeeze(repmat(mask,[1,1,Km,Ntime]));
    ind=MASK > 0;
    f(ind)=v;
  end

else

  if (url),
    f=read_nc_url(Fname,Vname,Tindex,ReplaceValue,PreserveType);
  else
    f=read_nc(Fname,Vname,Tindex,ReplaceValue,PreserveType);
  end
  
end

return

function f=read_nc(Fname,Vname,Tindex,ReplaceValue,PreserveType)

%
% READ_NC:  Internal routine to read requested NetCDF variable
%
% f=read_nc(Fname,Vname,Tindex,ReplaceValue,PreserveType)
%
% This function reads in a generic multi-dimensional field from a NetCDF
% file.
%
%    Fname         NetCDF file name or URL name (character string)
%    Vname         NetCDF variable name to read (character string)
%    Tindex        Time record index to read (integer)
%    ReplaceValue  Value to use when "_FillValue" or "missing_value"
%                    attribute is found in variable.
%    PreserveType  Switch to preserve numerical data type. If false,
%                    convert numerical data to double precision. It
%                    has not effect data is already in double
%                    precision.
%
% On Output:
%
%    f             Field (scalar or array)
%

%  Set-up printing information switch.

global IPRINT

if (isempty(IPRINT)),
  IPRINT=1;
end,

%  Activate switch for reading specific record.

time_rec=0;
if (~isempty(Tindex)),
  time_rec=1;
end,

% Open NetCDF file.

[ncid]=mexnc('ncopen',Fname,'nc_nowrite');
if (ncid == -1),
  error(['READ_NC: ncopen - unable to open file: ' Fname])
  return
end,

% Supress all error messages from NetCDF.

[status]=mexnc('setopts',0);

% Define NetCDF parameters.

[NC_DOUBLE]=mexnc('parameter','nc_double');
[NC_FLOAT ]=mexnc('parameter','nc_float');
[NC_INT   ]=mexnc('parameter','nc_int');
[NC_SHORT ]=mexnc('parameter','nc_short');
[NC_BYTE  ]=mexnc('parameter','nc_byte');
[NC_CHAR  ]=mexnc('parameter','nc_char');

%----------------------------------------------------------------------------
% Inquire about requested variable.
%----------------------------------------------------------------------------

% Get variable ID.

[varid]=mexnc('ncvarid',ncid,Vname);
if (varid < 0),
  [status]=mexnc('ncclose',ncid);
  nc_inq(Fname);
  disp('  ');
  error(['READ_NC: ncvarid - cannot find variable: ',Vname])
  return
end,

% Inquire about unlimmited dimension.

[ndims,nvars,natts,recdim,status]=mexnc('ncinquire',ncid);
if (status == -1),
  error(['READ_NC: ncinquire - cannot inquire file: ',Fname])
end,

% Get information about requested variable.

[Vname,nctype,nvdims,dimids,nvatts,status]=mexnc('ncvarinq',ncid,varid);
if (status == -1),
  error(['READ_NC: ncvarinq - unable to inquire about variable: ',Vname])
end,

% Inquire about the _FillValue attribute.

got.add_offset=false;
got.FillValue=false;
got.missing_value=false;
got.scale_factor=false;

for i = 0:nvatts-1,
  [attnam,status]=mexnc('inq_attname',ncid,varid,i);
  if (status == -1)
    error(['READ_NC: inq_attname: error while inquiring attribute ' num2str(i)])
  end,
  lstr=length(attnam);
  [atype,status]=mexnc('inq_atttype',ncid,varid,attnam(1:lstr));
  if (status == -1)
    error(['READ_NC: inq_atttype: error while inquiring attribute ' num2str(i)])
  end,
  switch (attnam(1:lstr))
    case 'add_offset'
      if (atype == NC_DOUBLE),
        [offset,status]=mexnc('get_att_double',ncid,varid,attnam(1:lstr)); 
      elseif (atype == NC_FLOAT),
        [offset,status]=mexnc('get_att_float' ,ncid,varid,attnam(1:lstr));
      elseif (atype == NC_INT),
        [offset,status]=mexnc('get_att_int'   ,ncid,varid,attnam(1:lstr));
      elseif (atype == NC_SHORT),
        [offset,status]=mexnc('get_att_short' ,ncid,varid,attnam(1:lstr));
      elseif (atype == NC_BYTE),
        [offset,status]=mexnc('get_att_schar' ,ncid,varid,attnam(1:lstr));
      elseif (atype == NC_CHAR),
        [offset,status]=mexnc('get_att_text'  ,ncid,varid,attnam(1:lstr));
      else
        [offset,status]=mexnc('ncattget'      ,ncid,varid,attnam(1:lstr));
      end
      if (status == -1),
        error(['READ_NC: ncattget error while reading: ', attnam(1:lstr)])
      end
      got.add_offset=true;
    case {'_FillValue', '_fillvalue', 'missing_value'}
      if (atype == NC_DOUBLE),
        [spval,status]=mexnc('get_att_double',ncid,varid,attnam(1:lstr)); 
      elseif (atype == NC_FLOAT),
        [spval,status]=mexnc('get_att_float' ,ncid,varid,attnam(1:lstr));
      elseif (atype == NC_INT),
        [spval,status]=mexnc('get_att_int'   ,ncid,varid,attnam(1:lstr));
      elseif (atype == NC_SHORT),
        [spval,status]=mexnc('get_att_short' ,ncid,varid,attnam(1:lstr));
      elseif (atype == NC_BYTE),
        [spval,status]=mexnc('get_att_schar' ,ncid,varid,attnam(1:lstr));
      elseif (atype == NC_CHAR),
        [spval,status]=mexnc('get_att_text'  ,ncid,varid,attnam(1:lstr));
      else
        [spval,status]=mexnc('ncattget'      ,ncid,varid,attnam(1:lstr));
      end
      if (status == -1),
        error(['READ_NC: ncattget error while reading: ', attnam(1:lstr)])
      end
      if (strcmp(attnam(1:lstr),'missing_value')),
        got.missing_value=true;
      else
        got.FillValue=true;
      end
    case 'scale_factor'
      if (atype == NC_DOUBLE),
        [scale,status]=mexnc('get_att_double',ncid,varid,attnam(1:lstr)); 
      elseif (atype == NC_FLOAT),
        [scale,status]=mexnc('get_att_float' ,ncid,varid,attnam(1:lstr));
      elseif (atype == NC_INT),
        [scale,status]=mexnc('get_att_int'   ,ncid,varid,attnam(1:lstr));
      elseif (atype == NC_SHORT),
        [scale,status]=mexnc('get_att_short' ,ncid,varid,attnam(1:lstr));
      elseif (atype == NC_BYTE),
        [scale,status]=mexnc('get_att_schar' ,ncid,varid,attnam(1:lstr));
      elseif (atype == NC_CHAR),
        [scale,status]=mexnc('get_att_text'  ,ncid,varid,attnam(1:lstr));
      else
        [scale,status]=mexnc('ncattget'      ,ncid,varid,attnam(1:lstr));
      end
      if (status == -1),
        error(['READ_NC: ncattget error while reading: ', attnam(1:lstr)])
      end
      got.scale_factor=true;
  end
end

% Inquire about dimensions.

index=0;
for n=1:nvdims,
  [name,dsize,status]=mexnc('ncdiminq',ncid,dimids(n));
  if (status == -1),
    error(['READ_NC: ncdiminq - unable to inquire about dimension ID: ',...
          num2str(dimids(n))])
  else
    lstr=length(name);
    dimnam(n,1:lstr)=name(1:lstr);
    dimsiz(n)=dsize;
    start(n)=0;
    count(n)=dsize;
    if ((dimids(n) == recdim) || ~isempty(findstr(name,'time'))),
      index=n;
    end
  end
end

%  It reading specific time record, reset variable bounds.

nvdim=nvdims;
if (time_rec && (index > 0)),
  start(index)=Tindex-1;
  count(index)=1;
  nvdims=nvdims-1;
else
  index=0;
end

%----------------------------------------------------------------------------
% Read in requested variable.
%----------------------------------------------------------------------------

%  Read in scalar.

if (nvdim == 0),
  switch nctype
    case NC_DOUBLE
      [f,status]=mexnc('get_var_double',ncid,varid);
      myfunc='get_var_double';
    case NC_FLOAT
      [f,status]=mexnc('get_var_float' ,ncid,varid);
      myfunc='get_var_float';
    case NC_INT
      [f,status]=mexnc('get_var_int'   ,ncid,varid);
      myfunc='get_var_int';
    case NC_SHORT
      [f,status]=mexnc('get_var_short' ,ncid,varid);
      myfunc='get_var_short';
    case NC_BYTE
      [f,status]=mexnc('get_var_schar' ,ncid,varid);
      myfunc='get_var_schar';
    case NC_CHAR
      [f,status]=mexnc('get_var_text'  ,ncid,varid);
      myfunc='get_var_text';
    otherwise
      [f,status]=mexnc('ncvarget1'     ,ncid,varid,[0]);
      myfunc='ncvarget1';
  end
  if (status == -1),
    error(['READ_NC: ',myfunc,' - error while reading: ',Vname])
  end

%  Read in a multidemensional array.

else

  switch nctype
    case NC_DOUBLE
      [f,status]=mexnc('get_vara_double',ncid,varid,start,count);
      myfunc='get_vara_double';
    case NC_FLOAT
      [f,status]=mexnc('get_vara_float' ,ncid,varid,start,count);
      myfunc='get_vara_float';
    case NC_INT
      [f,status]=mexnc('get_vara_int'   ,ncid,varid,start,count);
      myfunc='get_vara_int';
    case NC_SHORT
      [f,status]=mexnc('get_vara_short' ,ncid,varid,start,count);
      myfunc='get_vara_short';
    case NC_BYTE
      [f,status]=mexnc('get_vara_schar' ,ncid,varid,start,count);
      myfunc='get_vara_schar';
    otherwise
      [f,status]=mexnc('ncvarget'       ,ncid,varid,start,count);
      myfunc='ncvarget';
  end
  if (status == -1),
    error(['READ_NC: ',myfunc,' - error while reading: ',Vname])
  end,

  if (nvdims == 3),
    if (length(start) == 3),
      f=reshape(f,[count(3),count(2),count(1)]);
    elseif (length(start) == 4),
      f=reshape(f,[count(4),count(3),count(2)]);
    end
  end
  
  if (nvdims == 4),
    if (length(start) == 4),
      f=reshape(f,[count(4),count(3),count(2),count(1)]);
    elseif (length(start) == 5),
      f=reshape(f,[count(5),count(4),count(3),count(2)]);
    end  
  end
  
end

% Print information about variable.

if (IPRINT),
  if (nvdims > 0),
    disp('  ')
    disp([Vname ' has the following dimensions (input order):']);
    disp('  ')
    for n=1:nvdim,
      s=['           '  int2str(n) ') ' dimnam(n,:) ' = ' int2str(dimsiz(n))];
      disp(s);
    end,
    disp('  ')
    disp([Vname ' loaded into an array of size:  [' int2str(size(f)) ']']);
    disp('  ')
  else
    disp('  ')
    disp([Vname ' is a scalar and has a value of ',num2str(f)]);
    disp('  ')
  end,
end,

%----------------------------------------------------------------------------
%  Post-process read data.
%----------------------------------------------------------------------------

%  Search for fill value or missing values.

ind=[];

if (got.FillValue || got.missing_value),
  if (iscellstr(f) || ischar(f)),
    ind=find(f == spval);
    f(ind)=spval;
  else
    ind=find(abs(f-spval) < 4*eps(double(spval)));
  end
end

%  Scale data and/or add offset value.

if (isnumeric(f)),
  if (got.add_offset),
    if (got.scale_factor),
      switch nctype
        case {NC_FLOAT, NC_DOUBLE}
          f=f.*scale+offset;
        case {NC_INT, NC_SHORT, NC_BYTE}
          f=double(f).*scale+offset;
      end
    else
      switch nctype
        case {NC_FLOAT, NC_DOUBLE}
          f=f+offset;
        case {NC_INT, NC_SHORT, NC_BYTE}
          f=double(f)+offset;
      end
    end
  elseif (got.scale_factor);
    switch nctype
      case {NC_FLOAT, NC_DOUBLE}
        f=f.*scale;
      case {NC_INT, NC_SHORT, NC_BYTE}
        f=double(f).*scale;
    end
  end
end

%  Set fill values or missing values with specified ReplaceValue.

if (~isempty(ind) && isnumeric(f)),
  f(ind)=ReplaceValue;
end

%----------------------------------------------------------------------------
%  If desired, convert data to double precision.
%----------------------------------------------------------------------------

if (~PreserveType) && ~(iscellstr(f) || ischar(f)),
  f=double(f);
end

%----------------------------------------------------------------------------
% Close NetCDF file.
%----------------------------------------------------------------------------

[status]=mexnc('ncclose',ncid);
if (status == -1),
  error(['READ_NC: ncclose - unable to close NetCDF file.'])
end

return

function f=read_nc_url(Fname,Vname,Tindex,ReplaceValue,PreserveType)

%
% READ_NC_URL:  Internal routine to read requested NetCDF variable
%
% f=read_nc_url(Fname,Vname,Tindex,ReplaceValue,PreserveType)
%
% This function reads in a generic multi-dimensional field from a URL
% NetCDF file (OpenDAT).  It uses the SCNTOOLS java interface.  The
% internal parameter PRESERVE_FVD is set to TRUE so the data is always
% processed in column-major order (Fortran, Matlab). The data is not
% transposed by the SCNTOOLS interface.
%
% On Input:
%
%    Fname         NetCDF file name or URL name (character string)
%    Vname         NetCDF variable name to read (character string)
%    Tindex        Time record index to read (integer)
%    ReplaceValue  Value to use when "_FillValue" or "missing_value"
%                    attribute is found in variable.
%    PreserveType  Switch to preserve numerical data type. If false,
%                    convert numerical data to double precision. It
%                    has not effect data is already in double
%                    precision.
%
% On Output:
%
%    f             Field (scalar or array)
%

%  Initialize.

got.add_offset=false;
got.FillValue=false;
got.missing_value=false;
got.RecDim=false;
got.scale_factor=false;
nvdims=0;

%  Clear and set persistent switch to process data in column-major order.

clear nc_getpref                          % clear defining function
setpref('SNCTOOLS','PRESERVE_FVD',true);

%  Activate switch for reading specific record.

time_rec=false;
if (~isempty(Tindex)),
  time_rec=true;
end,

%  Inquire information about requested variable in NetCDF file.

Info=nc_varinfo(Fname,Vname);

nvdims=length(Info.Dimension);
nvatts=length(Info.Attribute);
nctype=Info.Nctype;

%  Check if there is an unlimited dimension or a time dimension.

if (nvdims > 0),
  for n=1:nvdims,
    dimnam=Info.Dimension{n};
    lstr=length(dimnam);
    vdnames(n,1:lstr)=dimnam(1:lstr);
    vdsizes(n)=Info.Size(n);
    if (Info.Unlimited || ~isempty(findstr(dimnam,'time'))),
      got.RecDim=true;
      TimeDimName=dimnam(1:lstr);
    end
  end
end

%  Inquire information about the attributes.

if (nvatts > 0),
  for n=1:nvatts,
    attnam=Info.Attribute(n).Name;
    lstr=length(attnam);
    switch attnam
      case 'add_offset'
        offset=Info.Attribute(n).Value;
        got.add_offset=true;
      case {'_FillValue', '_fillvalue', 'missing_value'}
        spval=Info.Attribute(n).Value;
        if (strcmp(attnam(1:lstr),'missing_value')),
          got.missing_value=true;
        else
          got.FillValue=true;
        end
      case 'scale_factor'
        scale=Info.Attribute(n).Value;
        got.scale_factor=true;
    end
  end
end

%  Set start and count indices to process.

for n=1:nvdims,
  name=deblank(vdnames(n,:));
  lstr=length(name);
  start(n)=0;
% count(n)=vdsizes(n)-1;
  count(n)=Inf;
  if (time_rec && got.RecDim),
    if (strcmp(name(1:lstr),TimeDimName)),
      start(n)=Tindex-1;
      count(n)=1;
    end
  end
end

%----------------------------------------------------------------------------
% Read in requested variable. Use SNCTOOLS function "nc_vargetr" to read
% variable in its native (raw) data type.
%----------------------------------------------------------------------------

%  Read in scalar.

if (nvdims == 0),

  f=nc_vargetr(Fname,Vname);

%  Read in a multidemensional array.

else

  f=nc_vargetr(Fname,Vname,start,count);

end

%----------------------------------------------------------------------------
%  Post-process read data.
%----------------------------------------------------------------------------

%  Search for fill value or missing values.

ind=[];

if (got.FillValue || got.missing_value),
  if (iscellstr(f) || ischar(f)),
    ind=find(f == spval);
    f(ind)=spval;
  else
    ind=find(abs(f-spval) < 4*eps(double(spval)));
  end
end

%  Scale data and/or add offset value.

if (isnumeric(f)),
  if (got.add_offset),
    if (got.scale_factor),
      switch nctype
        case {NC_FLOAT, NC_DOUBLE}
          f=f.*scale+offset;
        case {NC_INT, NC_SHORT, NC_BYTE}
          f=double(f).*scale+offset;
      end
    else
      switch nctype
        case {NC_FLOAT, NC_DOUBLE}
          f=f+offset;
        case {NC_INT, NC_SHORT, NC_BYTE}
          f=double(f)+offset;
      end
    end
  elseif (got.scale_factor);
    switch nctype
      case {NC_FLOAT, NC_DOUBLE}
        f=f.*scale;
      case {NC_INT, NC_SHORT, NC_BYTE}
        f=double(f).*scale;
    end
  end
end

%  Set fill values or missing values with specified ReplaceValue.

if (~isempty(ind) && isnumeric(f)),
  f(ind)=ReplaceValue;
end

%----------------------------------------------------------------------------
%  Convert data to double precision.
%----------------------------------------------------------------------------

if (~PreserveType) && ~(iscellstr(f) || ischar(f)),
  f=double(f);
end

return
