function [status]=nc_write(fname,vname,f,tindex);

%
% NC_WRITE:  Writes a variable into a NetCDF file
%
% [status]=nc_write(fname,vname,f,tindex)
%
% This routine writes in a generic multi-dimensional field into a NetCDF
% file.
%
% On Input:
%
%    fname       NetCDF file name (character string)
%    vname       NetCDF variable name to read (character string)
%    f           Field (scalar, matrix or array)
%    tindex      Optional, time index to write (integer):
%                  *  If tindex is not provided as an argument during
%                     function call, it is assumed that entire variable
%                     is to be written. 
%                  *  If variable has an unlimitted record dimension,
%                     tindex can be used to increase that dimension or
%                     replace an already existing record.
%                  *  If variable has the word "time" in its dimension
%                     name, tindex can be use to write at the specified
%                     the time record.
%
% On Output:
%
%    status      Error flag
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Activate switch for writing specific record.

time_rec=0;
if (nargin > 3),
  time_rec=1;
end

%  Open NetCDF file.

[ncid]=mexnc('ncopen',fname,'nc_write');
if (ncid == -1),
  error(['NC_WRITE: ncopen - unable to open file: ', fname])
  return
end

%  Supress all error messages from NetCDF.

[ncopts]=mexnc('setopts',0);

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

[varid]=mexnc('ncvarid',ncid,vname);
if (varid < 0),
  [status]=mexnc('ncclose',ncid);
  nc_inq(fname);
  disp('  ');
  error(['NC_WRITE: ncvarid - cannot find variable: ',vname])
  return
end

% Inquire about unlimmited dimension.

[ndims,nvars,natts,recdim,status]=mexnc('ncinquire',ncid);
if (status == -1),
  error(['NC_WRITE: ncinquire - cannot inquire file: ',fname])
end,

% Get information about requested variable.

[vname,nctype,nvdims,dimids,nvatts,status]=mexnc('ncvarinq',ncid,varid);
if (status == -1),
  error(['NC_WRITE: ncvarinq - unable to inquire about variable: ',vname])
end

% Inquire about dimensions.

index=0;
unlimited=0;
for n=1:nvdims,
  [name,size,status]=mexnc('ncdiminq',ncid,dimids(n));
  if (status == -1),
    error(['NC_WRITE: ncdiminq - unable to inquire about dimension ID: ',...
          num2str(dimids(n))])
  else
    lstr=length(name);
    dimnam(n,1:lstr)=name(1:lstr);
    dimsiz(n)=size;
    start(n)=0;
    count(n)=size;
    if ((dimids(n) == recdim) | ~isempty(findstr(name,'time'))),
      unlimited=1;
      index=n;
    end
  end
end

% Inquire about the _FillValue attribute.

got_FillValue=0;

for i = 0:nvatts-1,
  [attnam,status]=mexnc('inq_attname',ncid,varid,i);
  if (status == -1)
    error(['NC_WRITE: inq_attname: error while inquiring attribute ' ...
	    num2str(i)])
  end
  lstr=length(attnam);
  [atype,status]=mexnc('inq_atttype',ncid,varid,attnam(1:lstr));
  if (status == -1)
    error(['READ_NC: inq_atttype: error while inquiring attribute ' num2str(i)])
  end,
  if (strcmp(attnam(1:lstr),'_FillValue')     ||   ...
      strcmp(attnam(1:lstr),'missing_value')),
    switch atype
      case NC_DOUBLE
        [spval,status]=mexnc('get_att_double',ncid,varid,attnam(1:lstr));
        myfunc='get_att_double';
      case NC_FLOAT
        [spval,status]=mexnc('get_att_float' ,ncid,varid,attnam(1:lstr));
        myfunc='get_att_float';
      case NC_INT
        [spval,status]=mexnc('get_att_int'   ,ncid,varid,attnam(1:lstr));
        myfunc='get_att_int';
      case NC_SHORT
        [spval,status]=mexnc('get_att_short' ,ncid,varid,attnam(1:lstr));
        myfunc='get_att_short';
      case NC_BYTE
        [spval,status]=mexnc('get_att_schar' ,ncid,varid,attnam(1:lstr));
        myfunc='get_att_schar';
      otherwise
        [spval,status]=mexnc('ncattget'      ,ncid,varid,attnam(1:lstr));
        myfunc='ncattget';
    end
    if (status == -1),
      error(['NC_WRITE: ',myfunc,' - error while reading ' ...
             attnam(1:lstr),' attribute'])
    end
    got_FillValue=1;
  end
end

%  It writing specific time record, reset variable bounds.  Check
%  for for files having the more than one datum in the unlimited
%  dimension (like variational data assimilation files).

datum=0;
if (unlimited & (nvdims == 1)),
  if (length(f) > 1)
    datum=1;
    if (time_rec),
      start(index)=tindex-1;
    else
      start(index)=0;
    end
    count(index)=length(f);
  end
end

if (~datum & time_rec & (index > 0)),
  start(index)=tindex-1;
  count(index)=1;
end

%   Compute the minimum and maximum of the data to write.

fmin=min(f(:));
fmax=max(f(:));

%----------------------------------------------------------------------------
%  If _FillValue attribute, replace NaNs with fill value.
%----------------------------------------------------------------------------

if (got_FillValue),
  ind=isnan(f);
  if (~isempty(ind))
    switch nctype
      case NC_DOUBLE
        f(ind)=double(spval);
      case NC_FLOAT
        f(ind)=single(spval);
      case NC_INT
        f(ind)=int32(spval);
      case NC_SHORT
        f(ind)=int16(spval);
      case NC_BYTE
        f(ind)=int8(spval);
      otherwise
        f(ind)=spval;
    end
  end
end

%----------------------------------------------------------------------------
%  Write out variable into NetCDF file.
%----------------------------------------------------------------------------

if (nvdims > 0),
  switch nctype
    case NC_DOUBLE
      [status]=mexnc('put_vara_double',ncid,varid,start,count,f);
      myfunc='put_vara_double';
    case NC_FLOAT
      [status]=mexnc('put_vara_float' ,ncid,varid,start,count,f);
      myfunc='put_vara_float';
    case NC_INT
      [status]=mexnc('put_vara_int'   ,ncid,varid,start,count,f);
      myfunc='put_vara_int';
    case NC_SHORT
      [status]=mexnc('put_vara_short' ,ncid,varid,start,count,f);
      myfunc='put_vara_short';
    case NC_BYTE
      [status]=mexnc('put_vara_schar' ,ncid,varid,start,count,f);
      myfunc='put_vara_schar';
    case NC_CHAR
      [status]=mexnc('put_vara_text'  ,ncid,varid,start,count,f);
      myfunc='put_var_text';
    otherwise
      [status]=mexnc('ncvarput'       ,ncid,varid,start,count,f);
  end
else,
  switch nctype
    case NC_DOUBLE
      [status]=mexnc('put_var_double',ncid,varid,f);
      myfunc='put_var_double';
    case NC_FLOAT
      [status]=mexnc('put_var_float' ,ncid,varid,f);
      myfunc='put_var_float';
    case NC_INT
      [status]=mexnc('put_var_int'   ,ncid,varid,f);
      myfunc='put_var_int';
    case NC_SHORT
      [status]=mexnc('put_var_short' ,ncid,varid,f);
      myfunc='put_var_short'; 
    case NC_BYTE
      [status]=mexnc('put_var_schar' ,ncid,varid,f);
      myfunc='put_var_schar';
    case NC_CHAR
      [status]=mexnc('put_var_text'  ,ncid,varid,f);
      myfunc='put_var_text';
   otherwise
      [status]=mexnc('ncvarput1'     ,ncid,varid,f);
      myfunc='ncvarput1';
  end
end,

if (status ~= -1 & nvdims > 1),
  text(1:19)=' ';
  text(1:length(vname))=vname;
  if (nargin > 3),
    disp(['Wrote ',sprintf('%19s',text), ...
          ' into record: ',num2str(tindex,'%4.4i'), ...
          ', Min=',sprintf('%12.5e',fmin),...
          ' Max=',sprintf('%12.5e',fmax)]);
  else
    disp(['Wrote ',sprintf('%19s',text), ...
          ' Min=',sprintf('%12.5e',fmin),...
          ' Max=',sprintf('%12.5e',fmax)]);
  end,
end,

if (status == -1),
  error(['NC_WRITE: ',myfunc,' - error while writting variable: ', vname])
end

%----------------------------------------------------------------------------
%  Close NetCDF file.
%----------------------------------------------------------------------------

[status]=mexnc('ncclose',ncid);
if (status == -1),
  error(['NC_WRITE: ncclose - unable to close NetCDF file: ', fname])
end

return
