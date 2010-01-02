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
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Activate switch for writing specific record.

time_rec=0;
if (nargin > 3),
  time_rec=1;
end,

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
[NC_FLOAT]=mexnc('parameter','nc_float');
[NC_INT]=mexnc('parameter','nc_int');
[NC_CHAR]=mexnc('parameter','nc_char');

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
end,

% Inquire about unlimmited dimension.

[ndims,nvars,natts,recdim,status]=mexnc('ncinquire',ncid);
if (status == -1),
  error(['NC_WRITE: ncinquire - cannot inquire file: ',fname])
end,

% Get information about requested variable.

[vname,nctype,nvdims,dimids,nvatts,status]=mexnc('ncvarinq',ncid,varid);
if (status == -1),
  error(['NC_WRITE: ncvarinq - unable to inquire about variable: ',vname])
end,

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
    end,
  end,
end,

% Inquire about the _FillValue attribute.

got_FillValue=0;

for i = 0:nvatts-1
  [attnam,status]=mexnc('ncattname',ncid,varid,i);
  if (status == -1)
    error(['NC_WRITE: ncattname: error while inquiring attribute ' ...
	    num2str(i)])
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
      error(['NC_WRITE: ncattget error while reading _FillValue attribute'])
    end,
    got_FillValue=1;
  end,
end,

%  It writing specific time record, reset variable bounds.  Check
%  for for files having the more than one datum in the unlimited
%  dimension (like variational data assimilation files).

datum=0;
if (unlimited & (nvdims == 1)),
  if (length(f) > 1)
    datum=1;
    if (time_rec),
      start(index)=tindex-1;
    else,
      start(index)=0;
    end,
    count(index)=length(f);
  end,
end,

if (~datum & time_rec & (index > 0)),
  start(index)=tindex-1;
  count(index)=1;
end,

%   Compute the minimum and maximum of the data to write.

if (time_rec),
  if (nvdims == 2),
    fmin=min(f);
    fmax=max(f);
  elseif (nvdims == 3),
    fmin=min(min(f));
    fmax=max(max(f));
  elseif (nvdims == 4),
    fmin=min(min(min(f)));
    fmax=max(max(max(f)));
  elseif (nvdims == 5),
    fmin=min(min(min(min(f))));
    fmax=max(max(max(max(f))));
  end,
else,
  if (nvdims == 1),
    fmin=min(f);
    fmax=max(f);
  elseif (nvdims == 2),
    fmin=min(min(f));
    fmax=max(max(f));
  elseif (nvdims == 3),
    fmin=min(min(min(f)));
    fmax=max(max(max(f)));
  elseif (nvdims == 4),
    fmin=min(min(min(min(f))));
    fmax=max(max(max(max(f))));
  end,
end,

%----------------------------------------------------------------------------
%  If _FillValue attribute, replace NaNs with fill value.
%----------------------------------------------------------------------------

if (got_FillValue),
  ind=isnan(f);
  if (~isempty(ind)),
    if (nctype == NC_FLOAT | nctype == NC_INT),
      f(ind)=double(spval);
    else
      f(ind)=spval;
    end,      
  end,
end,

%----------------------------------------------------------------------------
%  Write out variable into NetCDF file.
%----------------------------------------------------------------------------

if (nvdims > 0),
  [status]=mexnc('ncvarput',ncid,varid,start,count,f);
else,
  [status]=mexnc('ncvarput1',ncid,varid,0,f);
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
  error(['NC_WRITE: ncvarput - error while writting variable: ', vname])
end


%----------------------------------------------------------------------------
%  Close NetCDF file.
%----------------------------------------------------------------------------

[status]=mexnc('ncclose',ncid);
if (status == -1),
  error(['NC_WRITE: ncclose - unable to close NetCDF file: ', fname])
end

return
