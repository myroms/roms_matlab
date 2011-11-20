function [status]=nc_drename(fname,dname_old,dname_new);

%
% NC_DRENAME:  Changes the name of a dimension in a NetCDF file
%
% [status]=nc_drename(fname,dname_old,dname_new);
%
% This function renames requested dimension in a NetCDF file.
%
% On Input:
%
%    fname      NetCDF file name (character string)
%    dname_old  old dimension name (character string)
%    dname_new  new dimension name (character string)
%
% On Output:
%
%    status     error flag
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2011 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Initialize error flag.

status=-1;

%  Inquire about the names of all variables in the NetCDF.  Then, check
%  if the old and new variable names exist.

[Dnames,Dsizes]=nc_dim(fname);

got_old=strmatch(dname_old,Dnames,'exact');
if isempty (got_old),
  got_old=0;
end,
if (~got_old),
  disp(' ');
  disp([setstr(7),'*** Error:  NC_DRENAME - dimension not found = ', ...
	dname_old,setstr(7)]);
  disp(' ');
  return
end,

got_new=strmatch(dname_new,Dnames,'exact');
if isempty (got_new),
  got_new=0;
end,
if (got_new),
  disp(' ');
  disp([setstr(7),'*** Error:  NC_DRENAME - new name already exits = ', ...
	dname_new,setstr(7)]);
  disp(' ');
  return
end,

%  Open NetCDF file.

[ncid]=mexnc('open',fname,'nc_write');
if (ncid < 0),
  disp('  ');
  error(['NC_DRENAME: open - unable to open file: ', fname]);
  return
end

%  Get dimension ID.

[dimid]=mexnc('inq_dimid',ncid,dname_old);
if (dimid < 0),
  [status]=mexnc('close',ncid);
  disp('  ');
  error(['NC_VRENAME: ncvarid - cannot find variable: ',dname_old]);
  return
end,

%  Put open file into define mode.

[status]=mexnc('redef',ncid);
if (status < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  [status]=mexnc('close',ncid);
  error(['NC_DRENAME: redef - unable to put in definition mode: ',fname]);
  return
end,

%  Change name of requested variable.

[status]=mexnc('rename_dim',ncid,dimid,dname_new);
if (status < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_DRENAME: rename_dim - unable to rename dimension: ',dname_old]);
end,

%  Exit definition mode.

[status]=mexnc('enddef',ncid);
if (status < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_DRENAME: enddef - unable to exit definition mode: ',fname]);
  return
end,

%  Close NetCDF file.

[cstatus]=mexnc('ncclose',ncid);
if (cstatus < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_DRENAME: ncclose - unable to close NetCDF file: ', fname])
end

return
