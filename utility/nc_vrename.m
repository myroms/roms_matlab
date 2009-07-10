function [status]=nc_vrename(fname,vname_old,vname_new);

%
% NC_VRENAME:  Changes the name of a variable in a NetCDF file
%
% [status]=nc_vrename(fname,vname_old,vname_new);
%
% This function renames requested variable in a NetCDF file.
%
% On Input:
%
%    fname      NetCDF file name (character string)
%    vname_old  old variable name (character string)
%    vname_new  new variable name (character string)
%
% On Output:
%
%    status     error flag
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2009 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Initialize error flag.

status=-1;

%  Inquire about the names of all variables in the NetCDF.  Then, check
%  if the old and new variable names exist.

[Vnames,nvars]=nc_vname(fname);

got_old=strmatch(vname_old,Vnames,'exact');
if isempty (got_old),
  got_old=0;
end,
if (~got_old),
  disp(' ');
  disp([setstr(7),'*** Error:  NC_VRENAME - variable not found = ', ...
	vname_old,setstr(7)]);
  disp(' ');
  return
end,

got_new=strmatch(vname_new,Vnames,'exact');
if isempty (got_new),
  got_new=0;
end,
if (got_new),
  disp(' ');
  disp([setstr(7),'*** Error:  NC_VRENAME - new name already exits = ', ...
	vname_new,setstr(7)]);
  disp(' ');
  return
end,

%  Open NetCDF file.

[ncid]=mexnc('open',fname,'nc_write');
if (ncid < 0),
  disp('  ');
  error(['NC_VRENAME: open - unable to open file: ', fname]);
  return
end

%  Get variable ID.

[varid]=mexnc('ncvarid',ncid,vname_old);
if (varid < 0),
  [status]=mexnc('close',ncid);
  disp('  ');
  error(['NC_VRENAME: ncvarid - cannot find variable: ',vname_old]);
  return
end,

%  Change name of requested variable.

[status]=mexnc('rename_var',ncid,varid,vname_new);
if (status < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_VRENAME: rename_var - unable to rename variable: ',vname_old]);
end,

%  Close NetCDF file.

[cstatus]=mexnc('ncclose',ncid);
if (cstatus < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_VRENAME: ncclose - unable to close NetCDF file: ', fname])
end

return
