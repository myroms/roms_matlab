function [vname,nvars]=nc_vname(fname);

%
% NC_VNAME:  Get the names of all variables in a NetCDF file
%
% [vname,nvars]=nc_vname(fname)
%
% This function gets the number and name of variables in requested
% NetCDF file.
%
% On Input:
%
%    fname      NetCDF file name or URL name (character string)
%
% On Output:
%
%    vname      Variable names (character array)
%    nvars      Number of variables
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2011 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Check if input file is URL from an OpenDAP server and process with the
%  appropriate interface.

url=nc_url(fname);

if (url),
  [vname,nvars]=nc_vname_java(fname);  
else
  [vname,nvars]=nc_vname_mexnc(fname);  
end

return

function [vname,nvars]=nc_vname_java(fname);

%
% NC_VNAME_JAVA:  Get the names of all variables in a NetCDF file
%
% [vname,nvars]=nc_vname_java(fname)
%
% This function gets the number and name of variables in requested
% URL OpenDAP NetCDF file(s).
%

%  Inquire information from URL NetCDF file.

Info=nc_info(fname); 

%  Extract requested variable information.

nvars=length(Info.Dataset);

for n=1:nvars,
  name=Info.Dataset(n).Name;
  lstr=length(name);
  vname(n,1:lstr)=name(1:lstr);
end

return

function [vname,nvars]=nc_vname_mexnc(fname);

%
% NC_VNAME:  Get the names of all variables in a NetCDF file
%
% [vname,nvars]=nc_vname_mexnc(fname)
%
% This function gets the number and name of variables in requested
% NetCDF file.  It cannot process URL OpenDAP file(s).

%  Open NetCDF file.
 
[ncid]=mexnc('ncopen',fname,'nc_nowrite');
if (ncid == -1),
  error(['NC_VNAME: ncopen - unable to open file: ', fname]);
  return
end
 
%  Supress all error messages from NetCDF.
 
[ncopts]=mexnc('setopts',0);

%  Inquire about the contents of NetCDf file. Display information.

[ndims,nvars,ngatts,recdim,status]=mexnc('ncinquire',ncid);
if (status == -1),
  error([ 'NC_VNAME: ncinquire - error while inquiring file: ' fname]);
end,

%  Inquire information about all the variables.

n=0;
for i=0:nvars-1,
  [name,vartyp,nvdims,vdims,nvatts,status]=mexnc('ncvarinq',ncid,i);
  if (status == -1),
    error(['NC_VNAME: ncvarinq - unable to inquire about variable ID: ',...
          num2str(i)]);
  else,
    n=n+1;
    lstr=length(name);
    vname(n,1:lstr)=name(1:lstr);
  end,
end,

%  Close NetCDF file.

[status]=mexnc('ncclose',ncid);
if (status == -1),
  error(['NC_VNAME: ncclose - unable to close NetCDF file.']);
end

return
