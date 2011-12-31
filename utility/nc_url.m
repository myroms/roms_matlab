function [url]=nc_url(fname);

%
% NC_URL:  Detects if the NetCDF file is URL (OpenDAP server file)
%
% [url]=nc_url(fname)
%
% This function return true when the NetCDF dataset is from an URL
% indicating an OpenDAT server.
%
% On Input:
%
%    fname      NetCDF file name (string)
%
% On Output:
%
%    url        Switch (true/false) indicating a URL NetCDF file

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2011 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

if (~isempty(findstr(fname,'http:'))    || ...
    ~isempty(findstr(fname,'https:'))   || ...
    ~isempty(findstr(fname,'thredds')))
  url=true;
else
  url=false;
end

return
