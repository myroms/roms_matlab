function [dnames,dsizes,recdim]=nc_dim(fname)

%
% NC_DIM:  Inquire about the dimensions in a NetCDF file
%
% [dnames,dsizes,recdim]=nc_dim(fname)
%
% This function gets dimensions information about requested NetCDF
% file.
%
% On Input:
%
%    fname      NetCDF file name or URL file name (string)
%
% On Output:
%
%    dnames     Dimension names (string array)
%    dsizes     Dimension sizes
%    recdim     Unlimited record dimension ID (-1 if not found)
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Check if input file is URL from an OpenDAP server and process with the
%  appropriate interface.

url=nc_url(fname);

if (url),
  [dnames,dsizes,recdim]=nc_dim_java(fname);
else
  [dnames,dsizes,recdim]=nc_dim_mexnc(fname);
end

return

function [dnames,dsizes,recdim]=nc_dim_java(fname)

%
% NC_DIM_JAVA:  Inquire about the dimensions in a NetCDF file
%
% [dnames,dsizes,recdim]=nc_dim_java(fname)
%
% This function gets dimensions information about requested URL
% OpenDAP file(s). It uses SNCTOOLS function "nc_info".
%

%  Inquire information from URL NetCDF file.

Info=nc_info(fname); 

%  Extract dimension information.

ndims=length(Info.Dimension);
recdim=-1;
  
for n=1:ndims,
  lstr=length(Info.Dimension(n).Name);
  dnames(n,1:lstr)=Info.Dimension(n).Name;
  dsizes(n)=Info.Dimension(n).Length;
  if (Info.Dimension(n).Unlimited),
    recdim=n-1;                            % IDs start from zero
  end
end

return

function [dnames,dsizes,recdim]=nc_dim_mexnc(fname)

%
% NC_DIM_MEXNC:  Inquire about the dimensions in a NetCDF file
%
% [dnames,dsizes,recdim]=nc_dim_mexnc(fname)
%
% This function gets dimensions information about requested NetCDF
% file. It uses MEXNC functions. Therefore, it cannot process a URL
% OpenDAP file.
%
  
%  Open NetCDF file.

[ncid,status]=mexnc('open',fname,'nc_nowrite');
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_DIM: ncopen - unable to open file: ', fname]);
  return
end
 
% Inquire about contents.

[ndims,nvars,natts,recdim,status]=mexnc('inq',ncid);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_DIM: INQ - cannot inquire file: ',fname]);
end

% Inquire about dimensions

for n=1:ndims;
  [name,size,status]=mexnc('inq_dim',ncid,n-1);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['NC_DIM: INQ_DIM - unable to inquire about dimension ID: ',...
          num2str(n)]);
  else,
    lstr=length(name);
    dnames(n,1:lstr)=name(1:lstr);
    dsizes(n)=size;
  end
end

% Close NetCDF file.

[status]=mexnc('close',ncid);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_DIM: CLOSE - unable to close file: ', fname]);
  return
end

return
