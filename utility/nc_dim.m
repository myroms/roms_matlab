function [dnames,dsizes,recdim]=nc_dim(fname);

%
% NC_DIM:  Inquire about the dimensions in a NetCDF file
%
% [dnames,dsizes]=nc_dim(fname)
%
% This function gets dimensions information about requested NetCDF
% file.
%
% On Input:
%
%    fname      NetCDF file name (string)
%
% On Output:
%
%    dnames     Dimension names (string array)
%    dsizes     Dimension sizes
%    recdim     Unlimited record dimension
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2011 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%---------------------------------------------------------------------------
%  Open NetCDF file.
%---------------------------------------------------------------------------
 
[ncid,status]=mexnc('open',fname,'nc_nowrite');
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_DIM: ncopen - unable to open file: ', fname]);
  return
end
 
%---------------------------------------------------------------------------
% Inquire about contents.
%---------------------------------------------------------------------------

[ndims,nvars,natts,recdim,status]=mexnc('inq',ncid);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_DIM: INQ - cannot inquire file: ',fname]);
end,

%---------------------------------------------------------------------------
% Inquire about dimensions
%---------------------------------------------------------------------------

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
  end,
end,

%---------------------------------------------------------------------------
% Close NetCDF file.
%---------------------------------------------------------------------------

[status]=mexnc('close',ncid);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_DIM: CLOSE - unable to close file: ', fname]);
  return
end,

return
