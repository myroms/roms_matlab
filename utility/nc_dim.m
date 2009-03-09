function [dnames,dsizes]=nc_dim(fname);

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
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2009 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%---------------------------------------------------------------------------
%  Open NetCDF file.
%---------------------------------------------------------------------------
 
[ncid]=mexnc('ncopen',fname,'nc_nowrite');
if (ncid == -1),
  error(['NC_DIM: ncopen - unable to open file: ', fname]);
  return
end
 
%---------------------------------------------------------------------------
%  Supress all error messages from NetCDF.
%---------------------------------------------------------------------------
 
[ncopts]=mexnc('setopts',0);

%---------------------------------------------------------------------------
% Inquire about contents.
%---------------------------------------------------------------------------

[ndims,nvars,natts,recdim,status]=mexnc('ncinquire',ncid);
if (status == -1),
  error(['NC_DIM: ncinquire - cannot inquire file: ',fname]);
end,

%---------------------------------------------------------------------------
% Inquire about dimensions
%---------------------------------------------------------------------------

for n=1:ndims;
  [name,size,status]=mexnc('ncdiminq',ncid,n-1);
  if (status == -1),
    error(['NC_DIM: ncdiminq - unable to inquire about dimension ID: ',...
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

[status]=mexnc('ncclose',ncid);
if (status == -1),
  error(['NC_DIM: ncclose - unable to close file: ', fname]);
  return
end,

return
