function [Favg]=average(Fname,Vname,Tstr,Tend);

%
% AVERAGE:  Computes the time average of requested NetCDF variable
%
% [Favg]=average(Fname,Vname)
%
% This function reads in requested variable from input NetCDF file
% and computes its time average.
%
% On Input:
%
%    Fname       NetCDF file name (character string)
%    Vname       NetCDF variable name to process (character string)
%    Tstr        Starting time record to process (integer, OPTIONAL)
%    Tend        Ending   time record to process (integer, OPTIONAL)
%
% On Output:
%
%    Favg        Requested time-averaged field (array)
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Inquire number of time records.

[dnames,dsizes,igrid]=nc_vinfo(Fname,Vname);
ndims=length(dsizes);

for n=1:ndims,
  name=deblank(dnames(n,:));
  switch name
    case 'time',
      Nrec=dsizes(n);
  end,
end,

if (nargin < 3),
  Tstr=1;
  Tend=Nrec;
end,

% Read in field and compute its time mean.

Favg=nc_read(Fname,Vname,Tstr);

ic=1;

for n=Tstr+1:Tend,
  f=nc_read(Fname,Vname,n);
  Favg=Favg+f;
  ic=ic+1;
end,

Favg=Favg./ic;

return
