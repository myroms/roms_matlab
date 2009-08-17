function [Fvar]=variance(Fname,Vname,Favg,Tstr,Tend);

%
% VARIANCE:  Computes the variance of requested NetCDF variable
%
% [Fvar]=variance(Fname,Vname,Favg)
%
% This function computes the variance of requested NetCDF variable from
% its specified time mean.
%
% On Input:
%
%    Fname       NetCDF file name (character string)
%    Vname       NetCDF variable name to process (character string)
%    Favg        Variable time mean (array)
%    Tstr        Starting time record to process (integer, OPTIONAL)
%    Tend        Ending   time record to process (integer, OPTIONAL)
%
% On Output:
%
%    Fvar        Requested variable variance (squared field units, array)
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2009 The ROMS/TOMS Group                              %
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

if (nargin < 4),
  Tstr=1;
  Tend=Nrec;
end,

% Read in field and compute the variance from its time mean.

Fvar=nc_read(Fname,Vname,Tstr);
Fvar=zeros(size(Fvar));

ic=1;

for n=Tstr:Tend,
  f=nc_read(Fname,Vname,n);
  Fvar=Fvar+(f-Favg).^2;
  ic=ic+1;
end,

Fvar=Fvar./(ic-1);

return
