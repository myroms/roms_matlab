function [lon,lat,h]=x_etopo5(Llon, Rlon, Blat, Tlat, Fname);

%
% X_ETOPO5:  Extract requested bathymetry from ETOPO-5 dataset
%
% [lon,lat,h]=x_etopo5(Llon, Rlon, Blat, Tlat, Fname)
%
% This function extract bathymetry data from ETOPO-5 dataset in the
% specified geographical area.
%
%
% On Input:
%
%    Llon          Left   corner longitude (degrees_east)
%    Rlon          Right  corner longitude (degrees_east)
%    Blat          Bottom corner latitude (degrees_north)
%    Tlat          Top    corner latitude (degrees_north)
%
% On Output:
%
%    h             Bathymetry (m, matrix)
%    lon           Longitude (degree_east, matrix)
%    lat           Latitude (degree_north, matrix)
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2013 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%----------------------------------------------------------------------------
%  Read in ETOPO5 longitude and latitude.
%----------------------------------------------------------------------------

x=nc_read(Fname,'topo_lon');
y=nc_read(Fname,'topo_lat');

%----------------------------------------------------------------------------
%  Determine indices to extract.
%----------------------------------------------------------------------------

I=find(x >= Llon & x <= Rlon);
X=x(I);
[X,Isort]=sort(X);
I=I(Isort);

J=find(y >= Blat & y <= Tlat);
Jmin=min(J);
Jmax=max(J);
Y=y(J); Y=Y';

%----------------------------------------------------------------------------
%  Read in bathymetry.
%----------------------------------------------------------------------------

topo=nc_read(Fname,'topo');
h=topo(I,Jmin:Jmax);
[Im,Jm]=size(h);

%----------------------------------------------------------------------------
%  Set coordinates of extracted bathymetry.
%----------------------------------------------------------------------------

lon=repmat(X,[1 Jm]);
lat=repmat(Y,[Im 1]);

return
