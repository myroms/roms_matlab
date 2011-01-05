function [x,y,f]=nc_slice(gname, fname, vname, depth, tindex, zflag);

%
% NC_SLICE:  Interpolate requested slice from a 3D NetCDf variable
%
% [x,y,f]=nc_slice(gname, fname, vname, depth, tindex, zflag)
%
% This function computes a horizontal slice from a NetCDF file generated
% by ROMS. This function should only be used when the grid has variable
% bathymetry.  The field slice is interpolated at the requested depth
% from the input terrain-following coordinate data.
%
% On Input:
%
%    gname       Grid NetCDF file name (character string).
%                  If no grid file, put field file name instead.
%    fname       Field NetCDF file name (character string).
%    vname       NetCDF variable name to process (character string).
%    depth       Slice depth (scalar; meters, negative).
%    tindex      Time index to read (integer).
%                  tindex=0  => All time records are processed.
%                  tindex>0  => only "tindex" record is processed.
%    zflag       Variable depths computation switch (integer):
%                  zflag=0  => use zero free-surface.
%                  zflag=1  => read in free-surface at "tindex".
%
% On Output:
%
%    x           Slice X-positions (matrix).
%    y           Slice Y-positions (matrix).
%    f           Field slice (array).
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2011 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Deactivate printing information switch from "nc_read".
 
global IPRINT

IPRINT=0;

if (nargin < 5),
  tindex=0;
end,

if (nargin < 6),
  zflag=0;
end,

%----------------------------------------------------------------------------
% Determine positions and Land/Sea masking variable names.
%----------------------------------------------------------------------------

[dnames,dsizes,igrid]=nc_vinfo(fname,vname);

if (~isempty(findstr(dnames(1,:),'time'))),
  if ((length(dsizes)-1) < 3),
    error(['nc_slice - cannot vertically interpolate: ',vname]);
  end
end,

if (tindex > 0)
  ndims=length(dsizes)-1;
else
  ndims=length(dsizes);
end,

switch ( igrid ),
  case 1
    Xname='lon_rho';
    Yname='lat_rho';
    Mname='mask_rho';
  case 2
    Xname='lon_psi';
    Yname='lat_psi';
    Mname='mask_psi';
  case 3
    Xname='lon_u';
    Yname='lat_u';
    Mname='mask_u';
  case 4
    Xname='lon_v';
    Yname='lat_v';
    Mname='mask_v';
  case 5
    Xname='lon_rho';
    Yname='lat_rho';
    Mname='mask_rho';
end,

%----------------------------------------------------------------------------
% Read in variable positions and Land/Sea mask.
%----------------------------------------------------------------------------

x=nc_read(gname,Xname);
y=nc_read(gname,Yname);

mask=nc_read(gname,Mname);

%----------------------------------------------------------------------------
% Compute grid depths.
%----------------------------------------------------------------------------

switch ( zflag ),
  case 0
    z=depths(fname,gname,igrid,0,0);
  case 1
    z=depths(fname,gname,igrid,0,tindex);
end,

[Im Jm Km]=size(z);

%----------------------------------------------------------------------------
% Compute requested slice.
%----------------------------------------------------------------------------

switch (ndims),
  case 3
    F=nc_read(fname,vname,tindex);
    for j=1:Jm,
      for i=1:Im,
        Zwrk=reshape(z(i,j,:),1,Km);
        Fwrk=reshape(F(i,j,:),1,Km);
        f(i,j)=interp1(Zwrk,Fwrk,depth);
      end,
    end,
  case 4
    Nm=dsizes(1);
    for n=1:Nm
      disp(['nc_slice - processing time record: ',num2str(n)]);
      F=nc_read(fname,vname,n);
      for j=1:Jm,
        for i=1:Im,
          Zwrk=reshape(z(i,j,:),1,Km);
          Fwrk=reshape(F(i,j,:),1,Km);
          f(i,j,n)=interp1(Zwrk,Fwrk,depth);
        end,
      end,
    end,
  otherwise
    error(['nc_slice - cannot vertically interpolate: ',vname]);
end,

return
