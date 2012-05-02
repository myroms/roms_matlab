function [x,y,f] = nc_slice(Gfile, Hfile, Vname, depth, Tindex, zflag)

%
% NC_SLICE:  Interpolate requested slice from a 3D NetCDf variable
%
% [x,y,f] = nc_slice(Gfile, Hfile, Vname, depth, Tindex, zflag)
%
% This function computes a horizontal slice from a NetCDF file generated
% by ROMS. This function should only be used when the grid has variable
% bathymetry.  The field slice is interpolated at the requested depth
% from the input terrain-following coordinate data.
%
% On Input:
%
%    Gfile       Grid NetCDF file name (string):
%                  If no grid file, put field file name instead
%    Hfile       Field history NetCDF file name (string)
%    Vname       NetCDF variable name to process (string)
%    depth       Slice depth (scalar; meters, negative)
%    Tindex      Time index to read (integer):
%                  Tindex = 0  => All time records are processed
%                  Tindex > 0  => only "Tindex" record is processed
%    zflag       Variable depths computation switch (integer):
%                  zflag  = 0  => use zero free-surface
%                  zflag  = 1  => read in free-surface at "Tindex"
%
% On Output:
%
%    x           Slice X-positions (matrix)
%    y           Slice Y-positions (matrix)
%    f           Field slice (array)
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.
 
if (nargin < 5),
  Tindex=0;
end

if (nargin < 6),
  zflag=0;
end

%--------------------------------------------------------------------------
% Determine positions and Land/Sea masking variable names.
%--------------------------------------------------------------------------

[dnames,dsizes,igrid]=nc_vinfo(Hfile,Vname);

if (~isempty(strfind(dnames(1,:),'time'))),
  if ((length(dsizes)-1) < 3),
    error(['nc_slice - cannot vertically interpolate: ',Vname]);
  end
end

if (Tindex > 0)
  ndims=length(dsizes)-1;
else
  ndims=length(dsizes);
end

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
end

%--------------------------------------------------------------------------
% Read in variable positions and Land/Sea mask.
%--------------------------------------------------------------------------

x=nc_read(Gfile,Xname);
y=nc_read(Gfile,Yname);

mask=nc_read(Gfile,Mname);

%--------------------------------------------------------------------------
% Compute grid depths.
%--------------------------------------------------------------------------

switch ( zflag ),
  case 0
    z=depths(Hfile,Gfile,igrid,0,0);
  case 1
    z=depths(Hfile,Gfile,igrid,0,Tindex);
end,

[Im Jm Km]=size(z);

%--------------------------------------------------------------------------
% Compute requested slice.
%--------------------------------------------------------------------------

switch (ndims),
  case 3
    F=nc_read(Hfile,Vname,Tindex);
    for j=1:Jm,
      for i=1:Im,
        Zwrk=reshape(z(i,j,:),1,Km);
        Fwrk=reshape(F(i,j,:),1,Km);
        f(i,j)=interp1(Zwrk,Fwrk,depth);
      end
    end
  case 4
    Nm=dsizes(1);
    for n=1:Nm
      disp(['nc_slice - processing time record: ',num2str(n)]);
      F=nc_read(Hfile,Vname,n);
      for j=1:Jm,
        for i=1:Im,
          Zwrk=reshape(z(i,j,:),1,Km);
          Fwrk=reshape(F(i,j,:),1,Km);
          f(i,j,n)=interp1(Zwrk,Fwrk,depth);
        end
      end
    end
  otherwise
    error(['nc_slice - cannot vertically interpolate: ',Vname]);
end

return
