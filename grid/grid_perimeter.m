function S = grid_perimeter(G)

%
% GRID_PERIMETER:  Sets Nested Grids Perimeters and Boundary Edges
%
% S = grid_perimeter(G)
%
% This function creates a structure containing information about nested
% grids perimeters, boundary edges, and other parameters.
%
% On Input:
%
%    G          Information grids structure (1 x Ngrid struct array)
%
%                 G(ng) = get_roms_grid ( char(Gnames(ng)) )
%
% On Output:
%
%    S          Nested grids information structure (struct array)
%
%
% This function adds the following fields to the output structure:
%
%    S.Ngrids                         - Number of nested grids
%    S.Ncontact                       - Number of contact regions
%    S.Nweights = 4                   - Number of horizontal weights
%    S.Ndatum = 0                     - Total number of contact points
%
%    S.western_edge  = 1              - Western  boundary edge index
%    S.southern_edge = 2              - Southern boundary edge index
%    S.eastern_edge  = 3              - Eastern  boundary edge index
%    S.northern_edge = 4              - Northern boundary edge index
%
%    S.spherical                      - Spherical switch
%
%    S.grid(ng).filename              - Grid NetCDF file name
%
%    S.grid(ng).Lp                    - Number of I-points (RHO)
%    S.grid(ng).Mp                    - Number of J-points (RHO)
%    S.grid(ng).L                     - Number of I-points (PSI)
%    S.grid(ng).M                     - Number of J-points (PSI)
%
%    S.grid(ng).refine_factor         - Refinement factor (0,3,5,7)
%
%    S.grid(ng).I_psi(:,:)            - ROMS I-indices at PSI-points
%    S.grid(ng).J_psi(:,:)            - ROMS J-indices at PSI-points
%    S.grid(ng).I_rho(:,:)            - ROMS I-indices at RHO-points
%    S.grid(ng).J_rho(:,:)            - ROMS J-indices at RHO-points
%    S.grid(ng).I_u  (:,:)            - ROMS I-indices at U-points
%    S.grid(ng).J_u  (:,:)            - ROMS J-indices at U-points
%    S.grid(ng).I_v  (:,:)            - ROMS I-indices at V-points
%    S.grid(ng).J_v  (:,:)            - ROMS I-indices at V-points
%
%    S.grid(ng).perimeter.X_psi(:)    - Perimeter X-coordinates (PSI)
%    S.grid(ng).perimeter.Y_psi(:)    - Perimeter Y-coordinates (PSI)
%    S.grid(ng).perimeter.X_rho(:)    - Perimeter X-coordinates (RHO)
%    S.grid(ng).perimeter.Y_rho(:)    - Perimeter Y-coordinates (RHO)
%    S.grid(ng).perimeter.X_u(:)      - Perimeter X-coordinates (U)
%    S.grid(ng).perimeter.Y_u(:)      - Perimeter Y-coordinates (U)
%    S.grid(ng).perimeter.X_v(:)      - Perimeter X-coordinates (V)
%    S.grid(ng).perimeter.Y_v(:)      - Perimeter Y-coordinates (V)
%                                       (counterclockwise)
%
%    S.grid(ng).corners.index(:)      - Corners linear IJ-index
%    S.grid(ng).corners.X(:)          - Corners X-coordinates (PSI)
%    S.grid(ng).corners.Y(:)          - Corners Y-coordinates (PSI)
%    S.grid(ng).corners.I(:)          - Corners I-indices (PSI)
%    S.grid(ng).corners.J(:)          - Corners J-indices (PSI)
%
%    S.grid(ng).boundary(ib).index(:) - Boundary linear IJ-index
%    S.grid(ng).boundary(ib).X(:)     - Boundary X-coordinates (PSI)
%    S.grid(ng).boundary(ib).Y(:)     - Boundary Y-coordinates (PSI)
%                                       (without corner points)
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.

S.Ngrids   = length(G);
S.Ncontact = (S.Ngrids-1)*2;
S.Nweights = 4;
S.Ndatum   = 0;

S.western_edge  = 1;
S.southern_edge = 2;
S.eastern_edge  = 3;
S.northern_edge = 4;

S.spherical = G(1).spherical;

for ng=1:S.Ngrids,
  S.grid(ng).filename = G(ng).grid_name;
  
  S.grid(ng).Lp = G(ng).Lm+2;
  S.grid(ng).Mp = G(ng).Mm+2;

  S.grid(ng).L  = G(ng).Lm+1;
  S.grid(ng).M  = G(ng).Mm+1;
end

% Get grid refinement factor, if any. Otherwise, set to zero.

for ng=1:S.Ngrids,
  S.grid(ng).refine_factor = 0;
  if (isfield(G(ng),'refine_factor')),
    if (~isempty(G(ng).refine_factor)),
      S.grid(ng).refine_factor = G(ng).refine_factor;
    end
  end
end

%--------------------------------------------------------------------------
% Spherical grids: set grid indices, perimeters, corners, and boundary
%                  edges.
%--------------------------------------------------------------------------

if (S.spherical),

  for ng=1:S.Ngrids,   
    Im = S.grid(ng).L;
    Jm = S.grid(ng).M;
    
    IstrP = 1;      IstrR = 1;        IstrU = 1;        IstrV = 1;
    IendP = Im;     IendR = Im+1;     IendU = Im;       IendV = Im+1;
    JstrP = 1;      JstrR = 1;        JstrU = 1;        JstrV = 1;
    JendP = Jm;     JendR = Jm+1;     JendU = Jm+1;     JendV = Jm;
    
% C-type variables ROMS indices. 

    [S.grid(ng).J_psi, S.grid(ng).I_psi] = meshgrid(1:Jm, 1:Im);
    [S.grid(ng).J_rho, S.grid(ng).I_rho] = meshgrid(0:Jm, 0:Im);
    [S.grid(ng).J_u  , S.grid(ng).I_u  ] = meshgrid(0:Jm, 1:Im);
    [S.grid(ng).J_v  , S.grid(ng).I_v  ] = meshgrid(1:Jm, 0:Im);

% Grid perimeter at PSI-points (counterclockwise from south). This is the
% physical grid perimeter.

    Xbox = [squeeze(G(ng).lon_psi(IstrP:IendP,JstrP));                  ...
            squeeze(G(ng).lon_psi(IendP,JstrP+1:JendP))';               ...
            squeeze(flipud(G(ng).lon_psi(IstrP:IendP-1,JendP)));        ...
            squeeze(fliplr(G(ng).lon_psi(IstrP,JstrP:JendP-1)))'];
    
    Ybox = [squeeze(G(ng).lat_psi(IstrP:IendP,JstrP));                  ...
            squeeze(G(ng).lat_psi(IendP,JstrP+1:JendP))';               ...
            squeeze(flipud(G(ng).lat_psi(IstrP:IendP-1,JendP)));        ...
            squeeze(fliplr(G(ng).lat_psi(IstrP,JstrP:JendP-1)))'];

    S.grid(ng).perimeter.X_psi = Xbox;
    S.grid(ng).perimeter.Y_psi = Ybox;

% Grid perimeter at RHO-points (counterclockwise from south). Needed for
% lateral boundary condition switch at RHO-points.

    Xbox = [squeeze(G(ng).lon_rho(IstrR:IendR,JstrR));                  ...
            squeeze(G(ng).lon_rho(IendR,JstrR+1:JendR))';               ...
            squeeze(flipud(G(ng).lon_rho(IstrR:IendR-1,JendR)));        ...
            squeeze(fliplr(G(ng).lon_rho(IstrR,JstrR:JendR-1)))'];
    
    Ybox = [squeeze(G(ng).lat_rho(IstrR:IendR,JstrR));                  ...
            squeeze(G(ng).lat_rho(IendR,JstrR+1:JendR))';               ...
            squeeze(flipud(G(ng).lat_rho(IstrR:IendR-1,JendR)));        ...
            squeeze(fliplr(G(ng).lat_rho(IstrR,JstrR:JendR-1)))'];

    S.grid(ng).perimeter.X_rho = Xbox;
    S.grid(ng).perimeter.Y_rho = Ybox;

% Grid perimeter at U-points (counterclockwise from south). Needed for
% lateral boundary condition switch at U-points.

    Xbox = [squeeze(G(ng).lon_u(IstrU:IendU,JstrU));                    ...
            squeeze(G(ng).lon_u(IendU,JstrU+1:JendU))';                 ...
            squeeze(flipud(G(ng).lon_u(IstrU:IendU-1,JendU)));          ...
            squeeze(fliplr(G(ng).lon_u(IstrU,JstrU:JendU-1)))'];
    
    Ybox = [squeeze(G(ng).lat_u(IstrU:IendU,JstrU));                    ...
            squeeze(G(ng).lat_u(IendU,JstrU+1:JendU))';                 ...
            squeeze(flipud(G(ng).lat_u(IstrU:IendU-1,JendU)));          ...
            squeeze(fliplr(G(ng).lat_u(IstrU,JstrU:JendU-1)))'];

    S.grid(ng).perimeter.X_u = Xbox;
    S.grid(ng).perimeter.Y_u = Ybox;

% Grid perimeter at V-points (counterclockwise from south). Needed for
% lateral boundary condition switch at V-points.

    Xbox = [squeeze(G(ng).lon_v(IstrV:IendV,JstrV));                    ...
            squeeze(G(ng).lon_v(IendV,JstrV+1:JendV))';                 ...
            squeeze(flipud(G(ng).lon_v(IstrV:IendV-1,JendV)));          ...
            squeeze(fliplr(G(ng).lon_v(IstrV,JstrV:JendV-1)))'];
    
    Ybox = [squeeze(G(ng).lat_v(IstrV:IendV,JstrV));                    ...
            squeeze(G(ng).lat_v(IendV,JstrV+1:JendV))';                 ...
            squeeze(flipud(G(ng).lat_v(IstrV:IendV-1,JendV)));          ...
            squeeze(fliplr(G(ng).lat_v(IstrV,JstrV:JendV-1)))'];

    S.grid(ng).perimeter.X_v = Xbox;
    S.grid(ng).perimeter.Y_v = Ybox;
    
% Grid domain corners (PSI-points).

    Icorners = false([Im Jm]);
    Icorners(1 ,1 ) = true;
    Icorners(Im,1 ) = true;
    Icorners(Im,Jm) = true;
    Icorners(1 ,Jm) = true;
    Icorners = find(Icorners(:) == true);

    S.grid(ng).corners.index = Icorners;
    S.grid(ng).corners.X     = G(ng).lon_psi(Icorners);
    S.grid(ng).corners.Y     = G(ng).lat_psi(Icorners);
    S.grid(ng).corners.I     = S.grid(ng).I_psi(Icorners);
    S.grid(ng).corners.J     = S.grid(ng).J_psi(Icorners);

% Boundary edges at PSI-points, excluding corners. 

    Iwest  = false([Im Jm]);    Iwest (1 ,2:Jm-1) = true;
    Isouth = false([Im Jm]);    Isouth(2:Im-1, 1) = true;
    Ieast  = false([Im Jm]);    Ieast (Im,2:Jm-1) = true;
    Inorth = false([Im Jm]);    Inorth(2:Im-1,Jm) = true;

    B(S.western_edge ).index = find(Iwest (:) == true);
    B(S.southern_edge).index = find(Isouth(:) == true);
    B(S.eastern_edge ).index = find(Ieast (:) == true);
    B(S.northern_edge).index = find(Inorth(:) == true);

    for ib=1:4,
      S.grid(ng).boundary(ib).index = B(ib).index;
      S.grid(ng).boundary(ib).X     = G(ng).lon_psi(B(ib).index);
      S.grid(ng).boundary(ib).Y     = G(ng).lat_psi(B(ib).index);
    end
  end
end

%--------------------------------------------------------------------------
% Cartesian grids: set grid indices, perimeters, corners, and boundary
%                  edges.
%--------------------------------------------------------------------------

if (~S.spherical),

  for ng=1:S.Ngrids,
    Im = S.grid(ng).L;
    Jm = S.grid(ng).M;
    
    IstrP = 1;      IstrR = 1;        IstrU = 1;        IstrV = 1;
    IendP = Im;     IendR = Im+1;     IendU = Im;       IendV = Im+1;
    JstrP = 1;      JstrR = 1;        JstrU = 1;        JstrV = 1;
    JendP = Jm;     JendR = Jm+1;     JendU = Jm+1;     JendV = Jm;
    
% C-type variables ROMS indices. 

    [S.grid(ng).J_psi, S.grid(ng).I_psi] = meshgrid(1:Jm, 1:Im);
    [S.grid(ng).J_rho, S.grid(ng).I_rho] = meshgrid(0:Jm, 0:Im);
    [S.grid(ng).J_u  , S.grid(ng).I_u  ] = meshgrid(0:Jm, 1:Im);
    [S.grid(ng).J_v  , S.grid(ng).I_v  ] = meshgrid(1:Jm, 0:Im);

% Grid perimeter at PSI-points (counterclockwise from south). This is the
% physical grid perimeter.

    Xbox = [squeeze(G(ng).x_psi(IstrP:IendP,JstrP));                    ...
            squeeze(G(ng).x_psi(IendP,JstrP+1:JendP))';                 ...
            squeeze(flipud(G(ng).x_psi(IstrP:IendP-1,JendP)));          ...
            squeeze(fliplr(G(ng).x_psi(IstrP,JstrP:JendP-1)))'];

    Ybox = [squeeze(G(ng).y_psi(IstrP:IendP,JstrP));                    ...
            squeeze(G(ng).y_psi(IendP,JstrP+1:JendP))';                 ...
            squeeze(flipud(G(ng).y_psi(IstrP:IendP-1,JendP)));          ...
            squeeze(fliplr(G(ng).y_psi(IstrP,JstrP:JendP-1)))'];

    S.grid(ng).perimeter.X_psi = Xbox;
    S.grid(ng).perimeter.Y_psi = Ybox;

% Grid perimeter at RHO-points (counterclockwise from south). Needed for
% lateral boundary condition switch at RHO-points.

    Xbox = [squeeze(G(ng).x_rho(IstrR:IendR,JstrR));                    ...
            squeeze(G(ng).x_rho(IendR,JstrR+1:JendR))';                 ...
            squeeze(flipud(G(ng).x_rho(IstrR:IendR-1,JendR)));          ...
            squeeze(fliplr(G(ng).x_rho(IstrR,JstrR:JendR-1)))'];
    
    Ybox = [squeeze(G(ng).y_rho(IstrR:IendR,JstrR));                    ...
            squeeze(G(ng).y_rho(IendR,JstrR+1:JendR))';                 ...
            squeeze(flipud(G(ng).y_rho(IstrR:IendR-1,JendR)));          ...
            squeeze(fliplr(G(ng).y_rho(IstrR,JstrR:JendR-1)))'];

    S.grid(ng).perimeter.X_rho = Xbox;
    S.grid(ng).perimeter.Y_rho = Ybox;

% Grid perimeter at U-points (counterclockwise from south). Needed for
% lateral boundary condition switch at U-points.

    Xbox = [squeeze(G(ng).x_u(IstrU:IendU,JstrU));                      ...
            squeeze(G(ng).x_u(IendU,JstrU+1:JendU))';                   ...
            squeeze(flipud(G(ng).x_u(IstrU:IendU-1,JendU)));            ...
            squeeze(fliplr(G(ng).x_u(IstrU,JstrU:JendU-1)))'];

    Ybox = [squeeze(G(ng).y_u(IstrU:IendU,JstrU));                      ...
            squeeze(G(ng).y_u(IendU,JstrU+1:JendU))';                   ...
            squeeze(flipud(G(ng).y_u(IstrU:IendU-1,JendU)));            ...
            squeeze(fliplr(G(ng).y_u(IstrU,JstrU:JendU-1)))'];
    
    S.grid(ng).perimeter.X_u = Xbox;
    S.grid(ng).perimeter.Y_u = Ybox;

% Grid perimeter at V-points (counterclockwise from south). Needed for
% lateral boundary condition switch at V-points.

    Xbox = [squeeze(G(ng).x_v(IstrV:IendV,JstrV));                      ...
            squeeze(G(ng).x_v(IendV,JstrV+1:JendV))';                   ...
            squeeze(flipud(G(ng).x_v(IstrV:IendV-1,JendV)));            ...
            squeeze(fliplr(G(ng).x_v(IstrV,JstrV:JendV-1)))'];
    
    Ybox = [squeeze(G(ng).y_v(IstrV:IendV,JstrV));                      ...
            squeeze(G(ng).y_v(IendV,JstrV+1:JendV))';                   ...
            squeeze(flipud(G(ng).y_v(IstrV:IendV-1,JendV)));            ...
            squeeze(fliplr(G(ng).y_v(IstrV,JstrV:JendV-1)))'];

    S.grid(ng).perimeter.X_v = Xbox;
    S.grid(ng).perimeter.Y_v = Ybox;

% Grid domain corners (PSI-points).

    Icorners = false([Im Jm]);
    Icorners(1 ,1 ) = true;
    Icorners(Im,1 ) = true;
    Icorners(Im,Jm) = true;
    Icorners(1 ,Jm) = true;
    Icorners = find(Icorners(:) == true);

    S.grid(ng).corners.index = Icorners;
    S.grid(ng).corners.X     = G(ng).x_psi(Icorners);
    S.grid(ng).corners.Y     = G(ng).y_psi(Icorners);
    S.grid(ng).corners.I     = S.grid(ng).I_psi(Icorners);
    S.grid(ng).corners.J     = S.grid(ng).J_psi(Icorners);

% Boundary edges at PSI-points, excluding corners.

    Iwest  = false([Im Jm]);    Iwest (1 ,2:Jm-1) = true;
    Isouth = false([Im Jm]);    Isouth(2:Im-1, 1) = true;
    Ieast  = false([Im Jm]);    Ieast (Im,2:Jm-1) = true;
    Inorth = false([Im Jm]);    Inorth(2:Im-1,Jm) = true;

    B(S.western_edge ).index = find(Iwest (:) == true);
    B(S.southern_edge).index = find(Isouth(:) == true);
    B(S.eastern_edge ).index = find(Ieast (:) == true);
    B(S.northern_edge).index = find(Inorth(:) == true);

    for ib=1:4,
      S.grid(ng).boundary(ib).index = B(ib).index;
      S.grid(ng).boundary(ib).X     = G(ng).x_psi(B(ib).index);
      S.grid(ng).boundary(ib).Y     = G(ng).y_psi(B(ib).index);
    end

  end
end

return
