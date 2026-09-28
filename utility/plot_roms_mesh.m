function h = plot_roms_mesh(g,varargin)
% han = roms_plot_mesh(grd,[decimation_factor,color,cgridposition])
%
% Plot a mesh showing a ROMS grid over an existing plot
%
% Options can be in any order, or absent
%   decimation factor - integer > 0 controls density of mesh (default 10)
%   color - any 1-character color string supported by plot command
%         - real number < 1 is interpreted as grey scale
%         - vector assumed to be a RGB color
%   cgrid - string to specify which mesh to plot
%         'psi' (DEFAULT) or 'rho' mesh from respective points
%         'edge' or 'boundary' = perimeter of domain
%         'coast' = discrete land/sea boundary
%         'wet' =  discrete land/sea boundary from wet/dry mask - must
%                  have loaded wet/dry masks with roms_get_grid(f,f,tindex)
%         'line' = lon_coast,lat_coast in the grid structure
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu

if nargin==0
  help roms_plot_mesh
end

% get plot state
nextplotstatewas = get(gca,'nextplot');

% hold whatever is already plotted

set(gca,'nextplot','add')

n = 5;
color = 0.7*[1 1 1];
cgrid = 'psi';

for i=1:length(varargin)
  opt = varargin{i};
  if isnumeric(opt)
    if length(opt)>1
      color = opt;
    else
      if opt < 1
        color = opt*[1 1 1];
      else
        n = opt;
      end
    end
  end
  if ischar(opt)
    if length(opt)==1
      color = opt;
    else
      cgrid = opt(1:3);
    end
  end
end

switch cgrid(1)

  case 'r' % rho points
    x = g.lon_rho;
    y = g.lat_rho;
    han1=plot(x(1:n:end,1:n:end),y(1:n:end,1:n:end),'w-');
    han2=plot((x(1:n:end,1:n:end))',(y(1:n:end,1:n:end))','w-');
    han = [han1; han2];

  case 'p' % psi points DEFAULT
    % if you want to trim the plotted mesh then tinker with these masks
    % This is how I trimmed the mesh plotted in the DBOFS figure for ETOOFS
    m = ones(size(g.lon_psi));
    if 0
      m = g.mask_psi;
      m(m==0) = NaN;
    end
    x = m.*g.lon_psi;
    y = m.*g.lat_psi;
    han1=plot(x(1:n:end,1:n:end),y(1:n:end,1:n:end),'w-');
    han2=plot((x(1:n:end,1:n:end))',(y(1:n:end,1:n:end))','w-');
    % han1=plot(g.lon_psi(1:n:end,1:n:end),g.lat_psi(1:n:end,1:n:end),'w-');
    % han2=plot((g.lon_psi(1:n:end,1:n:end))',(g.lat_psi(1:n:end,1:n:end))','w-');
    han = [han1; han2];

  case {'c','w'} % discrete coast
    [Lp,Mp] = size(g.lon_psi);
    if strcmp(cgrid(1),'w')
      u = diff(g.wetdry_mask_rho,1,1);
      v = diff(g.wetdry_mask_rho,1,2);
    else
      u = diff(g.mask_rho,1,1);
      v = diff(g.mask_rho,1,2);
    end
    % trim to avoid bounds error and to avoid 'coast' coinciding with
    % points along perimeter of psi grid (trial and error - not certain
    % this is correct in all cases)
    [J,I] = find(v~=0);
    I(J==1) = [];
    J(J==1) = [];
    I(J>Mp-0) = [];
    J(J>Mp-0) = [];
    PXv(1,:) = diag(g.lon_psi(J-1,I));
    PXv(2,:) = diag(g.lon_psi(J,I));
    PYv(1,:) = diag(g.lat_psi(J-1,I));
    PYv(2,:) = diag(g.lat_psi(J,I));
    [J,I] = find(u~=0);
    J(I==1) = [];
    I(I==1) = [];
    I(J>Mp-1) = [];
    J(J>Mp-1) = [];
    J(I>Lp) = [];
    I(I>Lp) = [];
    PXu(1,:) = diag(g.lon_psi(J,I-1));
    PXu(2,:) = diag(g.lon_psi(J,I));
    PYu(1,:) = diag(g.lat_psi(J,I-1));
    PYu(2,:) = diag(g.lat_psi(J,I));
    PX = [PXv PXu];
    PY = [PYv PYu];
    han = plot(PX,PY,'k-');

  case 'l'
    han = plot(g.lon_coast,g.lat_coast,'k-');

  otherwise % we presume edge or boundary (i.e. perimeter)

    han = plot(g.perimeter(:,1),g.perimeter(:,2),'w-');

end

set(han,'linew',1,'color',color);

if nargout>0
  h = han;
end

% restore nextplotstate to what it was

set(gca,'nextplot',nextplotstatewas);

return
