function cmap = mpl_sstanom(varargin)

% MPL_SSTANOM: 128 color palette
%
% cmap = mpl_sstanom(M)
%
% MPL SST anomaly colormap from NCL NCAR Graphics.
%
% https://www.ncl.ucar.edu/Document/Graphics/ColorTables/MPL_sstanom.shtml
%
% On Input:
%
%    M        Number of colors (integer, OPTIONAL)
%
% On Ouput:
%
%    cmap     Mx3 colormap matrix
%
% Usage:
%
%    colormap(mpl_sstanom)
%    colormap(flipud(mpl_sstanom))
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2024 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.md                            Hernan G. Arango        %
%===========================================================================%

% Initialize.

switch numel(varargin)
  case 0
    M = 128;
  case 1
    M = varargin{1};
end

% Set 128 colormap.

cmap = [[0.429066, 0.000000, 0.855040],
        [0.447982, 0.000000, 0.847474],
        [0.466897, 0.000000, 0.839908],
        [0.486305, 0.000000, 0.832834],
        [0.506482, 0.000000, 0.826528],
        [0.526659, 0.000000, 0.820223],
        [0.547543, 0.000000, 0.813210],
        [0.570242, 0.000000, 0.804383],
        [0.592941, 0.000000, 0.795556],
        [0.604291, 0.000000, 0.791142],
        [0.574856, 0.035433, 0.800046],
        [0.559093, 0.050565, 0.805090],
        [0.510188, 0.096886, 0.820838],
        [0.481184, 0.119585, 0.830927],
        [0.394171, 0.187682, 0.861192],
        [0.365167, 0.210381, 0.871280],
        [0.278155, 0.279677, 0.901546],
        [0.220146, 0.326336, 0.921722],
        [0.191142, 0.349666, 0.931811],
        [0.111557, 0.418593, 0.959954],
        [0.062376, 0.463991, 0.977609],
        [0.013195, 0.509389, 0.995263],
        [0.000000, 0.532088, 1.000000],
        [0.000000, 0.600185, 1.000000],
        [0.000000, 0.645583, 1.000000],
        [0.000000, 0.695686, 1.000000],
        [0.000000, 0.722168, 1.000000],
        [0.000000, 0.801615, 1.000000],
        [0.000000, 0.850704, 1.000000],
        [0.000000, 0.896101, 1.000000],
        [0.000000, 0.918800, 1.000000],
        [0.014717, 0.972687, 0.985283],
        [0.051288, 0.982776, 0.948712],
        [0.087858, 0.992864, 0.912142],
        [0.127013, 1.000000, 0.872987],
        [0.172411, 1.000000, 0.827589],
        [0.217809, 1.000000, 0.782191],
        [0.240508, 1.000000, 0.759492],
        [0.310096, 1.000000, 0.691396],
        [0.356755, 1.000000, 0.645998],
        [0.402860, 1.000000, 0.601799],
        [0.441953, 1.000000, 0.572795],
        [0.481046, 1.000000, 0.543791],
        [0.520138, 1.000000, 0.514787],
        [0.547082, 1.000000, 0.525875],
        [0.560323, 1.000000, 0.532180],
        [0.600046, 1.000000, 0.551096],
        [0.624375, 1.000000, 0.562630],
        [0.648335, 1.000000, 0.573979],
        [0.672295, 1.000000, 0.585329],
        [0.695317, 1.000000, 0.596678],
        [0.718016, 1.000000, 0.608028],
        [0.740715, 1.000000, 0.619377],
        [0.749189, 0.997124, 0.625052],
        [0.751080, 0.964967, 0.642076],
        [0.752341, 0.943529, 0.653426],
        [0.752280, 0.921430, 0.664775],
        [0.751019, 0.898731, 0.676125],
        [0.749758, 0.876032, 0.687474],
        [0.749020, 0.853333, 0.698824],
        [0.749020, 0.830634, 0.710173],
        [0.749020, 0.819285, 0.715848],
        [0.753249, 0.792157, 0.728258],
        [0.767120, 0.792157, 0.724475],
        [0.780992, 0.792157, 0.720692],
        [0.796586, 0.796586, 0.713956],
        [0.819285, 0.819285, 0.695040],
        [0.841984, 0.841984, 0.676125],
        [0.864683, 0.864683, 0.657316],
        [0.887382, 0.887382, 0.639662],
        [0.910081, 0.910081, 0.622007],
        [0.932780, 0.932780, 0.604352],
        [0.954248, 0.953018, 0.583007],
        [0.975686, 0.973195, 0.561569],
        [0.997124, 0.993372, 0.540131],
        [1.000000, 0.987774, 0.518800],
        [1.000000, 0.953725, 0.443137],
        [1.000000, 0.931027, 0.392695],
        [1.000000, 0.908328, 0.343207],
        [1.000000, 0.885629, 0.294025],
        [1.000000, 0.862930, 0.244844],
        [1.000000, 0.837785, 0.189143],
        [1.000000, 0.811303, 0.129873],
        [1.000000, 0.784821, 0.070604],
        [1.000000, 0.760369, 0.035694],
        [1.000000, 0.737670, 0.021822],
        [1.000000, 0.714971, 0.007951],
        [1.000000, 0.692272, 0.000000],
        [1.000000, 0.669573, 0.000000],
        [1.000000, 0.646874, 0.000000],
        [1.000000, 0.623775, 0.000000],
        [1.000000, 0.611795, 0.000000],
        [1.000000, 0.575855, 0.000000],
        [1.000000, 0.551111, 0.000000],
        [1.000000, 0.523368, 0.000000],
        [1.000000, 0.495625, 0.000000],
        [1.000000, 0.466159, 0.000000],
        [1.000000, 0.420761, 0.000000],
        [1.000000, 0.375363, 0.000000],
        [1.000000, 0.329965, 0.000000],
        [1.000000, 0.275848, 0.000000],
        [1.000000, 0.221622, 0.000000],
        [1.000000, 0.167397, 0.000000],
        [1.000000, 0.120923, 0.000000],
        [1.000000, 0.075525, 0.000000],
        [1.000000, 0.030127, 0.000000],
        [0.991280, 0.014764, 0.035848],
        [0.985606, 0.011611, 0.059177],
        [0.968581, 0.002153, 0.129166],
        [0.956401, 0.000000, 0.179977],
        [0.943791, 0.000000, 0.232941],
        [0.931180, 0.000000, 0.285905],
        [0.919262, 0.000000, 0.339562],
        [0.907912, 0.000000, 0.393787],
        [0.896563, 0.000000, 0.448012],
        [0.881338, 0.000000, 0.485629],
        [0.861161, 0.000000, 0.502022],
        [0.840984, 0.000000, 0.518416],
        [0.814579, 0.000000, 0.508651],
        [0.775486, 0.000000, 0.445598],
        [0.736394, 0.000000, 0.382545],
        [0.698685, 0.000000, 0.321707],
        [0.665898, 0.000000, 0.268743],
        [0.649504, 0.000000, 0.242261],
        [0.600323, 0.000000, 0.162676],
        [0.567536, 0.000000, 0.108451],
        [0.534748, 0.000000, 0.054225],
        [0.501961, 0.000000, 0.000000]];

% Interpolate to requested number of colors.

P = size(cmap,1);

if (P ~= M)
  cmap = interp1(1:size(cmap,1), cmap, linspace(1,P,M), 'linear');
end

return