function cmap = gmt_topo(varargin)

% GMT_Accent: 256 color palette
%
% cmap = gmt_topo(M)
%
% MPL Accent colormap from NCL NCAR Graphics.
%
% https://www.ncl.ucar.edu/Document/Graphics/ColorTables/MPL_Accent.shtml
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
%    colormap(gmt_topo)
%    colormap(flipud(gmt_topo))
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

cmap = [[0.850000, 0.467500, 0.554459],
        [0.850000, 0.467500, 0.575377],
        [0.850000, 0.467500, 0.596295],
        [0.850000, 0.467500, 0.617213],
        [0.850000, 0.467500, 0.638131],
        [0.850000, 0.467500, 0.659049],
        [0.850000, 0.467500, 0.679967],
        [0.850000, 0.467500, 0.700885],
        [0.850000, 0.467500, 0.721803],
        [0.900000, 0.540000, 0.798337],
        [0.900000, 0.540000, 0.818103],
        [0.900000, 0.540000, 0.837870],
        [0.900000, 0.540000, 0.857636],
        [0.900000, 0.540000, 0.877403],
        [0.900000, 0.540000, 0.897169],
        [0.883064, 0.540000, 0.900000],
        [0.863297, 0.540000, 0.900000],
        [0.843531, 0.540000, 0.900000],
        [0.824625, 0.540000, 0.900000],
        [0.808875, 0.540000, 0.900000],
        [0.793125, 0.540000, 0.900000],
        [0.777375, 0.540000, 0.900000],
        [0.761625, 0.540000, 0.900000],
        [0.745875, 0.540000, 0.900000],
        [0.730125, 0.540000, 0.900000],
        [0.714375, 0.540000, 0.900000],
        [0.698625, 0.540000, 0.900000],
        [0.682819, 0.540000, 0.900000],
        [0.666281, 0.540000, 0.900000],
        [0.649744, 0.540000, 0.900000],
        [0.633206, 0.540000, 0.900000],
        [0.616669, 0.540000, 0.900000],
        [0.600131, 0.540000, 0.900000],
        [0.583594, 0.540000, 0.900000],
        [0.567056, 0.540000, 0.900000],
        [0.550519, 0.540000, 0.900000],
        [0.540000, 0.546019, 0.900000],
        [0.543047, 0.567937, 0.905078],
        [0.546328, 0.590490, 0.910547],
        [0.549609, 0.613272, 0.916016],
        [0.552891, 0.636285, 0.921484],
        [0.556172, 0.659527, 0.926953],
        [0.559453, 0.682999, 0.932422],
        [0.562734, 0.706701, 0.937891],
        [0.566016, 0.730632, 0.943359],
        [0.569297, 0.754793, 0.948828],
        [0.565918, 0.774463, 0.950000],
        [0.560723, 0.793377, 0.950000],
        [0.555527, 0.812859, 0.950000],
        [0.550332, 0.832910, 0.950000],
        [0.545137, 0.853529, 0.950000],
        [0.539941, 0.874716, 0.950000],
        [0.534746, 0.896471, 0.950000],
        [0.529551, 0.918795, 0.950000],
        [0.524355, 0.941687, 0.950000],
        [0.522500, 0.950000, 0.934971],
        [0.522500, 0.950000, 0.911592],
        [0.522500, 0.950000, 0.888213],
        [0.522500, 0.950000, 0.864834],
        [0.522500, 0.950000, 0.841455],
        [0.522500, 0.950000, 0.818076],
        [0.522500, 0.950000, 0.794697],
        [0.522500, 0.950000, 0.771318],
        [0.522500, 0.950000, 0.747939],
        [0.522500, 0.950000, 0.724093],
        [0.522500, 0.950000, 0.699779],
        [0.522500, 0.950000, 0.675465],
        [0.522500, 0.950000, 0.651151],
        [0.522500, 0.950000, 0.626837],
        [0.522500, 0.950000, 0.602523],
        [0.522500, 0.950000, 0.578209],
        [0.522500, 0.950000, 0.553895],
        [0.522500, 0.950000, 0.529580],
        [0.539066, 0.950000, 0.522500],
        [0.561509, 0.950000, 0.522500],
        [0.583953, 0.950000, 0.522500],
        [0.606397, 0.950000, 0.522500],
        [0.628841, 0.950000, 0.522500],
        [0.651284, 0.950000, 0.522500],
        [0.673728, 0.950000, 0.522500],
        [0.696172, 0.950000, 0.522500],
        [0.718616, 0.950000, 0.522500],
        [0.741260, 0.950000, 0.522500],
        [0.764639, 0.950000, 0.522500],
        [0.788018, 0.950000, 0.522500],
        [0.811396, 0.950000, 0.522500],
        [0.834775, 0.950000, 0.522500],
        [0.858154, 0.950000, 0.522500],
        [0.881533, 0.950000, 0.522500],
        [0.904912, 0.950000, 0.522500],
        [0.928291, 0.950000, 0.522500],
        [0.950000, 0.948330, 0.522500],
        [0.950000, 0.924951, 0.522500],
        [0.950000, 0.901572, 0.522500],
        [0.950000, 0.878193, 0.522500],
        [0.950000, 0.854814, 0.522500],
        [0.950000, 0.831436, 0.522500],
        [0.950000, 0.808057, 0.522500],
        [0.950000, 0.784678, 0.522500],
        [0.950000, 0.761299, 0.522500],
        [0.950000, 0.737920, 0.522500],
        [0.944922, 0.728309, 0.524505],
        [0.939453, 0.719968, 0.526608],
        [0.933984, 0.711842, 0.528650],
        [0.928516, 0.703929, 0.530632],
        [0.923047, 0.696226, 0.532555],
        [0.917578, 0.688729, 0.534418],
        [0.912109, 0.681437, 0.536221],
        [0.906641, 0.674345, 0.537964],
        [0.901172, 0.667452, 0.539647],
        [0.895703, 0.660753, 0.541271],
        [0.890234, 0.654247, 0.542834],
        [0.884766, 0.647929, 0.544338],
        [0.879297, 0.641798, 0.545782],
        [0.873828, 0.635850, 0.547167],
        [0.868359, 0.630082, 0.548491],
        [0.862891, 0.624492, 0.549756],
        [0.857422, 0.619076, 0.550961],
        [0.851953, 0.613832, 0.552106],
        [0.846671, 0.635003, 0.635003],
        [0.841180, 0.630885, 0.630885],
        [0.835690, 0.626767, 0.626767],
        [0.830199, 0.622649, 0.622649],
        [0.824708, 0.618531, 0.618531],
        [0.819217, 0.614413, 0.614413],
        [0.813727, 0.610295, 0.610295],
        [0.808236, 0.606177, 0.606177],
        [0.802745, 0.602059, 0.602059],
        [0.450520, 0.497821, 0.700000],
        [0.440853, 0.540090, 0.700000],
        [0.431186, 0.586098, 0.700000],
        [0.421519, 0.635843, 0.700000],
        [0.411934, 0.688860, 0.700000],
        [0.402363, 0.700000, 0.654541],
        [0.392793, 0.700000, 0.594278],
        [0.387793, 0.705078, 0.538057],
        [0.402832, 0.732422, 0.512060],
        [0.417871, 0.759766, 0.482564],
        [0.432910, 0.787109, 0.449569],
        [0.459688, 0.800000, 0.440000],
        [0.483312, 0.800000, 0.440000],
        [0.506937, 0.800000, 0.440000],
        [0.530563, 0.800000, 0.440000],
        [0.554188, 0.800000, 0.440000],
        [0.577812, 0.800000, 0.440000],
        [0.601437, 0.800000, 0.440000],
        [0.626966, 0.802344, 0.443170],
        [0.659280, 0.813281, 0.458106],
        [0.691247, 0.824219, 0.473282],
        [0.722828, 0.835156, 0.488697],
        [0.753983, 0.846094, 0.504351],
        [0.784675, 0.857031, 0.520245],
        [0.814862, 0.867969, 0.536378],
        [0.844507, 0.878906, 0.552750],
        [0.873569, 0.889844, 0.569361],
        [0.900195, 0.899888, 0.585567],
        [0.902930, 0.898397, 0.593518],
        [0.905664, 0.897050, 0.601506],
        [0.908398, 0.895848, 0.609532],
        [0.911133, 0.894793, 0.617595],
        [0.913867, 0.893887, 0.625696],
        [0.916602, 0.893130, 0.633834],
        [0.919336, 0.892524, 0.642009],
        [0.922070, 0.892071, 0.650222],
        [0.924805, 0.891773, 0.658472],
        [0.927539, 0.891631, 0.666759],
        [0.930273, 0.891646, 0.675084],
        [0.933008, 0.891820, 0.683446],
        [0.935742, 0.892155, 0.691846],
        [0.938477, 0.892652, 0.700283],
        [0.941211, 0.893313, 0.708758],
        [0.943945, 0.894139, 0.717269],
        [0.946680, 0.895131, 0.725819],
        [0.949414, 0.896292, 0.734405],
        [0.952148, 0.897622, 0.743029],
        [0.954883, 0.899124, 0.751690],
        [0.957617, 0.900799, 0.760389],
        [0.960352, 0.902648, 0.769125],
        [0.963086, 0.904672, 0.777899],
        [0.965820, 0.906875, 0.786710],
        [0.968555, 0.909256, 0.795558],
        [0.971289, 0.911818, 0.804444],
        [0.974023, 0.914562, 0.813367],
        [0.976758, 0.917489, 0.822327],
        [0.979492, 0.920601, 0.831325],
        [0.982227, 0.923901, 0.840360],
        [0.984961, 0.927388, 0.849433],
        [0.987695, 0.931065, 0.858543],
        [0.990430, 0.934933, 0.867690],
        [0.993164, 0.938994, 0.876875],
        [0.995898, 0.943249, 0.886097],
        [0.998633, 0.947701, 0.895356],
        [1.000000, 0.975001, 0.950362],
        [1.000000, 0.975013, 0.951144],
        [1.000000, 0.975037, 0.951926],
        [1.000000, 0.975073, 0.952707],
        [1.000000, 0.975122, 0.953489],
        [1.000000, 0.975182, 0.954271],
        [1.000000, 0.975255, 0.955052],
        [1.000000, 0.975340, 0.955834],
        [1.000000, 0.975438, 0.956616],
        [1.000000, 0.975547, 0.957398],
        [1.000000, 0.975669, 0.958179],
        [1.000000, 0.975803, 0.958961],
        [1.000000, 0.975949, 0.959743],
        [1.000000, 0.976108, 0.960524],
        [1.000000, 0.976278, 0.961306],
        [1.000000, 0.976461, 0.962088],
        [1.000000, 0.976656, 0.962869],
        [1.000000, 0.976864, 0.963651],
        [1.000000, 0.977083, 0.964433],
        [1.000000, 0.977315, 0.965214],
        [1.000000, 0.977559, 0.965996],
        [1.000000, 0.977815, 0.966778],
        [1.000000, 0.978083, 0.967560],
        [1.000000, 0.978364, 0.968341],
        [1.000000, 0.978657, 0.969123],
        [1.000000, 0.978962, 0.969905],
        [1.000000, 0.979279, 0.970686],
        [1.000000, 0.979609, 0.971468],
        [1.000000, 0.979951, 0.972250],
        [1.000000, 0.980304, 0.973031],
        [1.000000, 0.980671, 0.973813],
        [1.000000, 0.981049, 0.974595],
        [1.000000, 0.981440, 0.975377],
        [1.000000, 0.981843, 0.976158],
        [1.000000, 0.982258, 0.976940],
        [1.000000, 0.982685, 0.977722],
        [1.000000, 0.983124, 0.978503],
        [1.000000, 0.983576, 0.979285],
        [1.000000, 0.984040, 0.980067],
        [1.000000, 0.984516, 0.980848],
        [1.000000, 0.985005, 0.981630],
        [1.000000, 0.985505, 0.982412],
        [1.000000, 0.986018, 0.983194],
        [1.000000, 0.986543, 0.983975],
        [1.000000, 0.987080, 0.984757],
        [1.000000, 0.987630, 0.985539],
        [1.000000, 0.988192, 0.986320],
        [1.000000, 0.988766, 0.987102],
        [1.000000, 0.989352, 0.987884],
        [1.000000, 0.989950, 0.988665],
        [1.000000, 0.990561, 0.989447],
        [1.000000, 0.991184, 0.990229],
        [1.000000, 0.991819, 0.991010],
        [1.000000, 0.992466, 0.991792],
        [1.000000, 0.993125, 0.992574],
        [1.000000, 0.993797, 0.993356],
        [1.000000, 0.994481, 0.994137],
        [1.000000, 0.995177, 0.994919],
        [1.000000, 0.995886, 0.995701],
        [1.000000, 0.996606, 0.996482],
        [1.000000, 0.997339, 0.997264],
        [1.000000, 0.998084, 0.998046],
        [1.000000, 0.998841, 0.998827],
        [1.000000, 0.999611, 0.999609]];

% Inter,polate to, requeste],d number of colors.

P = size(cmap,1);

if (P ~= M)
  cmap = interp1(1:size(cmap,1), cmap, linspace(1,P,M), 'linear');
end

return