function cmap = mpl_rainbow200(varargin)

% MPL_Accent: 200 color palette
%
% cmap = mpl_rainbow200(M)
%
% MPL rainbow200 colormap from NCL NCAR Graphics.
%
% https://www.ncl.ucar.edu/Document/Graphics/ColorTables/BkBlAqGrYeOrReViWh200.shtml
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
%    colormap(mpl_rainbow200)
%    colormap(flipud(mpl_rainbow200))
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

% Set 200 colormap.

cmap = [[0.000000000000000, 0.000000000000000, 0.000000000000000],
        [0.000000000000000, 0.000000000000000, 0.061448201327775],
        [0.000000000000000, 0.000000000000000, 0.126817971283001],
        [0.000000000000000, 0.000000000000000, 0.188266172610777],
        [0.000000000000000, 0.000000000000000, 0.249714373938552],
        [0.000000000000000, 0.000000000000000, 0.376532345221553],
        [0.000000000000000, 0.000000000000000, 0.437980546549328],
        [0.000000000000000, 0.000000000000000, 0.501528485409912],
        [0.000000000000000, 0.000000000000000, 0.564798517832330],
        [0.000000000000000, 0.000000000000000, 0.626246719160105],
        [0.000000000000000, 0.000000000000000, 0.688621275281766],
        [0.000000000000000, 0.000000000000000, 0.753064690443106],
        [0.000000000000000, 0.000000000000000, 0.814512891770882],
        [0.000000000000000, 0.000000000000000, 0.875961093098657],
        [0.000000000000000, 0.000000000000000, 0.941330863053883],
        [0.000000000000000, 0.002779064381658, 1.000000000000000],
        [0.000000000000000, 0.064227265709433, 1.000000000000000],
        [0.000000000000000, 0.129597035664660, 1.000000000000000],
        [0.000000000000000, 0.191045236992435, 1.000000000000000],
        [0.000000000000000, 0.252493438320210, 1.000000000000000],
        [0.000000000000000, 0.317492666357882, 1.000000000000000],
        [0.000000000000000, 0.379311409603211, 1.000000000000000],
        [0.000000000000000, 0.440759610930987, 1.000000000000000],
        [0.000000000000000, 0.504585456229736, 1.000000000000000],
        [0.000000000000000, 0.567577582213988, 1.000000000000000],
        [0.000000000000000, 0.629025783541763, 1.000000000000000],
        [0.000000000000000, 0.691678246101590, 1.000000000000000],
        [0.000000000000000, 0.755843754824765, 1.000000000000000],
        [0.000000000000000, 0.817291956152540, 1.000000000000000],
        [0.000000000000000, 0.878771035973444, 1.000000000000000],
        [0.000000000000000, 0.944109927435541, 1.000000000000000],
        [0.000000000000000, 0.998332561371005, 0.994441871236683],
        [0.002779064381658, 0.979898100972672, 0.935772734290567],
        [0.003921568627451, 0.961463640574340, 0.876547784468118],
        [0.007225567392311, 0.943029180176007, 0.817940404508260],
        [0.007843137254902, 0.924594719777675, 0.758097884823221],
        [0.011672070402964, 0.906160259379342, 0.698965570480161],
        [0.011764705882353, 0.887725798981010, 0.641438937779836],
        [0.015686274509804, 0.869291338582677, 0.580423035355875],
        [0.015686274509804, 0.848201327775205, 0.522464103751737],
        [0.016643507796820, 0.828500849158561, 0.464937471051413],
        [0.019607843137255, 0.810066388760228, 0.403489269723637],
        [0.021090010807473, 0.791631928361896, 0.345962637023313],
        [0.023529411764706, 0.773197467963563, 0.288219854871082],
        [0.025536513818126, 0.754763007565231, 0.226987802995214],
        [0.027450980392157, 0.736328547166898, 0.169461170294890],
        [0.029983016828779, 0.717894086768566, 0.110545005403736],
        [0.039709742164582, 0.704461942257218, 0.057989810097267],
        [0.101157943492357, 0.722896402655550, 0.054901960784314],
        [0.162606144820133, 0.741330863053883, 0.050980392156863],
        [0.220472440944882, 0.759765323452215, 0.047398486953837],
        [0.281580978848232, 0.778199783850548, 0.043137254901961],
        [0.343029180176008, 0.796634244248881, 0.039215686274510],
        [0.404477381503783, 0.815068704647213, 0.036807163810406],
        [0.465925582831558, 0.833503165045546, 0.031372549019608],
        [0.524440327312027, 0.854871082291184, 0.028439092172302],
        [0.584900416859657, 0.874293654469662, 0.026215840666975],
        [0.646348618187432, 0.892728114867994, 0.020071020534198],
        [0.707796819515208, 0.911162575266327, 0.017847769028871],
        [0.769245020842983, 0.929597035664660, 0.015624517523545],
        [0.828408213679172, 0.948031496062992, 0.009479697390767],
        [0.888219854871082, 0.966465956461325, 0.007256445885441],
        [0.949668056198857, 0.984900416859657, 0.003921568627451],
        [1.000000000000000, 0.995553496989347, 0.000000000000000],
        [1.000000000000000, 0.974309093716227, 0.000000000000000],
        [1.000000000000000, 0.950316504554578, 0.000000000000000],
        [1.000000000000000, 0.925737224023468, 0.000000000000000],
        [1.000000000000000, 0.905079512119809, 0.000000000000000],
        [1.000000000000000, 0.880500231588698, 0.000000000000000],
        [1.000000000000000, 0.858607379959858, 0.000000000000000],
        [1.000000000000000, 0.835263239153929, 0.000000000000000],
        [1.000000000000000, 0.810683958622819, 0.000000000000000],
        [1.000000000000000, 0.790026246719160, 0.000000000000000],
        [1.000000000000000, 0.765446966188050, 0.000000000000000],
        [1.000000000000000, 0.742905666203489, 0.000000000000000],
        [1.000000000000000, 0.720209973753281, 0.000000000000000],
        [1.000000000000000, 0.695630693222171, 0.000000000000000],
        [1.000000000000000, 0.674972981318511, 0.000000000000000],
        [1.000000000000000, 0.650393700787401, 0.000000000000000],
        [1.000000000000000, 0.611270649992280, 0.000000000000000],
        [1.000000000000000, 0.570788945499459, 0.000000000000000],
        [1.000000000000000, 0.531696773197468, 0.000000000000000],
        [1.000000000000000, 0.488899181719932, 0.000000000000000],
        [1.000000000000000, 0.449590859966034, 0.000000000000000],
        [1.000000000000000, 0.410498687664042, 0.000000000000000],
        [1.000000000000000, 0.367484946734599, 0.000000000000000],
        [1.000000000000000, 0.328392774432607, 0.000000000000000],
        [1.000000000000000, 0.288868303226803, 0.000000000000000],
        [1.000000000000000, 0.246286861201173, 0.000000000000000],
        [1.000000000000000, 0.207194688899182, 0.000000000000000],
        [1.000000000000000, 0.166496834954454, 0.000000000000000],
        [1.000000000000000, 0.125088775667748, 0.000000000000000],
        [1.000000000000000, 0.084915856106222, 0.000000000000000],
        [1.000000000000000, 0.044125366682106, 0.000000000000000],
        [1.000000000000000, 0.003890690134321, 0.000000000000000],
        [0.962544387833874, 0.000000000000000, 0.029612474911224],
        [0.921784776902887, 0.000000000000000, 0.062590705573568],
        [0.878771035973445, 0.000000000000000, 0.097236374864906],
        [0.839678863671453, 0.000000000000000, 0.131882044156245],
        [0.796665122742010, 0.000000000000000, 0.165323452215532],
        [0.757572950440019, 0.000000000000000, 0.197251814111471],
        [0.714559209510576, 0.000000000000000, 0.231218156553960],
        [0.673089393237610, 0.000000000000000, 0.264165508723174],
        [0.632453296279141, 0.000000000000000, 0.297267253358036],
        [0.589439555349699, 0.000000000000000, 0.331912922649375],
        [0.550347383047707, 0.000000000000000, 0.363007565230817],
        [0.507333642118265, 0.000000000000000, 0.397282692604601],
        [0.468241469816273, 0.000000000000000, 0.431928361895940],
        [0.425227728886831, 0.000000000000000, 0.465771190365910],
        [0.383634398641346, 0.000000000000000, 0.498718542535124],
        [0.343121815655396, 0.000000000000000, 0.531943801142504],
        [0.380299521383357, 0.058669136946117, 0.562667901806392],
        [0.421337038752509, 0.124038906901344, 0.589470433842829],
        [0.464350779681951, 0.185487108229118, 0.620194534506716],
        [0.503442951983944, 0.246935309556894, 0.650146672842365],
        [0.546456692913386, 0.311378724718234, 0.677721167207040],
        [0.585548865215378, 0.373753280839895, 0.708445267870928],
        [0.628562606144820, 0.435201482167670, 0.735649220318048],
        [0.669754515979620, 0.498471514590088, 0.765971900571252],
        [0.710668519376254, 0.562019453450671, 0.796572487262621],
        [0.753682260305697, 0.623467654778447, 0.823498533271577],
        [0.792774432607688, 0.685564304461942, 0.854222633935464],
        [0.835788173537132, 0.750285626061448, 0.882075034738305],
        [0.874880345839123, 0.811733827389223, 0.911749266635788],
        [0.917894086768566, 0.873182028716999, 0.942473367299676],
        [0.959209510575884, 0.938551798672224, 0.969275899336112],
        [1.000000000000000, 1.000000000000000, 1.000000000000000]];
        
% Interpolate to requested number of colors.

P = size(cmap,1);

if (P ~= M)
  cmap = interp1(1:size(cmap,1), cmap, linspace(1,P,M), 'linear');
end

return