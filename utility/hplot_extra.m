% This script adds extra characteristics to the plot.

% EXTRA='WC13';
% EXTRA='WC13_Vname';
% EXTRA='ECCOFS';
% EXTRA='ECCOFS_Vname';
  EXTRA='None';

switch EXTRA
  case 'WC13'
    [x,y]=m_ll2xy(-128.292,37.918);     % WC13 T,S observation
    plot(x, y, 'o','MarkerSize',8,                                   ...
        'MarkerEdgeColor', 'r', 'MarkerFaceColor',[0.8,0.8,0.80]);
    x=[-134, -122.5];
    y=[37.666 37.666];
    [x,y]=m_ll2xy(x,y);                 % WC13 cross-section
    plot(x,y,'r:');
  case 'WC13_Vname'
    x=-121;
    y=47.5;
    [x,y]=m_ll2xy(x,y);                 % WC13 variable name
    text(x,y, untexlabel(P.Vname),'FontSize',20,                     ...
         'FontWeight','bold','Color','k'); 
  case 'ECCOFS'
    x1=[-76 -70]; y1=[35 35];           % Cape Hatteras along 35N
    x2=[-86.9 -82.1]; y2=[21.5 26.6];   % Loop Current
    x3=[-68 -68]; y3=[20 45];           % Basin along 68W

    [x,y]=m_ll2xy(x1,y1);   plot(x,y,'w-', 'LineWidth', 3);
    [x,y]=m_ll2xy(x2,y2);   plot(x,y,'w-', 'LineWidth', 3);
    [x,y]=m_ll2xy(x3,y3);   plot(x,y,'w-', 'LineWidth', 3);
  case 'ECCOFS_Vname'
    x=-96;
    y=38;
    [x,y]=m_ll2xy(x,y);                 % WC13 variable name
    text(x,y, untexlabel(P.Vname),'FontSize',20,                     ...
         'FontWeight','bold','Color','k'); 
  otherwise
    % skip
end