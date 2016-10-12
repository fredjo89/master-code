clc; clear all; close all;

% load SPE 10 data set
mrstModule add spe10;
rock = SPE10_rock(); p=rock.poro; K=rock.perm;
%{
% show p
slice( reshape(p,60,220,85), [1 220], 60, [1 85]);
shading flat, axis equal off, set(gca, 'zdir', 'reverse'), box on;
colorbar('horiz');
%}







% show Kx
FigHandle = figure('Position', [200, 1000, 1000, 600]);
slice( reshape(log10(K(:,1)),60,220,85), [1 220], 60, [1 85]);
shading flat, axis equal off, set(gca, 'zdir', 'reverse'), box on;
h=colorbar();
set(h, 'XTickLabel',10.^[get(h, 'XTick')]);
set(h, 'YTick',mean(get(h, 'YLim')), 'YTickLabel', 'mD' );
colormap(jet(128)); 
ax = gca; 

axpos = ax.Position;
cpos = h.Position;
cpos(1) = 0.85;
cpos(2) = 0.15;
cpos(3) = 1.5*cpos(3);
cpos(4) = 0.9*cpos(4);
h.Position = cpos;
ax.Position = axpos;

ticks = [  0.001 0.01 0.1 1 10 100 1000 10000]
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks);
colormap(jet(128)); 
set(gca, 'FontSize',20);
set(get(h,'title'),'string','mD');
view(-135-5,30)

%print -dpng SPE_3D




FigHandle = figure('Position', [1000, 1000, 1000, 600]);
mrstModule add spe10
layers = 25:25;
[G, ~, rock] = SPE10_setup(layers);
rock.perm = rock.perm.*1.0132*10^(15);
R = plotCellData(G,log10(rock.perm(:,1)));
rotate( R, [0 0 1], 90)
shading flat, axis equal off, box on;
h=colorbar( 'horiz');
set(h, 'XTickLabel',10.^[get(h, 'XTick')]);
set(h, 'YTick',mean(get(h, 'YLim')), 'YTickLabel', 'mD' );
colormap(jet(128)); 

ax = gca; 
axpos = ax.Position;
cpos = h.Position;
%cpos(1) = 0.2;
cpos(2) = cpos(2)*1.1;
%cpos(3) = 0.95*cpos(3);
cpos(4) = 1.5*cpos(4);
h.Position = cpos;
ax.Position = axpos;
ticks = [  0.001 0.01 0.1 1 10 100 1000 10000]
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks);
colormap(jet(128)); 
set(gca, 'FontSize',20);

set(get(h,'title'),'string','mD');



max(rock.perm(:,1))
min(rock.perm(:,1))

max(rock.perm(:,1))/min(rock.perm(:,1))


%print -dpng SPE_layer_25

%{
max(rock.perm)
min(rock.perm)

max(K);
min(K);


rock.perm = convertFrom(rock.perm,milli*darcy);

max(rock.perm);
min(rock.perm);

x = rock.perm(1,1)

y = convertFrom(x,milli*darcy)

x/y

%}
























