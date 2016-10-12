clc; close all; clear all; 


my_green_1 = [93 148 111] ./ 255;

lineW = 5; 

x = [0 1.015];
y = [0 0]; 
plot(x,y, 'color', my_green_1, 'linewidth',lineW); 
hold on; 

x = [ 1 1]; 
y = [ 0 1.01]; 
plot(x,y, 'color', my_green_1, 'linewidth',lineW); 

x = [ 1 2.015]; 
y = [ 1 1]; 
plot(x,y, 'color', my_green_1, 'linewidth',lineW); 

x = [ 2 2]; 
y = [ 1 0]; 
plot(x,y, 'color', my_green_1, 'linewidth',lineW); 

x = [ 2-0.015 3]; 
y = [ 0 0]; 
plot(x,y, 'color', my_green_1, 'linewidth',lineW); 

%set(gca, 'XTick', xt, 'XTickLabel', x)
%set(gca,'ytick',[0:0.2:2]);


x = [ 1 2 ];

x = ['a' 'b' ]
xt = [1 2  ];

yl = get(gca, 'YLim');
str = cellstr( num2str(x(:)) );      %# format x-ticks as 2^{xx}
hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
    'Interpreter','tex', ...                   %# specify tex interpreter
    'VerticalAlignment','top', ...             %# v-align to be underneath
    'HorizontalAlignment','center', 'fontsize',20 ); 

set(gca,'fontsize',20)

%LEG1 = legend('Smoothing','Sending', 'Normalizing');
%set(LEG1);
%xlabel('Number of cores');
%ylabel('Fraction');
%axis([0,11,0,1]);
%set(gca,'xtick',x);
%set(gca,'ytick',[0:0.2:1]);



set(gca, 'XTickLabel',[])   
set(gca,'ytick',[0 1]);
set(gca,'XLim',[0.3 2.7],'YLim',[0 1.1])
set(gca,'linewidth',2)


text (0.35,0.93,'f_{a,b}(x)', 'fontsize',25) 

text (2.65,-0.08,'x', 'fontsize',25) 
 
 


% the arrows


patch(...
    [2.7 0.24; 2.7 0.36; 2.8 0.3], ...
    [0.05 1.1; -0.05 1.1; 0 1.18], 'k', 'clipping', 'off')














%{
% Some bogus functions
f = @(x) 50* 1.6.^(-x-5);
g = @(x) 50* 1.6.^(+x-10);

% Point where they meet
xE = 2.5;
yE = f(xE);

% Plot the bogus functions
figure(1), clf, hold on
x = 0:0.2:5;
plot(x,f(x),'r',  x,g(x),'b', 'linewidth', 2)

% get rid of standard axes decorations
set(gca, 'Xtick', [], 'Ytick', [], 'box', 'off')

% Fix the axes sizes
axis([0 5 0 5])

% the equilibrium point
plot(xE, yE, 'k.', 'markersize', 20)

% the dashed lines
line([xE 0; xE xE], [0 yE; yE yE], 'linestyle', '--', 'color', 'k')

% the arrows
xO = 0.2;  
yO = 0.1;
patch(...
    [5-xO -yO; 5-xO +yO; 5.0 0.0], ...
    [yO 5-xO; -yO 5-xO; 0 5], 'k', 'clipping', 'off')

% the squishy wiggly line pointing to the "equilibrium" text
h = @(x)0.5*(x+0.2) + 0.1*sin((x+0.2)*14);
x = 2.7:0.01:3.5;
plot(x, h(x), 'k', 'linewidth', 2)

% the static texts
text(xE-yO, -0.2, 'Q^*', 'fontweight', 'bold')
text(-2*yO,   yE, 'P^*', 'fontweight', 'bold')
text(-2*yO,    4, 'Price', 'rotation', 90, 'fontsize', 14)
text(    4, -0.2, 'Quantity', 'fontsize', 14)
text(   .5,  4.2, 'Demand', 'fontsize', 14, 'rotation', -55)
text(   4.0,  3.3, 'Supply', 'fontsize', 14, 'rotation', +55)
text(   3.6,  2.1, 'Equilibrium', 'fontsize', 14)
%}






















