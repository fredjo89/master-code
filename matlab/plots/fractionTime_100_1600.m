clc; clear all; close all; 

% 20 x 20 blocks
% upscaling factor: 100


T_total = [
819.7
78.5000
39.29
19.77
3.97
1.931
0.8183
0.6506
];

T_jacobi = [
782.6
74.8600
37.62
18.42
3.332
1.753
0.5453
0.4244
];

T_send = [
0
0.1683
0.2641
0.6817
0.5303
0.1382
0.2528
0.2104
];

T_norm = T_total - T_jacobi - T_send;



T_j_frac = T_jacobi./T_total

T_s_frac = T_send./T_total

T_n_frac = T_norm./T_total;


x = [
1
16
32
64
320
640
1280
1600
];


my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;



%FigHandle = figure('Position', [1200, 200, 500, 500]);
FigHandle = figure('Position', [1200, 200, 13*29, 11.5*29]);


xt = [0 1 2 3 4 5 6 7  ];
hold on; 
plot(xt, T_j_frac,'--o','Color',my_blue_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_blue_1, 'MarkerFaceColor', my_blue_1); 
plot(xt,T_s_frac,'--o','Color',my_green_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1); 
plot(xt,T_n_frac,'--o','Color',my_red_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_red_1, 'MarkerFaceColor', my_red_1); 

set(gca, 'XTickLabel',[])                    

yl = get(gca, 'YLim');
str = cellstr( num2str(x(:)) );      %# format x-ticks as 2^{xx}
hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
    'Interpreter','tex', ...                   %# specify tex interpreter
    'VerticalAlignment','top', ...             %# v-align to be underneath
    'HorizontalAlignment','center' ); 


LEG1 = legend('Smoothing','Sending', 'Normalizing');
xlabel('Number of cores');
ylabel('Fraction');
axis([0,11,0,1]);
set(gca,'xtick',x);
set(gca,'ytick',[0:0.2:1]);
set(LEG1);



xt = [1:length(x)];

xt = [1 2 3 4 5 6 7 8.3];

newY = [ T_s_frac T_n_frac T_j_frac   ]

FigHandle = figure('Position', [1200, 200, 80*length(x), 400]);
bar(newY,'stacked');  %# Create a stacked histogram
set(gca, 'XTickLabel',[])                    
set(gca, 'XTick', x)
set(gca,'ytick',[0:0.2:1]);
xlabel('Number of cores');
ylabel('Fraction');
LEG1 = legend('Sending','Normalizing', 'Smoothing');
LEG1 = legend('Message-passing','Normalizing', 'Smoothing');
colormap([my_red_1; my_blue_1; my_green_1])
set(gca,'fontsize',17)


set(gca, 'XTickLabel',[])                    
%xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', x)
set(gca,'ytick',[0:0.2:2]);



set(gca,'XLim',[0 9],'YLim',[0 1])









%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters epsFig_large
%}


