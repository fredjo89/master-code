clc; clear all; close all; 

% 20 x 20 blocks
% upscaling factor: 100


T_total = [
1.742   
0.8827	
0.4422	
0.2305	
0.1259 
0.07193 
0.06431
];

T_jacobi = [
1.559000
0.771400
0.383900
0.192500
0.102000
0.050900
0.035370
];

T_send = [
0
0.01858
0.01167
0.01469
0.01172
0.01428
0.02396
];

T_norm = T_total - T_jacobi - T_send;



T_j_frac = T_jacobi./T_total

T_s_frac = T_send./T_total

T_n_frac = T_norm./T_total;


x = [1 2 4 8 16 32 48];


my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;



%FigHandle = figure('Position', [1200, 200, 500, 500]);
FigHandle = figure('Position', [1200, 200, 13*29, 11.5*29]);

hold on; 
plot(log2(x),T_j_frac,'--o','Color',my_blue_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_blue_1, 'MarkerFaceColor', my_blue_1); 
plot(log2(x),T_s_frac,'--o','Color',my_green_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1); 
plot(log2(x),T_n_frac,'--o','Color',my_red_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_red_1, 'MarkerFaceColor', my_red_1); 

set(gca, 'XTickLabel',[])                    
xt = get(gca, 'XTick');
yl = get(gca, 'YLim');
str = cellstr( num2str(x(:)) );      %# format x-ticks as 2^{xx}
hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
    'Interpreter','tex', ...                   %# specify tex interpreter
    'VerticalAlignment','top', ...             %# v-align to be underneath
    'HorizontalAlignment','center' ); 




LEG1 = legend('Smoothing','Sending', 'Normalizing');
xlabel('Number of cores');
ylabel('Fraction');
axis([0,6,0,1]);
set(gca,'xtick',x);
set(gca,'ytick',[0:0.2:1]);
set(LEG1);



set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters epsFig



