clc; clear all; close all; 

% 20 x 20 blocks
% upscaling factor: 100


MPI_iter = [
1.742   0.8827	0.4422	0.2305	0.1259 0.07193 0.06431
];



x = [1 2 4 8 16 32 48];

MPI_iter = MPI_iter(1)./MPI_iter

x2 = log(x); 







my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;



%FigHandle = figure('Position', [1200, 200, 500, 500]);
FigHandle = figure('Position', [1200, 200, 13*29, 11.5*29]);

hold on; 
plot(log2(x),MPI_iter,'--o','Color',my_green_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1); 
plot(log2(x),x,'Color',my_blue_1 , 'LineWidth', 1.5)              %# plot on log2 x-scale
set(gca, 'XTickLabel',[])                      %# suppress current x-labels
xt = get(gca, 'XTick');
yl = get(gca, 'YLim');
str = cellstr( num2str(x(:)) );      %# format x-ticks as 2^{xx}
hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
    'Interpreter','tex', ...                   %# specify tex interpreter
    'VerticalAlignment','top', ...             %# v-align to be underneath
    'HorizontalAlignment','center' ); 




LEG1 = legend('MPI speed-up','Perfect speed-up');
xlabel('Number of processors');
ylabel('Speed-up');
axis([0,6,0,32]);
set(gca,'xtick',x);
set(gca,'ytick',[0:8:48]);
set(LEG1);






set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters epsFig


