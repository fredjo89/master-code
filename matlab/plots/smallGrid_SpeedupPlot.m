clc; clear all; close all; 

% 20 x 20 blocks
% upscaling factor: 100


MPI_iter = [
1.742   0.8827	0.4422	0.2305	0.1259 0.07193 0.06431
];



MPI_iter = [
1.742   
0.8827	
0.4422	
0.2305	
0.1259 
0.07193 
0.06431
0.04225
0.0437
];





MPI_iter = MPI_iter(1)./MPI_iter

my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;






x = [ 1 2 4 8 16 32 48 64 80];



xt = [0 1 2 3 4 5 6 7 8  ];
%FigHandle = figure('Position', [1200, 200, 500, 500]);
FigHandle = figure('Position', [1200, 200, 13*29, 11.5*29]);
hold on; 
plot(x,MPI_iter,'--o','Color',my_green_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1); 
plot(x,x,'Color',my_blue_1 , 'LineWidth', 1.5)              %# plot on log2 x-scale

LEG1 = legend('MPI speedup','Perfect speed-up');
xlabel('Number of cores');
ylabel('Speedup');
axis([1,80,0,60]);


%{
set(gca, 'XTickLabel',[]) 
yl = get(gca, 'YLim');
str = cellstr( num2str(temp(:)) );      %# format x-ticks as 2^{xx}
hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
    'Interpreter','tex', ...                   %# specify tex interpreter
    'VerticalAlignment','top', ...             %# v-align to be underneath
    'HorizontalAlignment','center' ); 
LEG1 = legend('MPI speedup','Perfect speed-up');
xlabel('Number of cores');
ylabel('Speedup');
axis([0,8,0,50]);
set(gca,'xtick',x);
set(gca,'ytick',[0:8:80]);
set(LEG1);
%}




%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters smallGrid_speedup
%}


