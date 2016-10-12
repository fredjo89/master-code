clc; clear all; close all; 

MPI =[
819.70
413.90
223.30
156.80
78.50
39.29
19.77
3.97
1.93
0.82
0.65
];

cores_MPI=[
1
2
4
8
16
32
64
320
640
1280
1600
];

Hybrid_8 =[
78.92
19.88
4.026
1.338
0.7673
0.5693
];

cores_Hybrid=[
16
64
320
640
1280
1600
];

MPI_sp = MPI(1)./MPI;
Hybrid_sp = MPI(1)./Hybrid_8;

my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;




%FigHandle = figure('Position', [1200, 200, 500, 500]);
FigHandle = figure('Position', [1200, 200, 13*29, 11.5*29]);
hold on; 
plot(cores_MPI,MPI_sp,'--o','Color',my_green_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1); 
plot(cores_Hybrid,Hybrid_sp,'--o','Color',my_red_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_red_1, 'MarkerFaceColor', my_red_1); 
plot(cores_MPI,cores_MPI,'Color',my_blue_1 , 'LineWidth', 1.5) 

LEG1 = legend('Pure MPI','Hybrid','Perfect speedup');
xlabel('Number of cores');
ylabel('Speedup');
axis([1,32,1,32]);
set(gca,'xtick',[1  2 4 8 16 32]);
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
print -dpdf -painters largeGrid_speedup_closeup
%}






















