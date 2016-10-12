clc; clear all; close all; 

MPI =[
153.30
77.62
41.15
27.29
13.80
6.97
4.78
3.61
2.02
1.13
0.81
0.76
0.56
];

cores_MPI=[
1
2
4
8
16
32
48
64
128
320
640
960
1600
];

Hybrid_16 =[
17.96
17.2500
11.27
8.315
4.128
0.8576
0.5826
0.5267
0.4773
];

Hybrid_8 =[
27.16
13.48
6.7530
4.577
3.421
1.768
0.8509
0.6141
0.4686
0.3703
];

cores_Hybrid=[
8
16
32
48
64
128
320
640
960
1600
];


MPI_sp = MPI(1)./MPI;
Hybrid16_sp = MPI(1)./Hybrid_16;

Hybrid8_sp = MPI(1)./Hybrid_8;

my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;


TIC_ONE = [ 128 320 640 960 1600];

%FigHandle = figure('Position', [1200, 200, 500, 500]);
FigHandle = figure('Position', [1200, 700, 13*29, 11.5*29]);
hold on; 
plot(cores_MPI,MPI_sp,'--o','Color',my_green_1, 'LineWidth', 1, 'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1); 
plot(cores_Hybrid,Hybrid8_sp,'--o','Color',my_red_1, 'LineWidth', 1, 'MarkerSize', 7, 'MarkerEdgeColor', my_red_1, 'MarkerFaceColor', my_red_1); 
plot(cores_MPI,cores_MPI,'Color',my_blue_1 , 'LineWidth', 1.5) 

LEG1 = legend('Pure MPI','Hybrid','Perfect speedup');
xlabel('Number of cores');
ylabel('Speedup');
axis([1,1600,1,450]);
set(gca,'xtick', TIC_ONE);
set(gca,'ytick', [0:100:450]);

TIC_TWO = [ 1 8 16 32 48 64];

theEnd = 8; 
theEnd2 = 5; 

%FigHandle = figure('Position', [1200, 200, 500, 500]);
FigHandle = figure('Position', [1200, 200, 13*29, 11.5*29]);
hold on; 
plot(cores_MPI(1:theEnd),MPI_sp(1:theEnd),'--o','Color',my_green_1, 'LineWidth', 1,'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1); 
plot(cores_Hybrid(1:theEnd2),Hybrid8_sp(1:theEnd2),'--o','Color',my_red_1, 'LineWidth', 1, 'MarkerSize', 7, 'MarkerEdgeColor', my_red_1, 'MarkerFaceColor', my_red_1); 
plot(cores_MPI(1:theEnd),cores_MPI(1:theEnd),'Color',my_blue_1 , 'LineWidth', 1.5) 

LEG1 = legend('Pure MPI','Hybrid','Perfect speedup');
xlabel('Number of cores');
ylabel('Speedup');
axis([1,cores_MPI(theEnd),1,50]);
set(gca,'xtick', TIC_TWO);
set(gca,'ytick', [0:10:50]);




%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters SPE_speedup_2
%}









