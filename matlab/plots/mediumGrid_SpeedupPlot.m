clc; clear all; close all; 

MPI =[
32.73
16.7
8.807
6.321
3.0490
1.351
0.8432
0.517
0.4128
0.3198
0.2287
0.1226
0.09537
0.05636
0.06533
0.0608
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
80
128
160
320
640
960
1280
1600
];

Hybrid_8 =[
6.309000
3.085
1.35
0.7625
0.526
0.488
0.2712
0.256
0.132
0.07528
0.0591
0.05028
0.04779
];

cores_Hybrid=[
8
16
32
48
64
80
128
160
320
640
960
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


RAWR = 11;

%FigHandle = figure('Position', [1200, 200, 500, 500]);
FigHandle = figure('Position', [1200, 200, 13*29, 11.5*29]);
hold on; 
plot(cores_MPI(1:RAWR),MPI_sp(1:RAWR),'--o','Color',my_green_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1); 
plot(cores_Hybrid(1:RAWR),Hybrid_sp(1:RAWR),'--o','Color',my_red_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_red_1, 'MarkerFaceColor', my_red_1); 
plot(cores_MPI,cores_MPI,'Color',my_blue_1 , 'LineWidth', 1.5) 

LEG1 = legend('Pure MPI','Hybrid','Perfect speedup');
xlabel('Number of cores');
ylabel('Speedup');
axis([1,cores_MPI(RAWR),1,cores_MPI(RAWR)]);
set(gca,'xtick', [1  40 80 120 160]);
set(gca,'xtick', [1 40 80 120 160 ]);




RAWR=length(cores_MPI);

%FigHandle = figure('Position', [1700, 200, 500, 500]);
FigHandle = figure('Position', [1700, 200, 13*29, 11.5*29]);
hold on; 
plot(cores_MPI,MPI_sp,'--o','Color',my_green_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1); 
plot(cores_Hybrid,Hybrid_sp,'--o','Color',my_red_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_red_1, 'MarkerFaceColor', my_red_1); 
plot(cores_MPI,cores_MPI,'Color',my_blue_1 , 'LineWidth', 1.5) 



%set(gca,'ytick',[0:250:1500]);

LEG1 = legend('Pure MPI','Hybrid','Perfect speedup');
xlabel('Number of cores');
ylabel('Speedup');
axis([1,1600,1,1000]);
set(gca,'xtick', [1  250 500 750 1000 1250  1600]);
set(gca,'ytick', [1 200 400 600 800 1000]);








%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters mediumGrid_speedup_2
%}






















