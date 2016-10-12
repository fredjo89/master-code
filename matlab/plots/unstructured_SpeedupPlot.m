clc; clear all; close all; 

MPI =[
51.2
25.25
13.23
7.856
4.428
2.956
3.471
2.257
4.068
3.857
3.082
3.811
];

T_send = [
0
0.2308
0.4563
0.8975
0.9004
1.205
2.456
1.514
3.566
3.363
2.797
3.531
];

noSend = MPI - T_send;


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
96
128
160
];


Hybrid_8 =[
7.571
3.76
2.12
1.724
1.381
1.247
1.17
0.9145
0.8853
];

cores_Hybrid=[
8
16
32
48
64
80
96
128
160
];


MPI_sp = MPI(1)./MPI;
Hybrid8_sp = MPI(1)./Hybrid_8;

noSend_sp = noSend(1)./noSend

my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;
my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;
my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;


TIC_ONE = [ 1 40 80 120 160];

%FigHandle = figure('Position', [1200, 200, 500, 500]);
FigHandle = figure('Position', [1200, 700, 13*29, 11.5*29]);
hold on; 
plot(cores_MPI,MPI_sp,'--o','Color',my_green_1, 'LineWidth', 1, 'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1); 
plot(cores_Hybrid,Hybrid8_sp,'--o','Color',my_red_1, 'LineWidth', 1, 'MarkerSize', 7, 'MarkerEdgeColor', my_red_1, 'MarkerFaceColor', my_red_1); 
plot(cores_MPI,cores_MPI,'Color',my_blue_1 , 'LineWidth', 1.5) 

LEG1 = legend('Pure MPI','Hybrid','Perfect speedup');
xlabel('Number of cores');
ylabel('Speedup');
axis([1,160,1,60]);
set(gca,'xtick', TIC_ONE);
set(gca,'ytick', [1 20 40 60]);








%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters unstructured_speedup
%}









