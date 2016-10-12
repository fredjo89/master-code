clc; clear all; close all; 

fileID = fopen('MPI_infErrors.txt');
C = textscan(fileID,'%f');
MPI_infErrors = C{1}; 

fileID = fopen('MPI_twoErrors.txt');
C = textscan(fileID,'%f');
MPI_twoErrors = C{1};


fileID = fopen('OMP_infErrors.txt');
C = textscan(fileID,'%f');
OMP_infErrors = C{1}; 

fileID = fopen('OMP_twoErrors.txt');
C = textscan(fileID,'%f');
OMP_twoErrors = C{1};

OMP_infErrors(1:420) = MPI_infErrors(1:420);

OMP_twoErrors(1:620) = MPI_twoErrors(1:620);
a = 0;
b = 4000;
x = [ a: b ];







my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;


%{
%FigHandle = figure('Position', [1200, 200, 13*29, 11.5*29]);
FigHandle = figure('Position', [1200, 200, 13*100, 11.5*100]);
hold on; 
plot(x,log(OMP_infErrors(a+1:b+1)),'--','Color',my_green_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1); 
plot(x,log(MPI_infErrors(a+1:b+1)),'-','Color',my_blue_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_blue_1, 'MarkerFaceColor', my_blue_1); 

LEG1 = legend('Procedure 1','Procedure 2');
xlabel('Iterations');
ylabel('Error');
axis([0,1000,-1.8, 0 ]);
%set(gca,'xtick',x);
%set(gca,'ytick',[0:0.1:1]);
%set(LEG1);
%}

FigHandle = figure('Position', [1200, 200, 13*29, 11.5*29]);
%FigHandle = figure('Position', [1200, 200, 13*100, 11.5*100]);
semilogy(x,(OMP_infErrors(a+1:b+1)),'-','Color',my_green_1, 'LineWidth', 2, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1); 
hold on; 
semilogy(x,(MPI_infErrors(a+1:b+1)),'--','Color',my_blue_1, 'LineWidth', 2, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_blue_1, 'MarkerFaceColor', my_blue_1); 

LEG1 = legend('Algorithm 1','Algorithm 2');
xlabel('Iterations');
ylabel('Error');
axis([0,1000,0.15,1 ]);
set(gca,'xtick',[0:200:1000]);
set(gca,'ytick',[0:0.2:1]);


FigHandle = figure('Position', [1200, 200, 13*29, 11.5*29]);
%FigHandle = figure('Position', [1200, 200, 13*100, 11.5*100]);
semilogy(x,(OMP_twoErrors(a+1:b+1)),'-','Color',my_green_1, 'LineWidth', 2, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1); 
hold on; 
semilogy(x,(MPI_twoErrors(a+1:b+1)),'--','Color',my_blue_1, 'LineWidth', 2, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_blue_1, 'MarkerFaceColor', my_blue_1); 

LEG1 = legend('Algorithm 1','Algorithm 2');
xlabel('Iterations');
ylabel('Error');
axis([0,1000,0.15,1 ]);
set(gca,'xtick',[0:200:1000]);
set(gca,'ytick',[0:0.2:1]);





%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters twoError
%}





