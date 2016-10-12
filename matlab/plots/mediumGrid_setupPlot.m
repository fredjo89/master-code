clc; clear all; close all; 



cores = [
1
8
16
32
64
128
320
640
1600
];

MPI_times =[
0.1252
0.1273
0.124
0.1282
0.1383
0.1321
0.1635
0.1554
0.4259
]*1000;

Hybrid_times = [
0
0.031650
0.1394
0.1292
0.1312
0.1448
0.1447
0.1291
0.1404
]*1000;


Y = [MPI_times' ; Hybrid_times']';



MPI_speedup = MPI_times(1)./MPI_times;
Hybrid_speedup = MPI_times(1)./Hybrid_times;



height1 = [ 1 1 1 1 1 1 1 1 -5 1 ]*500/35;
height2 = [ 1 1 1 1 1 1 1 1 1 1]*500/35;

my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;


HEY = [MPI_speedup'; Hybrid_speedup'];

HEY= round(HEY,2)

FigHandle = figure('Position', [1200, 200, 800*11/9, 500]);
h = bar(Y, 1);
drawnow
opts = {'VerticalAlign','middle', 'HorizontalAlign','left', ...
    'FontSize',15, 'Rotation',90};
    clr = h(1).Face.ColorData(1:3);
    vd = h(1).Face.VertexData;
    xy = double(vd(1:2,2:4:end) + vd(1:2,4:4:end)) / 2
    for j=1:size(xy,2)
        if j ==9
            text(xy(1,j), xy(2,j)+height1(j)  , num2str(HEY(1,j)), 'Color', 'w', opts{:});
        else
            text(xy(1,j), xy(2,j)+height1(j)  , num2str(HEY(1,j)), 'Color', 'k', opts{:});
        end
    end
    clr = h(2).Face.ColorData(1:3);
    vd = h(2).Face.VertexData;
    xy = double(vd(1:2,2:4:end) + vd(1:2,4:4:end)) / 2
  
    for j=2:length(Y)
        text(xy(1,j), xy(2,j)+height2(j)  , num2str(HEY(2,j)), 'Color', 'k', opts{:});
    end
 
LEG1 = legend('Pure MPI', 'Hybrid');
set(gca, 'XTickLabel',[])                    
xt = get(gca, 'XTick');
set(gca, 'XTick', [1:1:16], 'XTickLabel', cores)
set(gca,'ytick',[0:100:500]);
xlabel('Number of cores');
ylabel('Time (ms)');
colormap([my_red_1; my_blue_1; my_green_1])
set(gca,'fontsize',17)
set(gca,'XLim',[0.5 9.5],'YLim',[0 430])


   
%}



%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters mediumGrid_setup
%}


%print -dpng fluxFig2












