clc; clear all; close all; 





cores = [1 8 16 32 64 128 320 640 1600 ];

MPI_times =[
0.5303
1.445
1.211
1.408
1.24
1.017
0.9148
0.88
1.457
]*1000;





Hybrid_times = [
0
0.1387
0.6814
1.03
1.458
1.254
1.343
1.27
0.9768
]*1000;


Y = [
    MPI_times(1) Hybrid_times(1);
    MPI_times(2) Hybrid_times(2);
    MPI_times(3) Hybrid_times(3);
    MPI_times(4) Hybrid_times(4);
    MPI_times(5) Hybrid_times(5);
    MPI_times(6) Hybrid_times(6);
    MPI_times(7) Hybrid_times(7);
    MPI_times(8) Hybrid_times(8);
    MPI_times(9) Hybrid_times(9);
];

MPI_speedup = MPI_times(1)./MPI_times;
Hybrid_speedup = MPI_times(1)./Hybrid_times;




my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;
my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;
my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;


height1 = [ 1 -7 1 -7 1 1 1 1 -7 1 1 1 ]*1100/35;
height2 = [ 1 1 1 1 -7 1 -7 1 1 1 1 1 ]*1100/35;

HEY = [MPI_speedup'; Hybrid_speedup'];

HEY= round(HEY,2)

FigHandle = figure('Position', [1200, 200, 800*13/9, 500]);
h = bar(Y, 1);
drawnow
opts = {'VerticalAlign','middle', 'HorizontalAlign','left', ...
    'FontSize',15, 'Rotation',90};
    clr = h(1).Face.ColorData(1:3);
    vd = h(1).Face.VertexData;
    xy = double(vd(1:2,2:4:end) + vd(1:2,4:4:end)) / 2
    for j=1:size(xy,2)
        if (j == 2 || j ==4 || j == 9)
             text(xy(1,j), xy(2,j)+height1(j)  , num2str(HEY(1,j)), 'Color', 'w', opts{:});
        else 
            text(xy(1,j), xy(2,j)+height1(j)  , num2str(HEY(1,j)), 'Color', 'k', opts{:});
        end
    end
    clr = h(2).Face.ColorData(1:3);
    vd = h(2).Face.VertexData;
    xy = double(vd(1:2,2:4:end) + vd(1:2,4:4:end)) / 2
  
    for j=2:length(Y)
        if (j ==5 || j==7)
            text(xy(1,j), xy(2,j)+height2(j)  , num2str(HEY(2,j)), 'Color', 'w', opts{:});
        else
        text(xy(1,j), xy(2,j)+height2(j)  , num2str(HEY(2,j)), 'Color', 'k', opts{:});
        end
    end

 
LEG1 = legend('Pure MPI', 'Hybrid');
set(gca, 'XTickLabel',[])                    
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', cores)
set(gca,'ytick',[0:500:1500]);
xlabel('Number of cores');
ylabel('Time (ms)');
colormap([my_red_1; my_blue_1; my_green_1])
set(gca,'fontsize',17)
set(gca,'XLim',[0.5 9.5],'YLim',[0 1500])
   
%}



%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters SPE_setup
%}


%print -dpng fluxFig2












