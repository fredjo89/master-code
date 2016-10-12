clc; clear all; close all; 



%{
1	&7.29		&1.00		&7.29		&1.00\\
2	&9.74		&1.34		&4.49		&0.62\\
4	&8.76		&1.20		&3.16		&0.43\\
8	&8.94		&1.23		&2.68		&0.37\\
16	&13.79		&1.89		&4.16		&0.57\\
32	&9.91		&1.36		&	-		&-\\
48	&11.11		&1.52		&	-		&-\\
64	&13.11		&1.80		&	-		&-\\
80	&18.33		&2.51		&	-		&-\\
%}

cores = [1  8 16 32 64 320 640 1600 ];

MPI_times =[
3576
3220
3077
3131
3240
3391
3570
4720
];



Hybrid_times = [
0
0.851800
3.436
3.259
3.34
3.251
3.33
3.349
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
];

MPI_speedup = MPI_times(1)./MPI_times;
Hybrid_speedup = MPI_times(1)./Hybrid_times;




height1 = [ 1 1 1 1 1 1 1 -5 -5 -4 ]*5000/35;
height2 = [ 1 1 1 1 1 1 1 1 1 1]*5000/35;

my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;
my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;
my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;

my_colour = [236 250 239] ./ 255;


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
        if (j==8 || j == 9)
        text(xy(1,j), xy(2,j)+height1(j)  , num2str(HEY(1,j)), 'Color', 'w', opts{:});
        else
        text(xy(1,j), xy(2,j)+height1(j)  , num2str(HEY(1,j)), 'Color', 'k', opts{:});
        end
        j
    end
    clr = h(2).Face.ColorData(1:3);
    vd = h(2).Face.VertexData;
    xy = double(vd(1:2,2:4:end) + vd(1:2,4:4:end)) / 2
  
    for j=2:length(Y)
        text(xy(1,j)+0.02, xy(2,j)+height2(j)  , num2str(HEY(2,j)), 'Color', 'k', opts{:});
    end
 
LEG1 = legend('Pure MPI', 'Hybrid');
set(gca, 'XTickLabel',[])                    
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', cores)
set(gca,'ytick',[0:1000:5000]);
xlabel('Number of cores');
ylabel('Time (ms)');
colormap([my_red_1; my_blue_1; my_green_1])
set(gca,'fontsize',17)
set(gca,'XLim',[0.5 8.5],'YLim',[0 5000])
   
%}

%{
figure();
Y = R;
h = bar(Y);
drawnow   % this is needed for some reason!

opts = {'VerticalAlign','middle', 'HorizontalAlign','left', ...
    'FontSize',12, 'Rotation',90};
for i=1:numel(h)
    clr = h(i).Face.ColorData(1:3);
    vd = h(i).Face.VertexData;
    xy = double(vd(1:2,2:4:end) + vd(1:2,4:4:end)) / 2;
    for j=1:size(xy,2)
        text(xy(1,j), xy(2,j), sprintf(' %.2g',xy(2,j)), ...
            'Color','k', opts{:})
    end
end

LEG1 = legend('Iteration', 'Setup');
set(gca, 'XTickLabel',[])                    
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', [100 250 1000])
set(gca,'ytick',[0:0.2:1]);
xlabel('Upscaling factor');
ylabel('Ratio');
colormap([my_red_1; my_blue_1; my_green_1])
set(gca,'fontsize',13)

%}

%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters largeGrid_setup
%}


%print -dpng fluxFig2












