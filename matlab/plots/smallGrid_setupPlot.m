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

cores = [1 2 4 8 16 32 64 80];

MPI_times =[
7.29
9.74
8.76
8.94
13.79
9.91
13.11
18.33
];



OMP_times = [
7.29
4.49
3.16
2.68
4.16
];

MPI_speedup = MPI_times(1)./MPI_times;
OMP_speedup = OMP_times(1)./OMP_times;

OMP_speedup(6) = 0; 
OMP_speedup(7) = 0;
OMP_speedup(8) = 0;

Y = [
    MPI_times(1) OMP_times(1);
    MPI_times(2) OMP_times(2);
    MPI_times(3) OMP_times(3);
    MPI_times(4) OMP_times(4);
    MPI_times(5) OMP_times(5);
    MPI_times(6) 0;
    MPI_times(7) 0;
    MPI_times(8) 0;
];


height1 = [ 1 1 1 1 1 1 1 1 1 -4 ]*20/35;
height2 = [ 1 1 1 1 1 1 1 1 1 1]*20/35;

my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;
my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;
my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;


HEY = [MPI_speedup'; OMP_speedup'];

HEY= round(HEY,2)

FigHandle = figure('Position', [1200, 200, 800, 500]);
h = bar(Y, 1);
drawnow
opts = {'VerticalAlign','middle', 'HorizontalAlign','left', ...
    'FontSize',15, 'Rotation',90};
    clr = h(1).Face.ColorData(1:3);
    vd = h(1).Face.VertexData;
    xy = double(vd(1:2,2:4:end) + vd(1:2,4:4:end)) / 2
    for j=1:size(xy,2)
        if (j==8)
            text(xy(1,j), xy(2,j) -2  , num2str(HEY(1,j)), 'Color', 'w', opts{:});
        else
            text(xy(1,j), xy(2,j) +height1(j)  , num2str(HEY(1,j)), 'Color', 'k', opts{:});
        end
    end
    clr = h(2).Face.ColorData(1:3);
    vd = h(2).Face.VertexData;
    xy = double(vd(1:2,2:4:end) + vd(1:2,4:4:end)) / 2
  
    for j=1:length(OMP_times)
        text(xy(1,j), xy(2,j)+ height1(j) , num2str(HEY(2,j)), 'Color', 'k', opts{:});
    end
 
LEG1 = legend('MPI', 'OpenMP');
set(gca, 'XTickLabel',[])                    
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', cores)
set(gca,'ytick',[0:5:20]);
xlabel('Number of cores');
ylabel('Time (ms)');
colormap([my_red_1; my_blue_1; my_green_1])
set(gca,'fontsize',17)
set(gca,'XLim',[0.5 8.25],'YLim',[0 20])
   
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

LEG1 = legend('MPI', 'Setup');
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
print -dpdf -painters smallGrid_setup
%}


%print -dpng fluxFig2












