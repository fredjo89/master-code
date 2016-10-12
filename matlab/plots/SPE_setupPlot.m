clc; clear all; close all; 





cores = [1 2 4 8 16 32 48 64 128 320 640 960 1600 ];

MPI_times =[
0.5303
0.759
1.056
1.445
1.211
1.408
1.351
1.24
1.017
0.9148
0.88
0.84
1.457
]*1000;



Hybrid_times = [
0
0
0
0
0.1283
0.7025
0.8919
1.049
1.491
1.595
1.356
1.253
1.262
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
    MPI_times(10) Hybrid_times(10);
    MPI_times(11) Hybrid_times(11);
    MPI_times(12) Hybrid_times(12);
    MPI_times(13) Hybrid_times(13);
];

MPI_speedup = MPI_times(1)./MPI_times;
Hybrid_speedup = MPI_times(1)./Hybrid_times;




my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;


HEY = [MPI_speedup'; Hybrid_speedup'];

HEY= round(HEY,2)

FigHandle = figure('Position', [1200, 200, 800*13/9, 500]);
h = bar(Y, 1);
drawnow
opts = {'VerticalAlign','middle', 'HorizontalAlign','left', ...
    'FontSize',14, 'Rotation',90};
    clr = h(1).Face.ColorData(1:3);
    vd = h(1).Face.VertexData;
    xy = double(vd(1:2,2:4:end) + vd(1:2,4:4:end)) / 2
    for j=1:size(xy,2)
        text(xy(1,j), 200  , num2str(HEY(1,j)), 'Color', 'k', opts{:});
    end
    clr = h(2).Face.ColorData(1:3);
    vd = h(2).Face.VertexData;
    xy = double(vd(1:2,2:4:end) + vd(1:2,4:4:end)) / 2
  
    for j=5:length(Y)
        text(xy(1,j), 200  , num2str(HEY(2,j)), 'Color', 'k', opts{:});
    end
 
LEG1 = legend('Pure MPI', 'Hybrid');
set(gca, 'XTickLabel',[])                    
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', cores)
set(gca,'ytick',[0:250:1500]);
xlabel('Number of cores');
ylabel('Time (ms)');
colormap([my_red_1; my_blue_1; my_green_1])
set(gca,'fontsize',15)

   
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
print -dpdf -painters SPE_setup
%}


%print -dpng fluxFig2












