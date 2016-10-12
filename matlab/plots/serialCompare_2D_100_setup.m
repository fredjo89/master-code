clc; clear all; close all; 


OMP_iter_100=[	
0.0132561	0.0536571	0.121182	0.2157	0.334763
]*1000;


OMP_iter_100_round= round(OMP_iter_100,1);

MPI_iter_100 =[
0.007287	0.03081	0.07002	0.127	0.2043
]*1000;


Speedup = OMP_iter_100./MPI_iter_100;
Speedup= round(Speedup,2)

Y = [
    OMP_iter_100(1) MPI_iter_100(1);
    OMP_iter_100(2) MPI_iter_100(2);
    OMP_iter_100(3) MPI_iter_100(3);
    OMP_iter_100(4) MPI_iter_100(4);
    OMP_iter_100(5) MPI_iter_100(5);
];


%% Plotting

height1 = [ 1 1 1 1 -8 ]*350/35;
height2 = [ 1 1 1 1 1 ]*350/35;
xString = [' 20 x 20 '; ' 40 x 40 ';' 60 x 60 ';' 80 x 80 ';'100 x 100']



my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;
my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;
my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;
my_red_3 = [241 36 35] ./ 255;

FigHandle = figure('Position', [1200, 200, 650, 500]);
h = bar(Y, 1, 'BarWidth', 4);

drawnow
opts = {'VerticalAlign','middle', 'HorizontalAlign','left', ...
    'FontSize',18, 'Rotation',90};
    clr = h(1).Face.ColorData(1:3);
    vd = h(1).Face.VertexData;
    xy = double(vd(1:2,2:4:end) + vd(1:2,4:4:end)) / 2;
        for j=1:length(Y)
        text(xy(1,j), xy(2,j)+height1(j)  , num2str(OMP_iter_100_round(j)), 'Color', 'k', opts{:});
    end

    clr = h(2).Face.ColorData(1:3);
    vd = h(2).Face.VertexData;
    xy = double(vd(1:2,2:4:end) + vd(1:2,4:4:end)) / 2;
    for j=1:length(Y)
        text(xy(1,j), xy(2,j)+height2(j)  , num2str(Speedup(j)), 'Color', 'b', opts{:});
    end

LEG1 = legend('Shared ', 'Distributed');
set(gca, 'XTickLabel',[])                    
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xString)
set(gca,'ytick',[0:50:350]);
xlabel('Number of coarse blocks');
ylabel('Time (ms)');
colormap([my_red_1; my_blue_1; my_green_1])
set(gca,'fontsize',15)

set(gca,'XLim',[0.5 5.5],'YLim',[0 350])
xlhand = get(gca,'xlabel')
ylhand = get(gca,'ylabel')

%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters serialCompare_setup_100
%}













