clc; clear all; close all; 

OMP_setup = [1.97754 3.8601 12.2648 ];
MPI_setup = [ 0.7139 0.6055 0.5303];

OMP_iter = [152.229  158.354 171.219 ];
MPI_iter = [146.7 150.3  153.3];

setup_speed = OMP_setup./MPI_setup;
iter_speed = OMP_iter./MPI_iter

iter_save_ratio = (OMP_iter-MPI_iter)./OMP_iter

setup_save_ratio = (OMP_setup-MPI_setup)./OMP_setup

y = [ 152.229 146.7; 158.354 150.3; 171.219 153.3  ];


%bar(y)
%bar(iter_save_ratio)
%bar(setup_save_ratio)

R = [ iter_save_ratio(1) setup_save_ratio(1);  iter_save_ratio(2) setup_save_ratio(2);iter_save_ratio(3) setup_save_ratio(3) ]

%bar(R)
%bar( [iter_speed; setup_speed] )
%bar(setup_speed)
%bar(iter_speed)

R = [ iter_save_ratio(1) setup_save_ratio(1);  iter_save_ratio(2) setup_save_ratio(2);iter_save_ratio(3) setup_save_ratio(3) ]




my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;


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


%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters fluxFig
%}


%print -dpng fluxFig2





