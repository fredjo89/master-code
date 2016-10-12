clc; clear all; close all; 

% 20 x 20 blocks
% upscaling factor: 100


x = [
1
16
32
64
128
320
640
1280
1600
]';

T_total = [
32.73
3.0490
1.351
0.517
0.3198
0.1226
0.09537
0.06533
0.0608
];

T_jacobi = [
30.56
2.8780
1.209
0.4319
0.2122
0.06521
0.04978
0.01676
0.01704
];

T_send = [
0
0.0466
0.09091
0.05953
0.09355
0.0527
0.04093
0.04638
0.04052
];

%T_total = T_total(1:7);
%T_jacobi = T_jacobi(1:7);
%T_send = T_send(1:7);


T_norm = T_total - T_jacobi - T_send;



T_j_frac = T_jacobi./T_total

T_s_frac = T_send./T_total

T_n_frac = T_norm./T_total;





my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;




xt = [1:length(x)];

xt = [1 2 3 4 5 6 7 8 9.3];

newY = [ T_s_frac T_n_frac T_j_frac   ]

FigHandle = figure('Position', [1200, 200, 80*length(x), 400]);
bar(newY,'stacked');  %# Create a stacked histogram
set(gca, 'XTickLabel',[])                    
set(gca, 'XTick', x)
set(gca,'ytick',[0:0.2:1]);
xlabel('Number of cores');
ylabel('Fraction');
LEG1 = legend('Message-passing','Normalizing', 'Smoothing');
colormap([my_red_1; my_blue_1; my_green_1])
set(gca,'fontsize',16)


set(gca, 'XTickLabel',[])                    
%xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', x)
set(gca,'ytick',[0:0.2:2]);


%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters epsFig_medium
%}






