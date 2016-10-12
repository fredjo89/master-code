clc; clear all; close all; 

% 20 x 20 blocks
% upscaling factor: 100


T_total = [
7.571
3.76
2.12
1.381
0.9145
0.8853
];

T_jacobi = [
5.676
2.785
1.338
0.5821
0.2781
0.201
];

T_send = [
0
0.1122
0.3164
0.5049
0.4767
0.5928
];

%T_total = T_total(1:7);
%T_jacobi = T_jacobi(1:7);
%T_send = T_send(1:7);


T_norm = T_total - T_jacobi - T_send;



T_j_frac = T_jacobi./T_total

T_s_frac = T_send./T_total

T_n_frac = T_norm./T_total;


x = [
8
16
32
64
128
160   
]';


my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;
my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;
my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;


xt = [0 1 2 3 4 5 6 7 8  ];



%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters unstuctured_fraction_hybrid
%}







newY = [ T_s_frac T_n_frac T_j_frac   ]

FigHandle = figure('Position', [1200, 200, 80*length(x), 400]);
bar(newY,'stacked');  %# Create a stacked histogram
set(gca, 'XTickLabel',[])                    
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', x)
set(gca,'ytick',[0:0.2:1]);
LEG1 = legend('Message-passing','Normalizing', 'Smoothing');
set(LEG1);
xlabel('Number of cores');
ylabel('Fraction');
colormap([my_red_1; my_blue_1; my_green_1])
set(gca,'fontsize',17)

set(gca,'XLim',[0.5 12.5],'YLim',[0 1])


set(gca,'XLim',[0 7],'YLim',[0 1])



