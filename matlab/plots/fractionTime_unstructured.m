clc; clear all; close all; 

% 20 x 20 blocks
% upscaling factor: 100


T_total = [
51.2
25.25
13.23
7.856
4.428
2.956
2.257
3.082
3.811
];

T_jacobi = [
37.91
19.05
9.747
5.322
2.729
1.313
0.5487
0.2052
0.2021
];

T_send = [
0
0.2308
0.4563
0.8975
0.9004
1.205
1.514
2.797
3.531
];

%T_total = T_total(1:7);
%T_jacobi = T_jacobi(1:7);
%T_send = T_send(1:7);


T_norm = T_total - T_jacobi - T_send;



T_j_frac = T_jacobi./T_total

T_s_frac = T_send./T_total

T_n_frac = T_norm./T_total;

T_s_frac(6)

x = [
1
2
4
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
%FigHandle = figure('Position', [1200, 200, 500, 500]);
FigHandle = figure('Position', [1200, 200, 13*29, 11.5*29]);
hold on; 
plot(log2(x),T_j_frac,'--o','Color',my_blue_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_blue_1, 'MarkerFaceColor', my_blue_1); 
plot(log2(x),T_s_frac,'--o','Color',my_green_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1); 
plot(log2(x),T_n_frac,'--o','Color',my_red_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_red_1, 'MarkerFaceColor', my_red_1);
yl = get(gca, 'YLim');
str = cellstr( num2str(x(:)) );      %# format x-ticks as 2^{xx}
hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
    'Interpreter','tex', ...                   %# specify tex interpreter
    'VerticalAlignment','top', ...             %# v-align to be underneath
    'HorizontalAlignment','center' ); 


LEG1 = legend('Smoothing','Sending', 'Normalizing');
xlabel('Number of cores');
ylabel('Fraction');
axis([0,6,0,1]);
set(gca,'xtick',x);
set(gca,'ytick',[0:0.2:1]);
set(LEG1);


%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters unstuctured_fraction
%}







newY = [ T_s_frac T_n_frac T_j_frac   ];

FigHandle = figure('Position', [1200, 200, 80*length(x), 400]);
bar(newY,'stacked');  %# Create a stacked histogram
set(gca, 'XTickLabel',[])        ;            
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', x)
set(gca,'ytick',[0:0.2:1]);
set(LEG1);
xlabel('Number of cores');
ylabel('Fraction');
LEG1 = legend('Message-passing','Normalizing', 'Smoothing');
colormap([my_red_1; my_blue_1; my_green_1]);
set(gca,'fontsize',17);

set(gca,'XLim',[0.5 12.5],'YLim',[0 1]);


set(gca,'XLim',[0 10],'YLim',[0 1]);



