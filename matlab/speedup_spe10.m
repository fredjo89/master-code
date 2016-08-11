clc; clear all; close all; 

t_iter = [ 184.5 95.56 64.72 49.48 40.84 35.2 31.84 29.51 26.58 23.54 21.52 19.91 18.05 17.17 16.1 15.04 13.5 12.31 11.13 10.17 9.62 8.68 8.339 7.78 6.922 6.498 5.884 5.913 5.544 4.914 4.846 4.519 4.21 3.744 3.582 3.97];


t_setup = [0.548 0.7546 0.7616 1.053 0.7016 1.57 1.37 1.303 0.9134 1.352 1.301 1.233 1.199 1.159 0.9283 0.9862 1.364 1.254 1.248 1.215 1.1 1.136 0.9476 1.049 1.226 0.984 1.145 0.9811 0.9251 1.086 0.828 0.7787 0.9068 1.061 1.03 0.7988];

t_send = [0 0.1164 0.2794 0.347 0.6421 0.435 0.6582 0.6895 0.7711 0.523 0.4511 1.283 0.4009 1.151 0.7546 1.752 0.8706 1.143 1.804 1.796 0.8551 0.64 1.398 1.287 0.9387 1.222 0.7689 1.345 1.383 1.066 1.119 1.066 1.325 0.8292 1.178 1.366];



size(t_iter)

t_speedup = 1./(t_iter./t_iter(1));




x = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 18 20 22 24 26 28 30 32 36 39 42 45 48 52 56 60 64 70 75 80];

figure(); 
plot(x,t_setup,'*')

figure(); 
plot(x,t_send./t_iter,'*')




figure();
hold on; 
plot(x,t_speedup,'--bo', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); 
plot(x,x, 'LineWidth', 1.5)

LEG1 = legend('yolo');

fSize = 20; 
legSize = 15; 
xlabel('Number of processors', 'FontSize',fSize)
ylabel('Speedup', 'FontSize',fSize)
axis([1,80,1,80])
set(gca,'xtick',0:5:80, 'FontSize',fSize)
set(gca,'ytick',0:5:80)
set(LEG1,'FontSize',legSize);



set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters epsFig

%}



