clc; clear all; close all; 

four = [ 0.00970 0.00668 0.00664 0.00593] ; 
onesix = [0.0527 0.0273 0.0211 0.0184] ; 
sixfour = [ 0.239 0.123 0.0858 0.0686]; 
twotwofive = [ 0.764 0.408 0.299 0.263]; 
oneonefour = [0.60 0.304 0.220 0.199];

s_four = 1./(four./four(1));
s_onesix = 1./(onesix./onesix(1));
s_sixfour = 1./(sixfour./sixfour(1));
s_twotwofive = 1./(twotwofive./twotwofive(1));
s_oneonefour = 1./(oneonefour./oneonefour(1)); 

eighty = [ 0.109 0.0563 0.0398 0.0325];
onetwenty = [ 0.239 0.122 0.0846 0.0684]; 
onesixty = [ 0.416 0.210 0.146 0.120]; 
twohundred = [ 0.678 0.351 0.251 0.215];

s_eighty = 1./(eighty./eighty(1)); 
s_onetwenty = 1./(onetwenty./onetwenty(1)); 
s_onesixty = 1./(onesixty./onesixty(1)); 
s_twohundred = 1./(twohundred./twohundred(1)); 

x = linspace(1,4,4);


subplot(1,2,1);
hold on; 
plot(x,s_four,'--bo', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); 
plot(x,s_onesix,'--ro', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); 
plot(x,s_sixfour,'--go', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g'); 
%plot(x,s_twotwofive,'--ko', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'); 
plot(x,s_oneonefour,'--ko', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'); 
plot(x,x, 'LineWidth', 1.5)

LEG1 = legend('4 blocks, 30 x 30 fine cells',...
    '16 blocks, 60 x 60 fine cells','64 blocks, 120 x 120 fine cells'...
    ,'144 blocks, 180 x 180 fine cells ', 'Perfect speedup');

fSize = 20; 
legSize = 15; 
xlabel('Number of processors', 'FontSize',fSize)
ylabel('Speedup', 'FontSize',fSize)
axis([1,4,1,4])
set(gca,'xtick',0:4, 'FontSize',fSize)
set(gca,'ytick',0:4)
set(LEG1,'FontSize',legSize);

%% Second plot
subplot(1,2,2);
hold on; 
plot(x,s_eighty,'--bo', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); 

plot(x,s_sixfour,'--ro', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); 

plot(x,s_onesixty,'--go', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g'); 

plot(x,s_twohundred,'--ko', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'); 

plot(x,x, 'LineWidth', 1.5)


LEG1 = legend('80 x 80 fine cells','120 x 120 fine cells',...
    '160 x 160 fine cells','200 x 200 fine cells', 'Perfect speedup');

fSize = 20; 
legSize = 15; 
xlabel('Number of processors', 'FontSize',fSize)
ylabel('Speedup', 'FontSize',fSize)
axis([1,4,1,4])
set(gca,'xtick',0:4, 'FontSize', fSize)
set(gca,'ytick',0:4)
set(LEG1,'FontSize',legSize);



set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters epsFig



