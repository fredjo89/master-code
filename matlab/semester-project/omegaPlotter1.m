clc; clear all; close all; 

%% Case for homo rock, 120 x 120, 16 blocks

omega = [
    0.05
0.1
0.15
0.2
0.25
0.3
0.35
0.4
0.45
0.5
0.55
0.6
0.65
0.7
0.75
0.8
0.85
0.9
0.91
0.92
0.93
0.94
0.95
0.96
0.97
0.98
0.981
0.982
0.983
0.984
0.985
0.986
0.987
0.988
0.989
0.99
0.991
0.992
0.993
0.994
0.995
0.9955
0.996
0.9965
0.997
0.99725
0.9975
0.99775
0.998
0.9982
0.9984
0.9986
0.9988
0.999
0.9992
0.9994
0.9996
0.9998
1
];


it = [
2020
1640
1346
1147
1004
897
813
745
689
643
602
568
537
511
487
465
446
429
425
422
419
416
413
410
407
404
403
403
403
402
402
402
402
401
401
401
401
400
400
401
401
402
404
407
411
413
417
421
426
431
439
447
457
471
487
507
536
570
608
];

% plotting
fSize = 60; 
axisSize = 40; 
Lwidth = 15; 
dotSize = 30; 
figure();
plot(omega,it, 'LineWidth', Lwidth); 


hold off;
%LEG1 = legend(' 20 x 20 ');
set(gca, 'FontSize',axisSize);
xlabel('Relaxation factor', 'FontSize',fSize)
ylabel('Iterations', 'FontSize',fSize)
%set(LEG1,'FontSize',legSize);
axis([0 1 0 2000])


figure(); 
dotSize = 30; 
Lwidth = 4; 
N = 30; 
plot(omega(N:length(omega)),it(N:length(omega)), '--ko', 'LineWidth', Lwidth,'MarkerSize', dotSize, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); 
axis([0.9883 0.999 390 450])

x1 = 0.992; 
x2 = 0.993;
min = 400; 
hold on
plot(x1,min, 'ro', 'LineWidth', Lwidth,'MarkerSize', dotSize, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); 
plot(x2,min, 'ro', 'LineWidth', Lwidth,'MarkerSize', dotSize, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); 
set(gca, 'FontSize',axisSize);
xlabel('Relaxation factor', 'FontSize',fSize)
ylabel('Iterations', 'FontSize',fSize)


%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters epsFig
%}

