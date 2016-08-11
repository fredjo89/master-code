close all; 
read = [ 58.6 53.8 50.1 47.1 43.7 44.1 41.5 35.3 30.3 25.5 22.6 18.9 16.3 ];
makeLocal = [ 0.198 0.908 2.44 4.43 8.42 13.4 20.6 28.1 37.7 45.1 52.8 58.4 63.7];
mainComp = [ 3.63 7.26 12.6 14.8 17.4 16.3 16 14.7 13.8 12.2 11 9.75 8.81];
send = [ 29.5 24.4 18.5 17.1 12 10.7 7.81 8 5.56 5.65 3.61 3.98 3.38];
norm = [ 5.32 9.5 10.52 10.25 10.97 8.761 7.274 7.45 6.78 6.357 5.15 4.693 4.022];

x = [ 10 20 30 40 50 60 80 100 120 140 160 180 200];
subplot(2,1,1);
hold on; 
plot(x,read,'--bo', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); 

plot(x,makeLocal,'--ro', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); 

plot(x,mainComp,'--go', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g'); 

plot(x,send,'--ko', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'); 

plot(x,norm,'--mo', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm'); 
LEG1 = legend('Read from file','Create local data-structure','Main computation','Sending and recieving','Normalization')

fSize = 20; 
legSize = 15; 
xlabel('Number cells along one axis', 'FontSize',fSize)
ylabel('Percent of total runtime', 'FontSize',fSize)


subplot(2,1,2); hold on;
total = read+makeLocal+mainComp+send+norm; 
linear = read+makeLocal; 
iteration = mainComp+send+norm; 
plot(x,total,'--bo', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); 
plot(x,linear,'--ro', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); 
plot(x,iteration,'--go', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g'); 
LEG2 = legend('Sum','Sequential computation','Computation in Jacobi iteration')


xlabel('Number cells along one axis', 'FontSize',fSize)
ylabel('Percent of total runtime', 'FontSize',fSize)
set(LEG1,'FontSize',legSize);
set(LEG2,'FontSize',legSize);


















%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters epsFig
%}



%%Source info
%{

1: Read
2: LocalGrid and matrix 
3: Main comutation in loop 
4: Modification one
5: Sending and recieving
6: Normalization one
7: Normalization two
8: updating basis functions
9: Renormalization


10 X 10 
1: 58.6
2: 0.198
3: 3.63
4: 1.18
5: 29.5
6: 4.42
7: 0.9
8: 0.788
9: 0.825


20 X 20
1: 53.8
2: 0.908
3:7.26
4:2.02
5:24.4
6: 8.37
7: 1.13
8: 1.28
9: 0.907

30 X 30 
1: 50.1
2: 2.44
3: 12.6
4: 2.53
5: 18.5
6: 9.42
7: 1.1
8: 1.94
9: 1.39

40 X 40
1: 47.1
2: 4.43
3: 14.8
4: 2.65
5: 17.1
6: 9.18
7: 1.07
8: 2.25
9: 1.35

50 X 50 
1: 43.7
2: 8.42
3: 17.4
4: 3.32
5: 12
6: 9.93
7: 1.04
8: 2.74
9: 1.52

60 X 60
1: 44.1
2: 13.4
3: 16.3
4: 2.82
5: 10.7
6: 7.94
7: 0.821
8: 2.68
9: 1.32


80 X 80
1: 41.5
2: 20.6
3: 16
4: 3.05
5: 7.81
6: 6.61
7: 0.664
8: 2.63
9: 1.23

100 X 100
1:35.3
2:28.1
3:14.7
4:2.81
5: 8
6: 6.94
7: 0.51
8: 2.48
9: 1.14

120 X 120
30.3
37.7
13.8
2.42
5.56
6.36
0.421
2.35
1.01

140 X 140 
25.5
45.1
12.2
2.14
5.65
6.03
0.327
2.09
0.887

160 X 160 
22.6
52.8
11
2.01
3.61
4.9
0.25
1.98
0.765

180 X 180
18.9
58.4
9.75
1.73
3.98
4.48
0.213
1.8
0.679

200 X 200
16.3
63.7
8.81
1.53
3.38
3.84
0.182
1.6
0.601

%}