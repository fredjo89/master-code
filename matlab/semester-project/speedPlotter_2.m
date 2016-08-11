close all; 
read = [ 
58.6
54
53.9
50.5
45.5
46.7
42.5
36.6
30.8
27.2
22.4
19.1
16.2
    ];
makeLocal = [ 
0.193
0.799
2.3
4.32
8.86
12.2
20.2
29.5
37.7
46
53.1
58.9
63.5
    ];
mainComp = [ 
2.91
6.89
11.3
14.3
17.1
15.6
16.1
15.5
13.9
12.5
11.1
9.81
8.77
    ];
send = [ 
30.3
25.5
17.7
14.9
11.7
11.1
6.96
6.07
6.88
4.19
3.66
3.65
3.65
    ];
norm = [ 
5.495
8.82
9.43
9.84
9.89
8.023
7.352
5.758
4.934
4.986
4.894
4.4323
4.12
    ];

reNorm = [
0.752
0.871
1.27
1.33
1.38
1.23
1.23
1.16
1.02
0.875
0.769
0.672
0.596
];

update = [
    0.713
1.29
1.83
2.21
2.68
2.48
2.67
2.6
2.35
2.15
2.04
1.8
1.61
];

modOne = [
1.03
1.92
2.34
2.6
2.96
2.62
3.05
2.82
2.45
2.13
2.01
1.75
1.55
];

modifyAndUpdate = modOne + update; 


x = [ 10 20 30 40 50 60 80 100 120 140 160 180 200];
subplot(2,1,1);
hold on; 
plot(x,read,'--bo', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); 

plot(x,makeLocal,'--ro', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); 

plot(x,mainComp,'--go', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g'); 

plot(x,send,'--ko', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'); 

plot(x,norm,'--mo', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm'); 

%plot(x,reNorm,'--yo', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'y'); 

%plot(x,modifyAndUpdate,'--co', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c'); 




LEG1 = legend('Read from file','Create local data-structure','Main computation','Sending and recieving','Normalization')
fSize = 20; 
legSize = 15; 
axisSize = 20;
xlabel('Number cells along one axis', 'FontSize',fSize)
ylabel('Percent of total runtime', 'FontSize',fSize)
set(gca, 'FontSize',axisSize);

subplot(2,1,2); hold on;
linear = read+makeLocal; 
iteration = mainComp+send+norm+reNorm+modifyAndUpdate; 
%total = linear + iteration; 
%plot(x,total,'--bo', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); 
plot(x,linear,'--bo', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); 
plot(x,iteration,'--ro', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); 
LEG2 = legend('Sequential computation','Computation in Jacobi iteration')
xlabel('Number cells along one axis', 'FontSize',fSize)
ylabel('Percent of total runtime', 'FontSize',fSize)
set(LEG1,'FontSize',legSize);
set(LEG2,'FontSize',legSize);
axis([0 200 0 100]);
set(gca, 'FontSize',axisSize);




clear all; 
x = [ 10 20 30 40 50 60 80 100 120 140 160 180 200];
MPI_TOTAL = [
0.00775
0.0094
0.0119
0.0168
0.0214
0.0289
0.0442
0.0705
0.112
0.17
0.25
0.361
0.51
];

OMP_TOTAL = [
0.00645
0.0109
0.017
0.0254
0.0328
0.0442
0.0705
0.1035
0.146
0.197
0.2505
0.315
0.401
];
figure();
hold on; 
plot(x,MPI_TOTAL,'--bo', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); 

plot(x,OMP_TOTAL,'--ro', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); 


fSize = 20; 
legSize = 15; 
axisSize = 20; 
xlabel('Number cells along one axis', 'FontSize',fSize)
ylabel('RunTime in seconds', 'FontSize',fSize)
LEG = legend('MPI-code','openMP-code')
set(LEG,'FontSize',legSize);
set(gca, 'FontSize',axisSize);










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
58.6
0.193
2.91
1.03
30.3
4.67
0.825
0.713
0.752
MPI TOTAL:  0.00775 
OMPTOTAL:   0.00645

20 X 20
54
0.799
6.89
1.92
25.5
7.72
1.1
1.29
0.871
MPI TOTAL:  0.0094
OMPTOTAL:   0.0109

30 X 30 
53.9
2.3
11.3
2.34
17.7
8.35
1.08
1.83
1.27
MPI TOTAL:  0.0119
OMPTOTAL:   0.0170

40 X 40
50.5
4.32
14.3
2.6
14.9
8.76
1.08
2.21
1.33
MPI TOTAL:  0.0168
OMPTOTAL:   0.0254

50 X 50 
45.5
8.86
17.1
2.96
11.7
8.88
1.01
2.68
1.38
MPI TOTAL:  0.0214
OMPTOTAL:   0.0328

60 X 60
46.7
12.2
15.6
2.62
11.1
7.13
0.893
2.48
1.23
MPI TOTAL:  0.0289
OMPTOTAL:   0.0442

80 X 80
42.5
20.2
16.1
3.05
6.96
6.68
0.672
2.67
1.23
MPI TOTAL:  0.0442
OMPTOTAL:   0.0705

100 X 100
36.6
29.5
15.5
2.82
6.07
5.25
0.508
2.6
1.16
MPI TOTAL:  0.0705
OMPTOTAL:   0.1035

120 X 120
30.8
37.7
13.9
2.45
6.88
4.55
0.384
2.35
1.02
MPI TOTAL:  0.112
OMPTOTAL:   0.146

140 X 140 
27.2
46
12.5
2.13
4.19
4.69
0.299
2.15
0.875
MPI TOTAL:  0.170    
OMPTOTAL:   0.197

160 X 160 
22.4
53.1
11.1
2.01
3.66
4.64
0.254
2.04
0.769
MPI TOTAL:  0.250
OMPTOTAL:   0.2505

180 X 180
19.1
58.9
9.81
1.75
3.65
4.11
0.213
1.8
0.672
MPI TOTAL:  0.361
OMPTOTAL:   0.315

200 X 200
16.2
63.5
8.77
1.55
3.65
3.95
0.17
1.61
0.594
MPI TOTAL:  0.51
OMPTOTAL :  0.401




%}