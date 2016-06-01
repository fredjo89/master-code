clc; clear all; close all; 


FC=[8:4:40];           % Fine cells
CB = 2;                % Coarse blocks
l = length(FC); 
 
iterations = 1000; 
tol = 0.001; 



for (i = 1:l)
    FC(i)
    output(:,i) = homoRun(FC(i),CB, iterations, tol);
end

maxIter = max(output(1,:))



%% Plotting
markSize = 8; 

figure(); 
hold on
plot(FC,output(1,:),'--bo', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
plot(FC,output(2,:),'--ro', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); 
plot(FC,output(3,:),'--go', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g'); 
LEG1 = legend('Jacobi, 0.95', 'Gauss-Seidel', 'New GS');

fSize = 10; 
legSize = 10; 
xlabel('Number of cells along one axis', 'FontSize',fSize)
ylabel('Iterations', 'FontSize',fSize)
axis([FC(1),FC(l),0,maxIter])
set(gca, 'FontSize',fSize)
%set(gca,'ytick',0:10:maxIter)
set(LEG1,'FontSize',legSize);


figure();
hold on
plot(FC,output(4,:),'--bo', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
plot(FC,output(5,:),'--ro', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); 
plot(FC,output(6,:),'--go', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g'); 
LEG1 = legend('Jacobi, 0.95', 'Gauss-Seidel', 'New GS');

fSize = 10; 
legSize = 10; 
xlabel('Number of cells along one axis', 'FontSize',fSize)
ylabel('Inf-norm', 'FontSize',fSize)
axis([FC(1),FC(l),0,.7])
%set(gca,'xtick',FC, 'FontSize',fSize)
set(gca, 'FontSize',fSize)
set(gca,'ytick',0:0.1:1)
set(LEG1,'FontSize',legSize);

%{
figure();
hold on
plot(FC,output(7,:),'--bo', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
plot(FC,output(8,:),'--ro', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); 
plot(FC,output(9,:),'--go', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g'); 
LEG1 = legend('Jacobi, 2/3','Jacobi, 0.95', 'Gauss-Seidel');

fSize = 10; 
legSize = 10; 
xlabel('Number of cells along one axis', 'FontSize',fSize)
ylabel('Two-norm', 'FontSize',fSize)
axis([FC(1),FC(l),0,1])
set(gca, 'FontSize',fSize)
set(gca,'ytick',0:0.1:1)
set(LEG1,'FontSize',legSize);

%}

%{

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

%}

%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters 4inf
%}



