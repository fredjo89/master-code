clc; clear all; close all; 


FC=9;                   % Fine cells
CB = 3;                 % Coarse blocks
iterations = [1:1:20]; 
tol = -1; 
l = length(iterations); 


 




for (i = 1:l)
    iterations(i)
    %output(:,i) = homoRun_relaxedGS(FC,CB, iterations(i), tol);
    output(:,i) = homoRun(FC,CB, iterations(i), tol);
end





maxIter = max(output(1,:))



%% Plotting

lowestINF = min([min(min(output(4,:))), min(min(output(5,:))), min(min(output(6,:)))])
largestINF = max([max(max(output(4,:))), max(max(output(5,:))), max(max(output(6,:)))])

lowestTWO = min([min(min(output(7,:))), min(min(output(8,:))), min(min(output(9,:)))])
largestTWO = max([max(max(output(7,:))), max(max(output(8,:))), max(max(output(9,:)))])

markSize = 8; 

jConv = ones(1,l)*output(4,l) ;


subplot(2,2,1)
hold on
title(['Blocks: ' num2str(CB) '  Cells: ' num2str(FC)  ])
plot(iterations,output(4,:),'bo', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
plot(iterations,output(5,:),'ro', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); 
plot(iterations,output(6,:),'go', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g'); 
%plot(iterations,jConv,'c', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c'); 
LEG1 = legend('GS Small', ' GS Mid', 'GS Large', 'Smallest Jacobi error');

fSize = 20; 
legSize = 20; 
xlabel('Iterations', 'FontSize',fSize)
ylabel('Inf-norm', 'FontSize',fSize)
axis([iterations(1),iterations(l),max(lowestINF, 0 ),min(largestINF, 1)])
%set(gca,'xtick',FC, 'FontSize',fSize)
set(gca, 'FontSize',fSize)
%set(gca,'ytick',0:0.05:1)
set(LEG1,'FontSize',legSize);


diff_1 = output(4,:)-output(6,:);

subplot(2,2,2)
hold on
title(['Blocks: ' num2str(CB) '  Cells: ' num2str(FC)  ])
plot(iterations,output(4,:)-output(6,:),'bo', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');


xlabel('Iterations', 'FontSize',fSize)
ylabel('Inf-norm difference', 'FontSize',fSize)
axis([iterations(1),iterations(l),min(min(diff_1), 0),max(diff_1)])
%set(gca,'xtick',FC, 'FontSize',fSize)
set(gca, 'FontSize',fSize)
%set(gca,'ytick',min(min(diff_1), 0):0.01:1)
set(LEG1,'FontSize',legSize);

ax = gca;
ax.XAxisLocation = 'origin';






subplot(2,2,3);
hold on
title(['Blocks: ' num2str(CB) '  Cells: ' num2str(FC)  ])
plot(iterations,output(7,:),'bo', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
plot(iterations,output(8,:),'ro', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); 
plot(iterations,output(9,:),'go', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g'); 
LEG1 = legend('GS Small', ' GS Mid', 'GS Large', 'Smallest Jacobi error');


xlabel('Iterations', 'FontSize',fSize)
ylabel('Two-norm', 'FontSize',fSize)
axis([iterations(1),iterations(l),max(lowestTWO, 0 ),min(largestTWO, 1)])
set(gca, 'FontSize',fSize)
%set(gca,'ytick',0:0.05:1)
set(LEG1,'FontSize',legSize);


diff_1 = output(7,:)-output(9,:);

subplot(2,2,4);
hold on
title(['Blocks: ' num2str(CB) '  Cells: ' num2str(FC)  ])
plot(iterations,diff_1,'bo', 'LineWidth', 1, 'MarkerSize', markSize, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');


xlabel('Iterations', 'FontSize',fSize)
ylabel('Two-norm difference', 'FontSize',fSize)
axis([iterations(1),iterations(l),min(min(diff_1), 0),max(diff_1)])
%set(gca,'xtick',FC, 'FontSize',fSize)
set(gca, 'FontSize',fSize)
%set(gca,'ytick',min(min(diff_1), 0):0.01:max(diff_1))
set(LEG1,'FontSize',legSize);

ax = gca;
ax.XAxisLocation = 'origin';




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
print -dpdf -painters largeHetero_8_60
%}


counter1 = 0; 
counter2 = 0; 
for i=1:l
    temp1 = output(7,i);
    temp2 = output(9,i);
    if abs(temp2)>2
       counter1=counter1 +1;
    end
    if abs(temp1)>2
        counter2 = counter2+1; 
    end
end

counter1
counter2

