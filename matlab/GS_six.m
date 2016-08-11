%% 1D example 
clc; clear all; close all; 


N = 50;                              % Number of cells
n = 5;                                 % Number of cells in each block

M = ceil(N / n);                    % Number of blocks

G = cartGrid([N N]); 
G = computeGeometry(G);
rock = makeRock(G, 1, 1);
pv = partitionUI(G, [M, M]);
CG = generateCoarseGrid(G,pv); 
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG);
%CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);
CG = setupMexInteractionMapping(CG);
hT = computeTrans(G,rock); 
A = getIncomp1PhMatrix(G, hT);
gravity reset off; 
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
initState = initResSol(G,0.0, 1.0); 



iterations = 10; 
tol = 10^(-4); 
w = 1;

Jbasis_CONV = getJBasis(CG, A, 20000, tol, 0.95);


Jerror1 = getJBasis_compare(CG, A, iterations,  2/3,Jbasis_CONV.B );
Jerror2 = getJBasis_compare(CG, A, iterations,  0.95,Jbasis_CONV.B );
Jerror3 = getJBasis_compare(CG, A, iterations,  1,Jbasis_CONV.B );
RBerror = GSbasisRedBlack2(CG, A, iterations, -1, 1, Jbasis_CONV.B);


x = [0:iterations]; 


F = 0;
for i = 2:iterations+1
    temp = -1;
    for j = 1:iterations+1
       if RBerror(2,j)<=Jerror3(2,i) 
           temp = j;
           break;
       end
    end
    F(i) = j-1;
end

for i = 2:iterations+1
    F(i) = F(i)/x(i); 
end

F(1) = 1; 
figure(); 
semilogx(x,F)









lineSize = 2; 
legSize = 15; 
axisSize = 10; 
figure();
semilogy(x,Jerror1(2,:),'b', 'LineWidth',lineSize);
hold on
semilogy(x,Jerror2(2,:),'r', 'LineWidth',lineSize);
semilogy(x,Jerror3(2,:),'g', 'LineWidth',lineSize);
semilogy(x,RBerror(2,:),'k', 'LineWidth', lineSize);
LEG1 = legend('Jacobi, 2/3','Jacobi, 0.95','Jacobi, 1', ' Red-black');
set(gca, 'FontSize',axisSize);
set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('Error', 'FontSize',legSize);





lineSize = 2; 
legSize = 15; 
axisSize = 10; 
figure();
semilogx(x,F,'b', 'LineWidth',lineSize);
LEG1 = legend('Jacobi, 2/3','Jacobi, 0.95','Jacobi, 1', ' Red-black');
set(gca, 'FontSize',axisSize);
set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('Error', 'FontSize',legSize);

























%{
figure();
plot(temp2(2,:),'b', 'LineWidth',lineSize);
hold on
plot(temp2(4,:),'r', 'LineWidth', lineSize);
LEG1 = legend('Jacobi', ' Red-black');
set(gca, 'FontSize',axisSize);
set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('Convergence factor', 'FontSize',legSize);
%}

%{
figure();
semilogy(RBerror(2,:),'--')
hold on
semilogy(Jerror(2,:),'--')



figure(); 
hold on; 
plot(temp2(1,:),'--')
%plot(temp2(2,:),'--')
plot(temp2(3,:),'--')
%plot(temp2(4,:),'--')
%}


%{
v = [ 3 4 1 2];
FigHandle = figure('Position', [700, 800, 800, 800]);
for i = 1:M*M
    subplot(M,M,v(i))
    plotCellData(G,full(Jbasis(:,i)));
    outlineCoarseGrid(G,pv,'r')
end

FigHandle = figure('Position', [700, 800, 800, 800]);
for i = 1:M*M
    subplot(M,M,v(i))
    plotCellData(G,full(GSbasis(:,i)));
    outlineCoarseGrid(G,pv,'r')
end
%}

%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters inffactor-spe35
%}



