clc; clear all; close all; 


N = 5*10;                            % Number of cells
n = 5;                             % Number of cells in each block
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
% ensuring that rows sum to zero
A = A - diag(sum(A, 2));
% Scaling matrix
A = A./full(max(max(A)));

FigHandle = figure('Position', [200, 200, 500, 500]);
plotCellData(G, zeros(N*N,1));
outlineCoarseGrid(G,pv);




%{
mrstModule add spe10
layers = 35:35;
[G, ~, rock] = SPE10_setup(layers);
pv = partitionUI(G, [6, 11, ceil(G.cartDims(3)./5)]);
CG = generateCoarseGrid(G,pv); 
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG);
%CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);
CG = setupMexInteractionMapping(CG);
hT = computeTrans(G,rock); 
A = getIncomp1PhMatrix(G, hT);
% ensuring that rows sum to zero
A = A - diag(sum(A, 2));
% Scaling matrix
A = A./full(max(max(A)));
min(min(A))
max(max(A))

plotCellData(G,log10(rock.perm(:,1)))
shading flat, axis equal off, set(gca, 'zdir' ,' reverse'), box on;
colorbar( 'horiz' );
outlineCoarseGrid(G,pv)

%}



iterations = 1000; 
tol = 10^(-12);


[Jbasis1, Jerror1] = getJBasis(CG, A, iterations, tol, 2/3);
[Jbasis2, Jerror2] = getJBasis(CG, A, iterations, tol, 0.95);
[GSbasis, GSerror] = getBLGSBasis(CG, A, iterations, tol, 1);

%[GSbasis, GSerror] = getRedBlackBasis(CG, A, iterations, tol, 1);


AST = min([ min(Jerror1) min(Jerror2) min(GSerror) ]);

disp(['Start Error: ', num2str(Jerror1(1))]);
disp(['AST: ', num2str(AST)]);
disp(['Jerror1: ', num2str(min(Jerror1)-AST)]);
disp(['Jerror2: ', num2str(min(Jerror2)-AST)]);
disp(['GSerror: ', num2str(min(GSerror)-AST)]);

Jerror1 = (Jerror1-AST)/(Jerror1(1)-AST);
Jerror2 = (Jerror2-AST)/(Jerror2(1)-AST);
GSerror = (GSerror-AST)/(GSerror(1)-AST);


F = 1;
for i = 2:length(Jerror2)
    temp = -1;
    for j = 1:length(GSerror)
       if GSerror(j)<=Jerror2(i) 
           temp = j;
           break;
       end
    end
    F(i) = j-1;
end

for i = 2:length(Jerror2)
    F(i) = F(i)/(i-1); 
end

L1 = length(Jerror1);
L2 = length(Jerror2);
L3 = length(GSerror);





lineSize = 2; 
legSize = 15; 
axisSize = 12; 
FigHandle = figure('Position', [200, 800, 500, 500]);
R3 = 1300; 
R1 = R3 + L1-L3; 
R2 = R3 + L2 - L3;
semilogy([0:L3-1 - R3],GSerror(1:L3-R3),'k','LineWidth',lineSize)
hold on
semilogy([0:L1-1-R1],Jerror1(1:L1-R1),'r','LineWidth',lineSize)
semilogy([0:L2-1-R2],Jerror2(1:L2-R2),'b','LineWidth',lineSize)
LEG1 = legend('BLGS','Jacobi, 2/3','Jacobi, 0.95');
set(gca, 'FontSize',axisSize);
set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('\lambda^k', 'FontSize',legSize+5);



set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters error_BLGS_400



L = length(GSerror)-1; 
start = 6; 
lineSize = 2; 
legSize = 15; 
axisSize = 12; 
FigHandle = figure('Position', [1000, 800, 500, 500]);
REMOVE = 0;
semilogx([start:L-REMOVE-1],F(start+1:L-REMOVE),'-','LineWidth',lineSize)
set(gca, 'FontSize',axisSize);
%set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('Fraction', 'FontSize',legSize);


set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters fraction_BLGS_400




%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters fraction_BLGS_25
%}




















%{
Jbasis = full(Jbasis);
GSbasis = full(Jbasis);
RBbasis = full(Jbasis);

AST = min([ min(Jerror) min(GSerror) min(RBerror) ]);

disp(['AST: ', num2str(AST)]);
disp(['Jerror: ', num2str(min(Jerror)-AST)]);
disp(['GSerror: ', num2str(min(GSerror)-AST)]);
disp(['RBerror: ', num2str(min(RBerror)-AST)]);

Jerror = (Jerror-AST)/(Jerror(1)-AST);
GSerror = (GSerror-AST)/(GSerror(1)-AST);
RBerror = (RBerror-AST)/(RBerror(1)-AST);

L = min ([ length(Jerror) length(GSerror) length(RBerror)  ]);
R = L-100; 

lineSize = 2; 
legSize = 15; 
axisSize = 10; 
FigHandle = figure('Position', [200, 800, 800, 800]);
semilogy([0:L-R-1],RBerror(1:L-R),'k','LineWidth',lineSize)
hold on
semilogy([0:L-R-1],Jerror(1:L-R),'r','LineWidth',lineSize)
semilogy([0:L-R-1],GSerror(1:L-R),'b','LineWidth',lineSize)
LEG1 = legend('Red-black','Jacobi, 1','BLGS');
set(gca, 'FontSize',axisSize);
set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('Error', 'FontSize',legSize);

%}












%{
FigHandle = figure('Position', [200, 800, 500, 500]);
plotCellData(G, Jbasis(:,1));
outlineCoarseGrid(G,pv);

FigHandle = figure('Position', [800, 800, 500, 500]);
plotCellData(G, GSbasis(:,1));
outlineCoarseGrid(G,pv);
%}














