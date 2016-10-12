clc; clear all; close all; 

%{
N = 100;                            % Number of cells
n = 10;                             % Number of cells in each block
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
%}


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
% Scaling matrix
A = A./full(max(max(A)));

plotCellData(G,log10(rock.perm(:,1)))
shading flat, axis equal off, set(gca, 'zdir' ,' reverse'), box on;
colorbar( 'horiz' );
outlineCoarseGrid(G,pv)



iterations = 20000000; 
tol = 10^(-8);
disp('Jacobi 1')
tic
Jerror1 = getJBasis_compare2(CG, A, iterations, tol,  2/3);
toc
disp('Jacobi 2')
tic
Jerror2 = getJBasis_compare2(CG, A, iterations, tol,  0.95);
toc
disp('RED-BLACK 1')
tic
RBerror = RedBlack_compare(CG, A, iterations,tol, 1);
toc
disp('BLGS 1')
tic
[BLGSbasis, BLGSerror] = getBLGSBasis(CG, A, iterations, tol, 1);
toc

AST = min([ min(Jerror1) min(Jerror2) min(RBerror) min(BLGSerror) ]);

Jerror1(length(Jerror1))
Jerror2(length(Jerror2))
RBerror(length(RBerror))
BLGSerror(length(BLGSerror))

disp(['Start Error: ', num2str(Jerror1(1))]);
disp(['AST: ', num2str(AST)]);
disp(['Jerror1: ', num2str(min(Jerror1)-AST)]);
disp(['Jerror2: ', num2str(min(Jerror2)-AST)]);
disp(['RBerror: ', num2str(min(RBerror)-AST)]);
disp(['BLGSerror: ', num2str(min(BLGSerror)-AST)]);

Jerror1 = (Jerror1-AST)/(Jerror1(1)-AST);
Jerror2 = (Jerror2-AST)/(Jerror2(1)-AST);
RBerror = (RBerror-AST)/(RBerror(1)-AST);
BLGSerror = (BLGSerror-AST)/(BLGSerror(1)-AST);

F = 1;
for i = 2:length(Jerror2)
    temp = -1;
    for j = 1:length(RBerror)
       if RBerror(j)<=Jerror2(i) 
           temp = j;
           break;
       end
    end
    F(i) = j-1;
end

for i = 2:length(Jerror2)
    F(i) = F(i)/(i-1); 
end

F2 = 1;
for i = 2:length(Jerror2)
    temp = -1;
    for j = 1:length(BLGSerror)
       if BLGSerror(j)<=Jerror2(i) 
           temp = j;
           break;
       end
    end
    F2(i) = j-1;
end

for i = 2:length(Jerror2)
    F2(i) = F2(i)/(i-1); 
end






L1 = length(Jerror1);
L2 = length(Jerror2);
L3 = length(RBerror);
L4 = length(BLGSerror);




my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;



lineSize = 1.5; 
legSize = 10; 
axisSize = 10; 
FigHandle = figure('Position', [1200, 200, 13*29, 11.5*29]);
L = 500;
loglog([0:L-1],Jerror1(1:L),'k','LineWidth',lineSize)
hold on
loglog([0:L-1],Jerror2(1:L),'r','LineWidth',lineSize)
loglog([0:L-1],RBerror(1:L),'color',my_green_1,'LineWidth',lineSize)
loglog([0:L-1],BLGSerror(1:L),'color',my_blue_1,'LineWidth',lineSize)
LEG1 = legend('Jacobi, 2/3','Jacobi, 0.95','Red-black','BLGS');
xlabel('k','FontSize',15);
ylabel('\lambda', 'FontSize',15);
%axis([0,1000,-1.8, 0 ]);
%set(gca,'xtick',x);
%set(gca,'ytick',[0:0.1:1]);
%set(LEG1);




L = min([length(RBerror) length(BLGSerror) length(Jerror2)]); 
L = L - 1; 
start = 3; 
lineSize = 0.5; 
legSize = 15; 
axisSize = 10; 
FigHandle = figure('Position', [400, 200, 13*29, 11.5*29]);
REMOVE = L-500;
%REMOVE = 1000;
semilogx([start:L-REMOVE-1],F(start+1:L-REMOVE),'--o','color',my_green_1,'LineWidth',lineSize,'MarkerSize', 4)
hold on;
semilogx([start:L-REMOVE-1],F2(start+1:L-REMOVE),'--o','color',my_blue_1,'LineWidth',lineSize,'MarkerSize', 4)
%set(gca, 'FontSize',axisSize);
%set(LEG1,'FontSize',legSize);
xlabel('k','FontSize',15);
ylabel('s', 'FontSize',15);
LEG1 = legend('Red-black','BLGS');







%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters fraction_homo
%}



