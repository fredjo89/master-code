%% 1D example 
clc; clear all; close all; 


N = 25;                            % Number of cells
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



%{
mrstModule add spe10
layers = 85:85;
[G, ~, rock] = SPE10_setup(layers);
pv = partitionUI(G, [6, 11, ceil(G.cartDims(3)./5)]);
CG = generateCoarseGrid(G,pv); 
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG);
%CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);
CG = setupMexInteractionMapping(CG);
hT = computeTrans(G,rock); 
A = getIncomp1PhMatrix(G, hT);

plotCellData(G,log10(rock.perm(:,1)))
shading flat, axis equal off, set(gca, 'zdir' ,' reverse'), box on;
colorbar( 'horiz' );
outlineCoarseGrid(G,pv)
%}




iterations = 10000; 

Jerror1 = getJBasis_compare2(CG, A, iterations,  2/3);
Jerror2 = getJBasis_compare2(CG, A, iterations,  0.95);
RBerror = GSbasisRedBlack3(CG, A, iterations, 1);

AST = min([ min(Jerror1) min(Jerror2) min(RBerror) ]);

disp(['Start Error: ', num2str(Jerror1(1))]);
disp(['AST: ', num2str(AST)]);
disp(['Jerror1: ', num2str(min(Jerror1))]);
disp(['Jerror2: ', num2str(min(Jerror2))]);
disp(['RBerror: ', num2str(min(RBerror))]);

Jerror1=Jerror1/Jerror1(1);
Jerror2=Jerror2/Jerror2(1);
RBerror=RBerror/RBerror(1);

x = [0:iterations]; 
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
    F(i) = F(i)/x(i); 
end













L = length(RBerror)-1; 
start = 0; 
lineSize = 2; 
legSize = 15; 
axisSize = 10; 
figure();
plot([start:L],RBerror(start+1:L+1),'k','LineWidth',lineSize)
hold on;
plot([start:L],Jerror1(start+1:L+1),'r','LineWidth',lineSize)
plot([start:L],Jerror2(start+1:L+1),'b','LineWidth',lineSize)
LEG1 = legend('Red-black','Jacobi, 2/3','Jacobi, 0.95');
set(gca, 'FontSize',axisSize);
set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('Error', 'FontSize',legSize);


L = length(RBerror)-1; 
start = 2; 
lineSize = 1; 
legSize = 15; 
axisSize = 10; 
figure();
REMOVE = 0;
plot([start:L-REMOVE-1],F(start+1:L-REMOVE),'--o','LineWidth',lineSize)
set(gca, 'FontSize',axisSize);
set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('Fraction', 'FontSize',legSize);


%{
AST = min([min(smallStep) min(Jerror1) min(Jerror2) min(Jerror3) min(RBerror) ]);

disp(['AST: ', num2str(AST)]);
disp(['smallStep: ', num2str(min(smallStep))]);
disp(['Jerror1: ', num2str(min(Jerror1))]);
disp(['Jerror2: ', num2str(min(Jerror2))]);
disp(['Jerror3: ', num2str(min(Jerror3))]);
disp(['RBerror: ', num2str(min(RBerror))]);



L1 = length(Jerror1);
L2 = length(Jerror2);
L3 = length(Jerror3);
L4 = length(RBerror);
L = min([L1 L2 L3]); 


Jerror1 = (Jerror1-AST)/(Jerror1(1)-AST);
Jerror2 = (Jerror2-AST)/(Jerror2(1)-AST);
Jerror3 = (Jerror3-AST)/(Jerror3(1)-AST);
RBerror = (RBerror-AST)/(RBerror(1)-AST);


x = [0:iterations]; 
F = 1;
for i = 2:L2
    temp = -1;
    for j = 1:L4
       if RBerror(j)<=Jerror2(i) 
           temp = j;
           break;
       end
    end
    F(i) = j-1;
end

for i = 2:L2
    F(i) = F(i)/x(i); 
end



lineSize = 2; 
legSize = 15; 
axisSize = 10; 
figure();
REMOVE = 15; 
semilogy([0:REMOVE-1],RBerror(1:REMOVE),'k','LineWidth',lineSize)
hold on
semilogy([0:REMOVE-1],Jerror1(1:REMOVE),'r')
semilogy([0:REMOVE-1],Jerror2(1:REMOVE),'b','LineWidth',lineSize)
semilogy([0:REMOVE-1],Jerror3(1:REMOVE),'g')
LEG1 = legend('Red-black','Jacobi, 2/3','Jacobi, 0.95', 'Jacobi, 1');
set(gca, 'FontSize',axisSize);
set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('Error', 'FontSize',legSize);



lineSize = 1; 
legSize = 15; 
axisSize = 10; 
figure();
REMOVE = 2;
plot([0:L-REMOVE-1],F(1:L-REMOVE),'-o','LineWidth',lineSize)
set(gca, 'FontSize',axisSize);
set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('Fraction', 'FontSize',legSize);




%}



%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters fraction-25
%}


