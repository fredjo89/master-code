clc; clear all; close all;

R = [1 10 20 30 40 50 60 70 80 90 100 120 140 160 180 200 250 300 350 400 500 600 700 800 900 1000 1500 2000 4000 6000 8000 10000];
R = (1:5);
R = [1 2 2^2 2^3 2^4 2^5 2^6 2^7 2^8 2^9 2^10 2^11 2^12 2^13 1000  ];

%R = 50; 
R_prc = zeros(length(R),1);
ratio = zeros(length(R),1);


R = 1000; 


for i = 1 : length(R)
    Nx = 2; 
    Ny = 2; 
    nx = round(Nx*sqrt(R(i)));
    ny = round(Ny*sqrt(R(i)));
    
    R_prc(i)  = nx*ny/(Nx*Ny);
   
    G = cartGrid([nx, ny]);
    G = computeGeometry(G);
    pv = partitionUI(G, [Nx, Ny]);
    CG = generateCoarseGrid(G,pv);
    CG = coarsenGeometry(CG); 
    CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);

    ratio(i) = boundaryRatio(CG,G);
end

for i = 1 : length(R)
    Nxyz = 2; 
    nxyz = round ( Nxyz*R(i)^(1/3) );
  
    R_prc_3D(i) = nxyz^3 / Nxyz^3;
   
    G = cartGrid([nxyz, nxyz,nxyz]);
    G = computeGeometry(G);
    pv = partitionUI(G, [Nxyz, Nxyz,Nxyz]);
    CG = generateCoarseGrid(G,pv);
    CG = coarsenGeometry(CG); 
    CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);

    ratio_3D(i) = boundaryRatio(CG,G);
end

[lol, boundary] = boundaryRatio(CG,G);

ratio(length(R))
ratio_3D(length(R))

my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;
my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;
my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;
my_red_3 = [241 36 35] ./ 255;

%FigHandle = figure('Position', [1200, 1000, 1000, 600]);
FigHandle = figure('Position', [1200, 200, 13*29, 11.5*29]);
loglog(R_prc,ratio,'--o','Color',my_blue_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_blue_1, 'MarkerFaceColor', my_blue_1); 
hold on; 
loglog(R_prc_3D,ratio_3D,'--o','Color',my_green_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1); 

LEG1 = legend('2D','3D');
set(LEG1);
xlabel('Upscaling factor');
ylabel('Ratio');
axis([0,10000,0.01,1]);
set(gca,'fontsize',15)


%{
FigHandle = figure('Position', [1200, 200, 1000, 600]);
semilogx(R_prc,ratio,'o--')
hold on; 
semilogx(R_prc_3D,ratio_3D,'--*')
%}

%{
my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;
my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;
my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;
my_red_3 = [241 36 35] ./ 255;

colorMatrix = [ 1 1 1; .55 .55 .55; my_red_3; my_blue_1; my_green_2];

plotCellData(G,boundary,'EdgeColor', 'y')
plotGrid(G,'FaceColor', 'none')
axis equal tight off; 
%caxis([0 max(pv)]);
colormap(colorMatrix)
%}

%{
FigHandle = figure('Position', [200, 1000, 1000, 600]);
plotCellData (G , boundary, 'EdgeColor','k'); view (45,45);
axis tight off , set ( gca , 'DataAspect',[1 1 1])
h= colorbar ('horiz');
shading flat, axis equal off, set(gca, 'zdir', 'reverse'), box on;
colorbar('horiz');

%}






