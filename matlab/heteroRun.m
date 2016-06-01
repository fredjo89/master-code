function [output] = heteroRun(Temp1,Temp2, iterations, tol)



    %% Setting up the model



    [nx ny] = deal(Temp1);
    theBasisDesider = ceil(Temp1/Temp2);        % Desides how how to devide up the coarse grid.

    %% Create Grid
    lengthX = 500; 
    lengthY = lengthX;
    G = cartGrid([nx, ny], [lengthX, lengthY]);
    G = computeGeometry(G);

    %% Create Coarse Grid
    cgxy =ceil(nx/theBasisDesider);
    pv = partitionUI(G, [cgxy, cgxy]);
    CG = generateCoarseGrid(G,pv);
    CG = coarsenGeometry(CG); 
    CG = storeInteractionRegionCart(CG);

    
    
    


    midWidth = 0;
    channelWidth = 0; 


    outerPoro = [ .3 0.5];
    channelPoro = [ 0.6 0.65];
    innerPoro = [0.1 0.2];
    
    rock =  getRock( G, nx, ny,midWidth, channelWidth,...
                    outerPoro, innerPoro, channelPoro );

    gravity reset off; 
    fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
    initState = initResSol(G,0.0, 1.0); 

     %% Create A
    hT = computeTrans(G,rock);
    hT = hT* 1/fluid.properties(initState);
    A = getIncomp1PhMatrix(G, hT);


    %% Create sourceterms
    flowRate = 10^(-16)*sum(poreVolume(G,rock));
    src = addSource([], ny,flowRate);
    src = addSource(src,G.cells.num-ny+1 ,-flowRate);
    
    %{
    for i = 1:nx-1
        src = addSource(src, i*ny+1,flowRate); 
        src = addSource(src, i*ny , -flowRate);
    end
    %}



    %% Get the basis functions
  %{
    [JbasisTwo,it1] = getJBasis(CG,A,iterations, tol, .95);
    [GSbasis,it2] = getGSBasis(CG,A,iterations, tol, 1);
    [newGSbasis,it3] = getRepPureGS(CG,A,iterations, tol, 1);
   %} 
    
    [JbasisTwo,it1] = getJBasis(CG,A,iterations, tol, 1);
    [GSbasis,it2] = getRepPureGS(CG,A,iterations, tol, 1);
    %[newGSbasis,it3] = getRepPureGS(CG,A,iterations, tol, 1.5);
    [newGSbasis,it3] = getOrderedGSBasis(CG,A,iterations, tol, 1);
    

    %% Solving

    % Solve with incompTPFA
    FSsolution   =  incompTPFA(initState, G, hT, fluid, 'src',src);
    JTwoSolution =  incompMultiscale(initState,CG,hT,fluid,JbasisTwo,'src',src);
    GSsolution   =  incompMultiscale(initState,CG,hT,fluid,GSbasis,'src',src);
    newGSsolution   =  incompMultiscale(initState,CG,hT,fluid,newGSbasis,'src',src);

    %% Computing discrepancy norms

    % Compute errorNorms


    JacobiTwo_Diff = abs(FSsolution.pressure - JTwoSolution.pressure); 
    JacobiTwo_InfNorm = max(JacobiTwo_Diff)/max(abs(FSsolution.pressure));
    JacobiTwo_twoNorm = sqrt(sum(JacobiTwo_Diff.^2)/sum(FSsolution.pressure.^2));

    GS_Diff = abs(FSsolution.pressure - GSsolution.pressure); 
    GS_InfNorm = max(GS_Diff)/max(abs(FSsolution.pressure));
    GS_twoNorm = sqrt(sum(GS_Diff.^2)/sum(FSsolution.pressure.^2));
    
    newGS_Diff = abs(FSsolution.pressure - newGSsolution.pressure); 
    newGS_InfNorm = max(newGS_Diff)/max(abs(FSsolution.pressure));
    newGS_twoNorm = sqrt(sum(newGS_Diff.^2)/sum(FSsolution.pressure.^2));


    
    output = [
            it1;
            it2;
            it3;
            JacobiTwo_InfNorm;
            GS_InfNorm;
            newGS_InfNorm;
            JacobiTwo_twoNorm;
            GS_twoNorm;
            newGS_twoNorm;
        ];
    
    
    
    
    %figure();
    %plotCellData(G,full( FSsolution.pressure ));
   
    
    
    %% Plotting results 
    
    
    %{
    lWidth = 1; 
    fSize = 1; 
    axisSize = 1; 


    figure(); 
    plotCellData(G,FSsolution.pressure);
    plotGrid(G,src.cell, 'FaceColor', 'w'); 
    axis equal tight; colormap(jet(32)); 
    c = colorbar('horiz');  
    set(gca, 'FontSize',axisSize);

    
       
    figure(); 
    plotCellData(G,JTwoSolution.pressure);
    plotGrid(G,src.cell, 'FaceColor', 'w'); 
    axis equal tight; colormap(jet(32)); 
    c=colorbar('horiz'); 
    set(gca, 'FontSize',axisSize);
    
    figure(); 
    plotCellData(G,GSsolution.pressure);
    plotGrid(G,src.cell, 'FaceColor', 'w'); 
    axis equal tight; colormap(jet(32)); 
    c=colorbar('horiz'); 
    set(gca, 'FontSize',axisSize);
    
    figure(); 
    plotCellData(G,newGSsolution.pressure);
    plotGrid(G,src.cell, 'FaceColor', 'w'); 
    axis equal tight; colormap(jet(32)); 
    c=colorbar('horiz'); 
    set(gca, 'FontSize',axisSize);


    figure(); 
    plotCellData(G,rock.poro);
    outlineCoarseGrid(G,pv,'r')
    plotGrid(G,src.cell, 'FaceColor', 'w'); 
    axis equal tight; colormap(jet(1000*128)); 
    c=colorbar('horiz'); 
    set(gca, 'FontSize',axisSize);

    
    

    
    %}
    
    
    
    
    
    

end
