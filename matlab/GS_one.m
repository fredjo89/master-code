clc; clear all; close all;

    

    

    Temp1 = 8;
    Temp2 = 4;
    iterations = 40 ;
    tol = -1;

    [nx ny] = deal(Temp1);
    theBasisDesider = ceil(Temp1/Temp2);        % Desides how how to devide up the coarse grid.

    %% Create Grid
    G = cartGrid([nx ny], [2, 2]);
    G = computeGeometry(G);

    %% Create Coarse Grid
    cgxy =ceil(nx/theBasisDesider);
    pv = partitionUI(G, [cgxy, cgxy]);
    CG = generateCoarseGrid(G,pv);
    CG = coarsenGeometry(CG); 
    CG = storeInteractionRegionCart(CG);

    rock = makeRock(G, 1, 1);

    gravity reset off; 
    fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
    initState = initResSol(G,0.0, 1.0); 

     %% Create A
    hT = computeTrans(G,rock);
    hT = hT* 1/fluid.properties(initState);
    A = getIncomp1PhMatrix(G, hT);


    %% Create sourceterms
    flowRate = 10^(-16)*sum(poreVolume(G,rock));
    
    src = addSource([], 1,flowRate);
    src = addSource(src,G.cells.num ,-flowRate);
   
    
    for i = 1:nx-1
        src = addSource(src, i*ny+1,flowRate); 
        src = addSource(src, i*ny , -flowRate);
    end
    



    %% Get the basis functions
    
    [Jbasis,it1] = getJBasis(CG,A,iterations, tol, 0.95);
    

    [Jbasis2,it2] = getOrderedGSBasis(CG,A,iterations, tol, 1);
    
    
    figure()
    plotToolbar(G, Jbasis.B-Jbasis2.B)
    figure()
    plotToolbar(G, Jbasis2.B)
    
    diff = sum(sum(abs(Jbasis.B-Jbasis2.B)))/Temp1^2
    
    
    
    
    
    
    
    
    
    
    %[Obasis,it3] = getOrderedGSBasis3(CG,A,iterations, tol, 1);
    
    %[Obasis_2,it4] = getOrderedGSBasis2(CG,A,iterations, tol, 1);
    
    
    
    %FSsolution   =  incompTPFA(initState, G, hT, fluid, 'src',src);
    %Osolution =  incompMultiscale(initState,CG,hT,fluid,Obasis,'src',src);
    
    %diff = abs(FSsolution.pressure - Osolution.pressure);
   % InfNorm = max(diff)/max(abs(FSsolution.pressure));
    %twoNorm = sqrt(sum(diff.^2)/sum(FSsolution.pressure.^2));
    
    
    %{
    plotCellData(G,full(Obasis.B(:,5) ));
    figure();
    plotCellData(G,full(GSbasis.B(:,5) ));
    outlineCoarseGrid(G,pv)
    %}
    

    %% Solving
    %{
    
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
    
    %}
    
    

