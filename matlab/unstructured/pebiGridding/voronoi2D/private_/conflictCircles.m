function [I] = conflictCircles(Pts, CC, CR)
    TOL = 100*eps;

    nc = size(CC,1);
    np = size(Pts,1);

    %CRSqr = CR.^2;
    I = cell(numel(CR),1);
    for i = 1:nc
        dist = sqrt(sum((repmat(CC(i,:),np,1)-Pts).^2,2));
        I{i} = find(dist<CR(i)-CR(i)*TOL);
    end
end
