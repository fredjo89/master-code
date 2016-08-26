function [F] = createFaultGridPoints(faultLines,faultGridSize, varargin) 
% Places fault grid points on both sides of faults.
%
% SYNOPSIS:
%   F = createFaultGridPoints(F,faultGridSize)
%   F = createFaultGridPoints(..., 'Name1', Value1, 'Name2', Value2, ...)
%
% Parameters:
%   faultLine       A cell of arrays. Each nx2 array in the cell contains 
%                   the piecewise linear approximation of a fault. The 
%                   values must be sorted along the line, e.g., a fault 
%                   consisting of two lines would be [x1,y1; x2,y2; x3,y3].
%                      .----------.-------------.
%                   (x1,y1)    (x2,y2)       (x3,y3)
%
%   faultGridSize   Desired distance between fault points along a fault.
%                   The true distance is set such that it is a factor of 
%                   the total fault length, and therefore might be slightly
%                   smaller than the desired distance.
%   
%   circleFactor    - OPTIONAL.
%                   Default value 0.6. The ratio between faultGridSize and 
%                   the radius of the circles used to create the points. If
%                   circleFactor is increased, the distance between the 
%                   fault and the points are increased. Valid values for
%                   circleFactor are in the interval (0.5,1.0).
%
%   distFun         - OPTIONAL.
%                   Default value @(x) faultGridSize*constFunc(x). Function
%                   handle that specify a relative fault grid size. If
%                   distFun = @(x) 1-0.5*x(:,1), in the unit square, then
%                   any any faults on the right side will have about twice
%                   as many grid points as equivalent faults placed on the
%                   left side. 
%
%   fCut            - OPTIONAL.
%                   Default value array of zeros. Array of length equal the
%                   number of faults. The value equals the output of the 
%                   function [~, fCut,~]=splitFaults. The value of element
%                   i tells if fault i share a start or end point with 
%                   other faults. If the value is 1 it share an end point.
%                   If the value is 2 it share a start point. If the value
%                   is 3 it share both a start point and an end point.
%
%   fwCut           - OPTIONAL.
%                   Default value array of zeros. Array of length equal the
%                   number of faults. The value equals the output of the 
%                   function [~,~, fwCut] = splitFaults. The value of 
%                   element i tells if the start or end point of fault i 
%                   should start half a step length from the start/end. If
%                   the value is 1 it starts half a step length from the 
%                   end. If the value is 2 it starts half a step length 
%                   from the start. If the value is 3 it starts half a step
%                   length from both the start and end.
%
% RETURNS:
%   F               Struct with elements:
%     F.f.pts       Point coordinates
%     F.f.Gs        Grid spacing for each fault point. This is the distance
%                   between the two points on oposite sides of a fault.
%     F.f.c         Map from fault points to circles
%     F.f.cPos      Fault point i was created using circle
%                   F.f.c(F.f.cPos(i):F.f.cPos(i+1)-1).
%     F.c.CC        Coordinates to circle centers
%     F.c.R         Radius of circles
%     F.c.f         Map from circles to fault points
%     F.c.fPos      Circle i created fault points
%                   F.c.f(F.c.fPos(i):F.c.fPos(i+1)-1)
%     F.l.fPos      Fault points F.l.fPos(i):F.l.fPos(i+1)-1 is created
%                   for fault i.
%     F.l.l         faultLine from input.
%
% EXAMPLE
%   fl = {[0.2,0.2;0.8,0.8]};
%   gS = 0.1;
%   F = createFaultGridPoints(fl,gS);
%   figure(); hold on
%   plot(fl{1}(:,1), fl{1}(:,2))
%   plot(F.f.pts(:,1),F.f.pts(:,2),'.r','markersize',20)
%
% SEE ALSO:
%   pebiGrid, compositePebiGrid, createWellGridPoints, splitFaults, pebi.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

% Load options

fh  = @(x) faultGridSize*constFunc(x);
opt = struct('distFun',     fh,  ...
             'circleFactor', 0.6, ...
             'fCut',         zeros(numel(faultLines),1), ...
             'fwCut',        zeros(numel(faultLines),1));
opt = merge_options(opt, varargin{:});

fh           = opt.distFun;
circleFactor = opt.circleFactor;
fCut         = opt.fCut;
fwCut        = opt.fwCut;


% Initialize variables.
F.f.Gs    = [];                % Fault point grid size
F.f.pts   = [];                % Fault points
F.f.c     = [];                % Map from a fault to the circle center
F.f.cPos  = 1;                 

F.c.CC    = [];                % Center of circle used to create fault pts
F.c.R     = [];                % Radius of the circle
F.c.f     = [];                % Map from the circle center to a fault  
F.c.fPos  = 1;
F.l.f     = [];
F.l.fPos  = 1;                  % Map frm fault lines to fault points
F.l.l     = faultLines;
F.l.nFault = numel(faultLines);



for i = 1:F.l.nFault
  faultLine     = F.l.l{i};
  sePtn         = .5*[fwCut(i)==2|fwCut(i)==3; fwCut(i)==1|fwCut(i)==3];
  [p, fracSpace, fCi, fRi, f2ci,cPos, c2fi,fPos] =   ...
                          faultLinePts(faultLine,     ... 
                                       faultGridSize,...
                                       circleFactor, ...
                                       fCut(i),sePtn,fh);
  nl = size(p,1)/2;
  if nl==0 % No fault points created
    F.l.fPos = [F.l.fPos; F.l.fPos(end)];
    continue
  end
  F.f.Gs   = [F.f.Gs;fracSpace];
  F.l.fPos = [F.l.fPos; size(F.f.pts,1)+1+size(p,1)];
  F.f.pts  = [F.f.pts;p];
  F.c.CC   = [F.c.CC; fCi];
  F.c.R    = [F.c.R; fRi]; 
  F.f.c    = [F.f.c; f2ci + size(F.c.fPos,1)-1];
  F.f.cPos = [F.f.cPos; cPos(2:end) + F.f.cPos(end)-1];
  F.c.f    = [F.c.f; c2fi + size(F.f.pts,1)-nl*2];
  F.c.fPos = [F.c.fPos; fPos(2:end) + F.c.fPos(end)-1];
end
F.l.f = (1:F.l.fPos(end)-1)';

% Add well-fault intersections
if ~isempty(F.c.CC)
  cutAtEnd = fwCut==1|fwCut==3;
  cutAtStr = fwCut==2|fwCut==3;
    
  endCirc = F.f.c(F.f.cPos(F.l.fPos([false;cutAtEnd]))-1); %F.l.f is has not changed
  strCirc = F.f.c(F.f.cPos(F.l.fPos(cutAtStr)));

  p       = circCircInt(F.c.CC(strCirc,:), F.c.R(strCirc),...
                        F.c.CC(endCirc,:), F.c.R(endCirc));
  l2fId1  = repmat(F.l.fPos(cutAtStr),1,2);
  l2fId2  = repmat(F.l.fPos([false;cutAtEnd]),1,2); %Not -1 because of how insertVec works
  l2fId   = reshape([l2fId1,l2fId2]',[],1);
  fId     = (size(F.f.pts,1)+1:size(F.f.pts,1) + size(p,1))';
  fId     = reshape(fId,2,[]);
  fId     = repmat(fId,2,1);
  cId     = reshape([strCirc,endCirc]',[],1);
  c2fId   = repmat(F.c.fPos(cId)',2,1);
  c2fId   = c2fId(:);
  
  
  F.l.f   = insertVec(F.l.f, fId(:), l2fId);
  F.f.pts = [F.f.pts;p];
  nGs     = repmat(sqrt(sum(diff(p).^2,2)),1,2)';
  F.f.Gs  = [F.f.Gs;reshape(nGs(:,1:2:end),[],1)];
  F.c.fPos= F.c.fPos + ...
    cumsum(accumarray([cId+1;size(F.c.fPos,1)],2*[ones(1,size(cId,1)),0]));
  F.c.f   = insertVec(F.c.f, fId(:), c2fId);
  cId     = repmat(reshape(cId,2,[]),2,1);
  F.f.cPos= [F.f.cPos;F.f.cPos(end)+2*cumsum(ones(size(p,1),1))];
  F.f.c   = [F.f.c;cId(:)];
  
  F.l.fPos(2:end) = F.l.fPos(2:end)+ 2*cumsum(cutAtEnd+cutAtStr);
end


% Merge fault intersections
if ~isempty(F.f.pts)
  % Remove duplicate fault centers
  [~, IA, IC] = uniquetol(F.c.CC,'byRows',true);
  F.c.CC      = F.c.CC(IA,:);
  F.c.R       = F.c.R(IA);
  [~,I]       = sort(IC);
  map         = [F.c.fPos(1:end-1), F.c.fPos(2:end)-1];
  map         = map(I,:);
  map         = arrayfun(@colon, map(:,1),map(:,2),'uniformOutput',false);
  F.c.f       = F.c.f(cell2mat(map'));
  fNum        = diff(F.c.fPos);
  F.c.fPos    = cumsum([1; accumarray(IC,fNum)]);
  F.f.c       = IC(F.f.c);
  % Merge intersections
  [F]         = fixIntersections(F);
end
end

function [Pts, gridSpacing, circCenter, circRadius, f2c,f2cPos, c2f,c2fPos] = ...
    faultLinePts(faultLine, fracDs, circleFactor, isCut,sePtn, fh) 

    assert(0.5<circleFactor && circleFactor<1)
    assert(size(faultLine,1)>1 && size(faultLine,2)==2);
    
    % interpolate fault line to get desired grid spacing. 
    circCenter = interFaultLine(faultLine, fh, fracDs,sePtn);
    numOfFracPts = size(circCenter,1)-1;
    
    % Test faultLine
    if numOfFracPts == 1
      d = sqrt(sum((circCenter(2,:)-circCenter(1,:)).^2, 2));
      if d < 0.8*fh((circCenter(2,:)+circCenter(1,:))/2);
        Pts         = [];
        gridSpacing = [];
        circCenter  = [];
        circRadius  = [];
        f2c         = [];
        f2cPos      = [];
        c2f         = [];
        c2fPos      = [];
        return
      end
    end
    % Calculate the line lenth and circle radiuses. If you experience
    % imaginary faultOffset you might want to try the max lineLength
    % instead of the mean.
    lineLength = sqrt(sum((circCenter(2:end,:)-circCenter(1:end-1,:)).^2, 2));
    circRadius = circleFactor*[lineLength(1); ...
                              (lineLength(1:end-1) + lineLength(2:end))/2; ...
                               lineLength(end)];
                             
    switch isCut
      case 1
        circRadius(end) = fh(circCenter(end,:))*circleFactor;
      case 2
        circRadius(1)   = fh(circCenter(1,:))*circleFactor;
      case 3
        circRadius(1)   = fh(circCenter(1,:))*circleFactor;
        circRadius(end) = fh(circCenter(end,:))*circleFactor;
    end
    
    % Calculate the crossing of the circles
    bisectPnt = (lineLength.^2 - circRadius(2:end).^2 + circRadius(1:end-1).^2)...
                ./(2*lineLength);
    faultOffset = sqrt(circRadius(1:end-1).^2 - bisectPnt.^2);
    n1 = (circCenter(2:end,:)-circCenter(1:end-1,:))./repmat(lineLength,1,2); %Unit vector
    n2 = [-n1(:, 2), n1(:,1)];                                                %Unit normal
    
    % Set fault points on left and right side of fault
    left   = circCenter(1:end-1,:) + bsxfun(@times, bisectPnt, n1)  ...
             + bsxfun(@times, faultOffset, n2);
    right  = circCenter(1:end-1,:) + bsxfun(@times, bisectPnt, n1)  ...
             - bsxfun(@times, faultOffset, n2);
         
    % Put together result
    Pts = [right;left];
    f2c = cumsum(accumarray((1:2:size(Pts,1)+1)',1));
    f2c = repmat(f2c(2:end),2,1); 
    f2cPos = (1:2:numel(f2c)+1)';
    nf  = size(left,1);
    c2f = [   nan,          nan,         1,       nf+1;...
           (1:nf-1)', (nf+1:2*nf-1)', (2:nf)', (nf+2:2*nf)';...
              nf,           2*nf,       nan,     nan]';
    c2f = c2f(3:end-2)';
    c2fPos = [1;(3:4:numel(c2f))';numel(c2f)+1];
    gridSpacing = 2*[faultOffset;faultOffset];
end

function [p] = interFaultLine(line, fh, lineDist,sePtn, varargin)
    % Interpolate a fault line. 
    % Arguments:
    %   line        Coordinates of the fault line. Must be ordered.
    %   fh          A function handle for the relative distance function 
    %               for the interpolation fh = 1 will give a equiv distant
    %               interpolation
    %   lineDist    Scalar which set the distance between interpolation
    %               points (Relative to fh = 1)
    % varargin:    
    %               Arguments passed to fh

    % Parameters
    TOL = 1e-4; maxIt = 10000;

    % Create initial points, equally distributed.
    p = eqInterpret(line, lineDist,sePtn);
    % add auxillary points
    if sePtn(1)~=0, p = [line(1,:);p];   end
    if sePtn(2)~=0, p = [p;line(end,:)]; end
    count=0;
    while count<maxIt
      count = count+1;
      % Calculate distances, and wanted distances
      d = distAlLine(line, p);
      pmid = (p(1:end-1,:) + p(2:end,:))/2;
      dw = fh(pmid,varargin{:});
      if sePtn(1)~=0, dw(1)   = dw(1).*sePtn(1);   end
      if sePtn(2)~=0, dw(end) = dw(end).*sePtn(2); end

      % Possible insert or remove points
      if sum(d - dw)>min(dw)
          [~, id] = max(d - dw);
          p = [p(1:id,:); pmid(id,:); p(id+1:end,:)];
          continue
      elseif sum(d - dw)<-max(dw)
          [~, id] = min(d - dw);
          if id == 1, id = 2; end
          p = p([1:id-1,id+1:end],:);
          continue
      end
      % If we only have external nodes, we can do nothing.
      if size(p,1)<=2, return, end
      % Move points based on desired length
      Fb = dw - d;                       % Bar forces
      Fn = Fb(1:end-1) - Fb(2:end);      % Force on internal nodes
      moveNode = Fn*0.2;                 % Movement of each internal node.
      d = d + [moveNode(1); moveNode(2:end) - moveNode(1:end-1); -moveNode(end)];
      p = interpLine(line,d);            % Update node positions

      % Terminate if Nodes have moved (relative) less  than TOL
      if all(abs(moveNode)<TOL*lineDist), break; end
    end
    
    if sePtn(1)~=0, p = p(2:end,:);end
    if sePtn(2)~=0, p = p(1:end-1,:);end
    if count == maxIt
        warning('Fault interpolation did not converge.')
    end

end


function [d] = distAlLine(line, p)
    % Calculates the distace between consecutive interpolation points along
    % line
    % Arguments:
    %   line    line that is interpolated
    %   p       Interpolation points
    % Returns:
    %   d       distance between consecutive points of p, along line

    TOL = 50*eps;
    
    N = size(p,1);
    d = zeros(N-1,1);
    jointDist = 0;
    for i = 1:size(line,1)-1
        lineStart = repmat(line(i,:), size(p,1),1);
        lineEnd = repmat(line(i+1,:), size(p,1),1);
        distA = eucDist(lineStart, p) + eucDist(p,lineEnd);
        distB = eucDist(lineStart,lineEnd);
        indx  = find(abs(distA - distB) < TOL*distB);
        if numel(indx)==0 
            jointDist = jointDist + eucDist(line(i,:), line(i+1,:));
            continue
        elseif numel(indx)>=2
            d(indx(1:end-1)) = sqrt(sum((p(indx(1:end-1),:) ... 
                             - p(indx(2:end),:)).^2,2));
            
        end
        if indx(1)>1 && eucDist(line(i,:),p(indx(1),:))>TOL
            d(indx(1)-1) = jointDist + eucDist(line(i,:), p(indx(1),:));
        end
        jointDist = eucDist(p(indx(end),:), line(i+1,:));
    end
end


function [d] = eucDist(a, b)
    d = sqrt(sum((a - b).^2,2));
end


function [newPoints] = interpLine(path, dt)
    distS = sqrt(sum(diff(path,[],1).^2,2));
    t = [0; cumsum(distS)];
    
    newPtsEval = [0; cumsum(dt)];
    newPtsEval(end) = t(end); % Last point can not move

    newX = interp1(t,path(:,1),newPtsEval);
    newY = interp1(t,path(:,2),newPtsEval);
    newPoints = [newX,newY];
end


function [F] = fixIntersections(F)
  assert(all(diff(F.f.cPos)==2),'all points must be created from exactly 2 circles');
  
  % Find conflict circles
  I = conflictCircles(F.f.pts, F.c.CC, F.c.R);
  
  circ    = find(~cellfun(@isempty,I));
  if isempty(circ)
    return
  end

  circNum = cellfun(@numel, I(circ));
  id      = zeros(sum(circNum),1);
  circPos = cumsum([1; circNum]);
  id(circPos(1:end-1)) = 1;
  id      = cumsum(id);
  circ    = [circ(id), F.f.c([F.f.cPos(vertcat(I{circ})),...
             F.f.cPos(vertcat(I{circ}))+1])];
  
  % Find shared circle
  [neigh,neighPos] = findNeighbors(circ(:,1),F.c.f,F.c.fPos, F.f.c,F.f.cPos);
  assert(all(diff(neighPos)==2));
  
  neigh  = reshape(neigh,2,[])';
  shared = 2*any(bsxfun(@eq, neigh, circ(:,2)),2) ...
          +3*any(bsxfun(@eq, neigh, circ(:,3)),2);
  keep   = find(shared);
  circ   = circ(keep,:);
  swap   = shared(keep)==2;                % Set shared circle at third row
  circ(swap,:) = [circ(swap,1),circ(swap,3),circ(swap,2)];
  
  % Remove duplicate pairs
  [~,IA] = unique(sort(circ,2),'rows');
  circ   = circ(IA,:);
  
  % Calculate new radiuses
  line = [F.c.CC(circ(:,3),:),reshape(mean(reshape(F.c.CC(circ(:,1:2)',:),2,[]),1),[],2)];
  int  = lineCircInt(F.c.CC(circ(:,3),:),F.c.R(circ(:,3)), line);

  % set radius to smallest
  R = sqrt(sum((F.c.CC(circ(:,1:2),:)-[int;int]).^2,2));
  I = false(size(circ(:,1:2)));
  for i = 1:numel(R)
    if R(i)<F.c.R(circ(i))
      F.c.R(circ(i)) = R(i);
      I(circ(:,1:2)==circ(i)) = false;
      I(i) = true;
    elseif R(i)==F.c.R(circ(i))
      I(i) = true;
    end
  end
  c = unique(circ(:,1:2));
  if size(c,2)>1
    c = c';
  end
  % Calculate new Pts
  map = arrayfun(@colon,F.c.fPos(c),F.c.fPos(c+1)-1,'uniformOutput',false)';
  fId = F.c.f(horzcat(map{:})');
%   fId = unique(fId(:));
%   fId = fId(~isnan(fId));

  [neigh,neighPos] = findNeighbors(c, F.c.f,F.c.fPos, F.f.c,F.f.cPos);
  assert(all(diff(neighPos)==2));
  neigh = reshape(neigh,2,[])';
  
  p = circCircInt(F.c.CC(c,:), F.c.R(c),...
                 reshape(F.c.CC(neigh',:)',4,[])',reshape(F.c.R(neigh),[],2));
  assert(isreal(p),'Failed to merge fault crossings. Possible too sharp intersections');
  F.f.pts(fId,:) = p;
  nGs            = repmat(sqrt(sum(diff(p).^2,2)),1,2)';
  F.f.Gs(fId)    = reshape(nGs(:,1:2:end),[],1);
  map            = [F.f.cPos(fId),F.f.cPos(fId)+1]';
  F.f.c(map(:))  = [c';neigh(:,1)';c';neigh(:,1)';c';neigh(:,2)';c';neigh(:,2)'];%reshape(repmat([c',c';neigh(:,1)',neigh(:,2)'],2,1),2,[]);

  [~, IA, IC] = uniquetol(F.f.pts,'byRows',true);
  F.f.pts = F.f.pts(IA,:);
  F.f.Gs = F.f.Gs(IA);
  [~,I] = sort(IC);
  
  mapc = [F.f.cPos(1:end-1), F.f.cPos(2:end)-1];
  mapc = mapc(I,:);
  mapc = arrayfun(@colon, mapc(:,1),mapc(:,2),'uniformOutput',false);
  F.f.c = F.f.c(cell2mat(mapc'));
  cNum = diff(F.f.cPos);
  F.f.cPos = cumsum([1;accumarray(IC,cNum)]);
  F.c.f = IC(F.c.f);
  
 
  F.l.f = IC(F.l.f);
end


function [neigh,neighPos] = findNeighbors(c, c2f,c2fPos, f2c,f2cPos)
map   = arrayfun(@colon, c2fPos(c),c2fPos(c+1)-1,'uniformoutput',false);
pId   = cellfun(@(c) c2f(c), map,'uniformOutput',false);
neighMap = cellfun(@(c) cell2mat(arrayfun(@colon, f2cPos(c),f2cPos(c+1)-1,'uniformOutput',false)')...
                    ,pId,'uniformOutput',false);
neigh = cellfun(@(c) f2c(c),neighMap,'uniformOutput',false);
neigh = cellfun(@unique, neigh,'uniformOutput',false);
neigh = arrayfun(@(i) neigh{i}(neigh{i}~=c(i)),1:numel(neigh),'uniformOutput',false)';
neighPos = cumsum([1;cellfun(@numel, neigh)]);
neigh = vertcat(neigh{:});
end

function [p] = circCircInt(CC1, CR1, CC2,CR2)
if isempty(CC1) || isempty(CC2)
  p = [];
  return
end
% Expand matrices for computation
CC1 = repmat(CC1, 1,size(CC2,2)/size(CC1,2));
CR1 = repmat(CR1, 1,size(CR2,2)/size(CR1,2));
CC1 = reshape(CC1',2,[])';
CC2 = reshape(CC2',2,[])';
CR1 = reshape(CR1',1,[])';
CR2 = reshape(CR2',1,[])';

d = sqrt(sum((CC1 - CC2).^2,2));              % Distance between centers
bisectPnt = (d.^2 - CR2.^2 + CR1.^2)./(2*d);  % Mid-Point
faultOffset = sqrt(CR1.^2 - bisectPnt.^2);    % Pythagoras
n1 = (CC2-CC1)./repmat(d,1,2);                % Unit vector
n2 = [-n1(:, 2), n1(:,1)];                    % Unit normal

% Set right left and right intersection points
left   = CC1 + bsxfun(@times, bisectPnt, n1)  ...
         + bsxfun(@times, faultOffset, n2);
right  = CC1 + bsxfun(@times, bisectPnt, n1)  ...
         - bsxfun(@times, faultOffset, n2);

% Put result together
p = reshape([right,left]',2,[])';

end

function [p] = lineCircInt(CC, CR, line)
vec = line(:,3:4) - line(:,1:2);
c2l = line(:,1:2) - CC;
a   = dot(vec,vec,2);
b   = 2*dot(c2l,vec,2);
c   = dot(c2l,c2l,2) - dot(CR,CR,2);

dist    = (b.*b - 4*a.*c);
lineHit = dist>=0;
distSqr = sqrt(dist(lineHit));

%t(:,1) = -b - distSqr./(2*a); % This is the intersection on wrong side.
t = -b + distSqr./(2*a);
p = bsxfun(@times,vec,t) + line(:,1:2);


end


