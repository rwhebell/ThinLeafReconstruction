function S = marchingTetra( nodes, elements, F )
% nodes: Nx3 list of points
% elements: Mx4 list of tetrahedra
% F: Nx1 vector of function evals at nodes

% Any elements where all vertices have the same sign needn't be processed
elemsSign = int8(sign(F(elements)));
elemsSign(elemsSign == 0) = 1; % zero is positive!
withoutSurface = all(elemsSign == -1, 2) | all(elemsSign == 1, 2);
elements = elements( ~withoutSurface, : );

% Can't have any NaNs in there either
elemsWithNanEvals = any( isnan(F(elements)), 2 );
elements = elements( ~elemsWithNanEvals, : );

TR = triangulation( elements, nodes );
edges = TR.edges;

% Find which edges between nodes have a sign change
% (we'll need to interpolate those)
edgesSign = sign(F(edges));
edgesSign(edgesSign == 0) = 1; % zero is positive!
signchange = edgesSign(:,1) ~= edgesSign(:,2);

E = edges(signchange,:);
FE = F(E);

U = ( abs(FE(:,2)).*nodes(E(:,1),:) + abs(FE(:,1)).*nodes(E(:,2),:) ) ./ ...
	( abs(FE(:,2)-FE(:,1)) );

elements = sort(elements,2); % number nodes in ascending order
E = sort(E,2); % keep the ascending order consistent, so its easier to search

nodeCombinations = nchoosek(1:4,2); % four nodes, two to each edge

faces = zeros(2*size(elements,1),3);
numfaces = 0;

tic
for el = 1:size(elements,1)
	
	% What nodes are in this element?
	tet = elements(el, :);
	
	% What edges are in this element?
	tetEdges = tet( nodeCombinations );
	
	% Get the indices of these edges in the big edge list, E
	ind = zeros(1,4);
	ptr = 1;
	for i = 1:6
		new_ind = find( all(E==tetEdges(i,:),2) );
		ind(ptr:ptr+length(new_ind)-1) = new_ind;
		ptr = ptr + length(new_ind);
	end
	
	ind = ind(ind~=0);
	
	newfaces = trisFromSurfaceIntersections( E(ind,:), ind );
	faces( numfaces+1:numfaces+size(newfaces,1), : ) = newfaces;
	numfaces = numfaces + size(newfaces,1);
	
	if mod(el, floor(size(elements,1)/10)) == 0
		estTimeRemaining = toc/el * size(elements,1) - toc;
		fprintf('Est time remaining: %.3g mins\n', estTimeRemaining/60)
	end
	
end

S = struct( 'vertices', U, 'faces', faces(1:numfaces,:) );

end

function [curve,I] = edgeListToCurve(edges)

curve = zeros( size(edges,1)+1, 1 );
curve(1) = edges(1,1);

if nargout > 1
	I = [];
end

edges_used = zeros( size(edges,1), 1, 'logical' );
for i = 2:size(edges,1)+1
    [next_edge,j] = find( ~edges_used & (edges == curve(i-1)), 1 );
    if j == 1
        curve(i) = edges( next_edge, 2 );
    elseif j == 2
        curve(i) = edges( next_edge, 1 );
    end
    edges_used(next_edge) = true;
	
	if nargout > 1
		I = [ I, next_edge ];
	end
end

curve = curve(curve~=0);

end

function F = trisFromSurfaceIntersections( E, U )
% E: edges from the tetrahedron, where there are surface intersections
% U: indices of intersection points in their larger set

assert( size(E,2) == 2 )
assert( ismember( length(U), [3,4] ) )
assert( size(E,1) == length(U) )

if size(E,1) == 3
	
	% easy!
	F = U;
	
elseif size(E,1) == 4
	
	% not so easy :(
	
	% need to establish a loop along the edges in E
	[~, edgeOrder] = edgeListToCurve(E);
	
	% Now we have an order for U that defines a square!
	% Making two triangles is simple now
	U = U(edgeOrder);
	
	F = U([1,2,3;3,4,1]);
	
end

end