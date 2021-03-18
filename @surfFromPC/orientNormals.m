function [] = orientNormals(obj,coarseGrid,graphNbrs)

if obj.hasOrientedNormals
	warning("Overwriting existing normal orientation.")
end

% Heavy downsample
pcDS = pcdownsample(obj.pc,"gridAverage",coarseGrid);

% Calc normals
pcDS.Normal = pcnormals( pcDS, graphNbrs );

% Preallocate sparse adjacency matrix
G = spalloc(pcDS.Count,pcDS.Count,graphNbrs*pcDS.Count);

% Calulate graph weights by how 'parallel' normals are
for i = 1:pcDS.Count
	
	[nbrs,dists] = ...
		pcDS.findNearestNeighbors(pcDS.Location(i,:),graphNbrs);
	nbrs = nbrs(dists < 2*coarseGrid);
	weights = 1 - abs(sum( pcDS.Normal(i,:) .* pcDS.Normal(nbrs,:), 2 ));
	G(i,nbrs) = weights;
	G(nbrs,i) = weights;
	G(i,i) = 0;
	
end

Gr = graph(G);

% Find minimal spanning tree of graph
T = minspantree( Gr, "Method","sparse", "Type","forest" );

% Search (breadth-first) along graph
search = bfsearch(T,1,'edgetonew','Restart',true);

% Align normals of points that are connected in the minimal spanning tree
for edge = 1:size(search,1)
	
	i = search(edge,1);
	j = search(edge,2);
	if dot( pcDS.Normal(i,:), pcDS.Normal(j,:) ) < 0
		pcDS.Normal(j,:) = -pcDS.Normal(j,:);
	end
	
end

nbrs = knnsearch(pcDS.Location,obj.pc.Location);
dots = sum( pcDS.Normal(nbrs,:) .* obj.pc.Normal, 2 );
obj.pc.Normal(dots<0,:) = -obj.pc.Normal(dots<0,:);

obj.normalEstOpts.coarseGrid = coarseGrid;
obj.normalEstOpts.graphNbrs = graphNbrs;

obj.hasOrientedNormals = true;
obj.pcChangedSinceFit = true;

end