function [T,V] = fixTriPinch(T,V)
% Small helper function that sometimes stops MATLAB from failing
% to mesh an alpha-shape

tet = triangulation(T,V); % make use of triangulation class' functionality

F = tet.freeBoundary();

tr = triangulation(F,V);

E = tr.edges;                             % get edge list
TI = tr.edgeAttachments(E(:,1), E(:,2));  % get lists of faces attached to each edge
edgeDegrees = cellfun(@(e) numel(e), TI); % count those faces
badEdges = E(edgeDegrees > 2, :);        % edges without exactly 2 attached faces are bad

newtetralist = [];

for k = 1:size(badEdges,1)
	
	edge = badEdges(k,:);
	
	faces = tr.edgeAttachments(edge);
	faces = faces{:};
	
	% find face normals
	N = tr.faceNormal(faces');
	
	% find pairwise dot products of face normals
	Nsim = N*N' - eye(length(faces));
	
	% find most similar face normals
	[~,mostSim] = max(Nsim);
	
	pairs = { [1,mostSim(1)], setdiff(1:4,[1,mostSim(1)]) };
	
	for j = 1:length(pairs)
		pair = faces(pairs{j});
		newedge = [ 0 0 ];
		for i = 1:length(pairs)
			% find which vertices are not on the bad edge
			face = pair(i);
			if isequal(F(face,:),[79156 82252 79156])
				keyboard
			end
			newedge(i) = setdiff(F(face,:),edge);
		end
		% make two triangles between those and each side of the bad edge
		newtetralist = [ newtetralist; [newedge, edge] ]; %#ok<AGROW>
	end
	
end

T = [ T; newtetralist ];

end

