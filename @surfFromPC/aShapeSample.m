function [S, ashape] = aShapeSample(OBJ, alpha, Hmax)
% Gets an isosurface with marching tetrahedra, after meshing an alpha
% shape of the augmented point cloud.
% Be warned: MATLAB sometimes fails to mesh the alpha shape.

ashape = alphaShape( OBJ.pcAug.Location, alpha );
S = struct("faces",[],"vertices",[]);

for regionID = 1:ashape.numRegions
	
	[elements, nodes] = ashape.alphaTriangulation(regionID);
	[elements, nodes] = OBJ.fixTriPinch(elements,nodes);
	
	tris = freeBoundary( triangulation( elements, nodes ) );
	
	model = createpde();
	geometryFromMesh(model, nodes', tris');
	MSH = generateMesh(model, "GeometricOrder", "Linear", "Hmax", Hmax );
	
	nodes = MSH.Nodes';
	elements = MSH.Elements';
	
	F = OBJ.eval(nodes);
	
	[newfaces, newvertices] = ...
		marchingTetra( nodes, elements, F );
	
	S.faces = [S.faces; newfaces + size(S.vertices,1)];
	S.vertices = [S.vertices; newvertices];
	
end