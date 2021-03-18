function [] = estimateNormals(obj,pcaNbrs)
%estimateNormals approximates surface normals from the point cloud of a
%	surfFromPC object.
%
%	Arguments:
%	pcaNbrs - number of neighbours to use in PCA

if ~isempty(obj.pc.Normal)
	warning("Overwriting existing normals.")
end

obj.normalEstOpts.pcaNbrs = pcaNbrs;
obj.pc.Normal = pcnormals(obj.pc, pcaNbrs);

obj.pcChangedSinceFit = true;

end