function [condition] = fit(obj, rbf, polydegree, rho, regFn)
%fit Fits an interpolant to the point cloud of an surfFromPC object.
%	[obj] = fit(obj, rbf, polydegree, rho, regFn)
%	For a global smoothing parameter, use a scalar for rho.
%	For GCV estimated rho, use a range [minrho,maxrho] for rho.

global waitbarToggle

if isempty(obj.patches) || isempty(obj.radii)
	error("This object doesn't have an octree yet.")
end

if ~obj.hasOrientedNormals
	warning("Normals have not been oriented.")
end

if ~obj.hasOffSurfPts
	warning("Point cloud has not been augmented with off-surface points.")
end

obj.pcChangedSinceFit = false;

numPatches = size(obj.patches,1);
obj.weights = cell( 1, numPatches );
obj.centres = cell( 1, numPatches );
condition = NaN( 1, numPatches );

obj.fitOpts = struct("rbf", rbf, "polydegree", polydegree, ...
	"regFn", regFn);

if length(rho) == 1
	obj.fitOpts.rho = rho;
elseif length(rho) == 2
	obj.fitOpts.rho = NaN( 1, numPatches );
else
	error("Invalid rho.")
end

if exist('waitbarToggle','var') && waitbarToggle
	wb = waitbar(0,sprintf('%d out of %d patches done.',0,numPatches),...
		'Name','Fitting RBFPUM Interpolant...');
	wbCleanup = onCleanup(@() delete(wb));
end

for i = 1:numPatches
	
	J = findNeighborsInRadius( obj.pcAug, obj.patches(i,:), obj.radii(i) );
	
	obj.centres{i} = obj.pcAug.Location(J,:);
	
	if length(rho) == 1
		[obj.weights{i}, condition(i)] = ...
			surfFromPC.fitLocal( obj.centres{i}, obj.f(J), ...
			rbf, regFn, rho, polydegree );
	elseif length(rho) == 2
		[obj.weights{i}, obj.fitOpts.rho(i), condition(i)] = ...
			surfFromPC.fitLocalGCV( obj.centres{i}, obj.f(J), ...
			rbf, regFn, rho, polydegree );
	end
	
	if exist('waitbarToggle','var') && waitbarToggle
		waitbar( i/numPatches, wb, ...
			sprintf('%d out of %d patches done.',i,numPatches) )
	end
	
end

obj.Ni = cellfun(@(X) size(X,1), obj.centres);

if exist('waitbarToggle','var') && waitbarToggle
	close(wb)
end

end