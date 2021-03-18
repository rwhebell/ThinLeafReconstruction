classdef surfFromPC < handle
	%surfFromPC Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
	end % public properties
	
	properties (SetAccess = private)
		
		pc pointCloud
		pcAug pointCloud
		f (:,1)
		patches (:,3)
		radii (:,1)
		centres
		Ni
		
		weights
		
		ashape alphaShape = alphaShape([0 0 0])
		
		hasBeenDenoised (1,1) logical = false
		hasBeenDownsampled (1,1) logical = false
		hasOrientedNormals (1,1) logical = false
		hasOffSurfPts (1,1) logical = false
		pcChangedSinceFit (1,1) logical = false
		
		denoiseOpts (1,1) struct
		downsampleOpts (1,1) struct
		normalEstOpts (1,1) struct
		offSurfPtOpts (1,1) struct
		octreeOpts (1,1) struct
		fitOpts (1,1) struct
		
	end % private properties
	
	% ------------------------
	methods
		function obj = surfFromPC(pc)
			%surfFromPC Construct an instance of this class
			%   Detailed explanation goes here
			obj.pc = pc;
		end % surfFromPC
		
		[inlierIndices, outlierIndices] = ...
			denoise(obj,NumNeighbors,Threshold)
		[] = downsample(obj,method,param)
		[] = estimateNormals(obj,pcaNbrs)
		[] = orientNormals(obj,coarseGrid,graphNbrs)
		[] = addOffSurfPoints(obj,h,el)
		counts = buildOctree(obj, Nmin, Nmax, expand)
		
		condition = fit(obj, rbf, polydegree, rho, regFn)
		
		[F, orphans, unorms, curv, norms, Hdd ] = eval(obj, x, rbfd, rbfdd)
		
		[X, Y, Z, INTERP] = gridSample(obj,Ngrid,alpha)
		S = aShapeSample(OBJ, alpha, Hmax)
		
		[] = makeLog(obj,filename)
		
	end%public methods
	
	
	% ------------------------
	methods (Static = true)
		
		varargout = fastPlotPoints(ptcld,varargin)
		
	end%public static methods
	
	
	% ------------------------
	methods (Access = private)
		
	end%private methods
	
	
	% ------------------------
	methods (Access = private, Static = true)
		
		[D,Dd,Ddd] = distMatrix(x,y)
		[z,zd,zdd] = distMatrixAct(x,y,v,rbf,rbfd,rbfdd)
		[weights, condition] = fitLocal( X, f, rbf, regFn, rho, polydegree )
		[weights, rho, condition] = ...
			fitLocalGCV( X, f, rbf, regFn, logrhoRange, polydegree )
		S = marchingTetra( nodes, elements, F )
		[T,V] = fixTriPinch(T,V)
		
	end%private static methods
	
	
end

