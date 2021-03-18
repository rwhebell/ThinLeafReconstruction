function [X, Y, Z, INTERP] = gridSample(obj,Ngrid,alpha)

if isempty(obj.weights)
	error("This surface does not have an indicator function to sample.")
end

if obj.pcChangedSinceFit
	warning("Point cloud has changed since this surface was fit.")
end

ranges = range(obj.pcAug.Location);
if length(Ngrid) == 1
	Ngrid = Ngrid * ranges./max(ranges);
end

ranges = num2cell(ranges);
mins = num2cell(min(obj.pcAug.Location));
maxs = num2cell(max(obj.pcAug.Location));
Ngrid = num2cell(Ngrid);

gridVecs = ...
	cellfun( @(r,x1,x2,N) single(linspace( x1-r/200, x2+r/200, N )), ...
	ranges, mins, maxs, Ngrid, 'UniformOutput', false );

[X,Y,Z] = meshgrid(gridVecs{:});
PTS = [X(:),Y(:),Z(:)];

obj.ashape = alphaShape(double(obj.pcAug.Location), alpha );
inshape = obj.ashape.inShape(double(PTS));

INTERP = NaN(size(PTS,1),1,'single');
INTERP(inshape) = obj.eval( PTS(inshape,:) );

INTERP = reshape(INTERP,size(X));

end