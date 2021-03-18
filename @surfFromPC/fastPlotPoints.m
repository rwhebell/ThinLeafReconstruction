function varargout = fastPlotPoints(ptcld,varargin)
%plotPoints plots the point cloud of an surfFromPC object.
%	Name-value pairs will be passed to scatter3.
%

gridStep = (1/200) * min( [ ...
	ptcld.XLimits(2) - ptcld.XLimits(1), ...
	ptcld.YLimits(2) - ptcld.YLimits(1), ...
	ptcld.ZLimits(2) - ptcld.ZLimits(1) ] );

pcDS = pcdownsample(ptcld, "gridAverage", gridStep);

P = pcDS.Location;

h = scatter3( P(:,1), P(:,2), P(:,3), [], P(:,3), varargin{:} );
axis equal

if nargout > 0
	varargout = {h};
end

end