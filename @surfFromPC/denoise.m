function [inlierIndices, outlierIndices] = ...
	denoise(obj,NumNeighbors,Threshold)
%denoise Removes outliers from the point cloud of an surfFromPC object.
%	Arguments:
%	NumNeighbors - num nbrs to calculate average dist to nbrs
%	Threshold - if a point's average dist to nbrs is this many standard
%		deviations above the mean, it is an outlier.

% Warn if this is not the first time pc has been downsampled
if obj.hasBeenDenoised
	warning("Point cloud had already been denoised.")
end

% Store options
obj.denoiseOpts.NumNeighbors = NumNeighbors;
obj.denoiseOpts.Threshold = Threshold;

% Denoise
[obj.pc, inlierIndices, outlierIndices] = ...
	pcdenoise(obj.pc, ...
	"NumNeighbors", NumNeighbors, ...
	"Threshold", Threshold );

obj.hasBeenDenoised = true;
obj.pcChangedSinceFit = true;

end