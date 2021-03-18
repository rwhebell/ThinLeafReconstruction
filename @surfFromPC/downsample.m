function [] = downsample(obj,method,param)
%downsample Downsamples the point cloud of an surfFromPC object.

% Warn if this is not the first time pc has been downsampled
if obj.hasBeenDownsampled
	warning("Point cloud had already been downsampled.")
end

% Store options
obj.downsampleOpts.method = method;
obj.downsampleOpts.param = param;

% Downsample
obj.pc = pcdownsample(obj.pc, method, param);

obj.hasBeenDownsampled = true;
obj.pcChangedSinceFit = true;

end