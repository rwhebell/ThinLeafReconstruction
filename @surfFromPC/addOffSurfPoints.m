function [] = addOffSurfPoints(obj,h,L)

if obj.hasOffSurfPts
	warning("Overwriting existing off-surface points.")
end

off_pts = ...
	[ obj.pc.Location(1:h:end,:) + L*obj.pc.Normal(1:h:end,:);
	obj.pc.Location(1:h:end,:) - L*obj.pc.Normal(1:h:end,:) ];

num_ext = size(off_pts,1);

obj.pcAug = pointCloud([ obj.pc.Location; off_pts ]);
obj.f = [ zeros(obj.pc.Count,1); ...
	L*ones(num_ext/2,1); ...
	-L*ones(num_ext/2,1) ];

obj.pcAug.Color = zeros( obj.pcAug.Count, 3, 'uint8' );
obj.pcAug.Color(obj.f == -L, 1) = 255;
obj.pcAug.Color(obj.f == L, 3) = 255;

obj.hasOffSurfPts = true;
obj.offSurfPtOpts = struct("h",h,"L",L);

end 