function [counts] = buildOctree(obj, Nmin, Nmax, expand)

arguments
	obj surfFromPC
	Nmin (1,1)
	Nmax (1,1)
	expand (1,1) = 1.01
end

obj.octreeOpts = struct("Nmin",Nmin,"Nmax",Nmax,"expand",expand);

ranges = range(obj.pcAug.Location);
mins = min(obj.pcAug.Location);

patches = mins + ranges/2; % right in the middle
side_lengths = max(ranges); % must envelop all the points
radii = sqrt(3)*side_lengths/2;

counts = obj.pcAug.Count;

while max( counts ) > Nmax
	
	% Find most popular cube
	[ ~, i ] = max( counts );
	
	% Cartesian product def'n of new cube centres
	sets = mat2cell( ...
		[patches(i,:)'-side_lengths(i)/4, ...
		patches(i,:)'+side_lengths(i)/4], ...
		[1,1,1], 2 );
	[ X,Y,Z ] = meshgrid( sets{:} );
	new_patches = [ X(:), Y(:), Z(:) ];
	
	% Choose new side length
	new_side_lengths = repmat( 0.5*side_lengths(i), 8, 1 );
	
	to_add = zeros(8,1,'logical');
	new_counts = zeros(8,1);
	new_radius = expand*sqrt(3)*side_lengths(i)/4;
	for j = 1:8 % for each new cube
		ind = findNeighborsInRadius( obj.pcAug, new_patches(j,:), new_radius );
		if isempty(ind)
			% trim this sphere, it's empty
		else
			to_add(j) = true;
			new_counts(j) = length(ind);
		end
	end%for
	
	% Delete old and add new
	patches = [ patches([1:i-1,i+1:end],:); new_patches(to_add,:) ];
	side_lengths = [ side_lengths([1:i-1,i+1:end]); new_side_lengths(to_add) ];
	radii = [ radii([1:i-1,i+1:end]); repmat(new_radius,sum(to_add),1) ];
	counts = [ counts([1:i-1,i+1:end]); new_counts(to_add) ];
	
end%while

for i = 1:length(counts)
	if counts(i) < Nmin
		[~,dists] = findNearestNeighbors(obj.pcAug,patches(i,:),Nmin);
		radii(i) = max(dists);
		counts(i) = Nmin;
	end
end

obj.patches = patches;
obj.radii = radii;

end%function