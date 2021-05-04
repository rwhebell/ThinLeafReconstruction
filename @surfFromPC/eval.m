function [F, orphans, unorms, curv, norms, Hdd ] = ...
	eval( obj, x, rbfd, rbfdd )

global waitbarToggle

rbf = obj.fitOpts.rbf;
polydegree = obj.fitOpts.polydegree;

[N,~] = size(x);
F = zeros(N, 1, 'like', x);
wsum = zeros(N,1);

if nargout>2
	Fd = zeros(N,3);
	norms = zeros(N,3);
	wsumd = zeros(N,3);
	assert(nargin > 2);
end
if nargout>3
	Fdd = zeros(N,6);
	curv = zeros(N,1);
	wsumdd = zeros(N,6);
	assert(nargin > 3);
end

[M,~] = size(obj.patches);

tree = KDTreeSearcher( x, "Distance", "euclidean" );

[ xpowers, ypowers, zpowers ] = ...
	meshgrid( 0:polydegree, 0:polydegree, 0:polydegree );
powers = [ xpowers(:), ypowers(:), zpowers(:) ];
powers = powers( sum(powers,2) <= polydegree, : );

if exist('waitbarToggle','var') && waitbarToggle
	wb = waitbar(0,sprintf('%d out of %d patches done.',0,M),...
		'Name','Evaluating RBFPUM Interpolant...');
	wbCleanup = onCleanup(@() delete(wb));
end

for i = 1:M
	
	% Find the points within the bounding box of this subdomain
	[ J, dists ] = ...
		rangesearch( tree, obj.patches(i,:), obj.radii(i), ...
		"SortIndices", false );
	J = J{1}';
	dists = dists{1}';
	
	if isempty(J) || isempty(obj.centres{i})
		continue
	end
	
	distsVec = x(J,:) - obj.patches(i,:);
	if nargout>2
		distsd = distsVec./dists;
	end
	if nargout>2
		distsdd = [1,0,0,1,0,1]./dists-distsVec(:,[1,1,1,2,2,3]).*distsVec(:,[1,2,3,2,3,3])./dists.^3;
	end
	
	% Calculate the weights for this local interpolant
	w = phi( dists ./ obj.radii(i) );
	if nargout>2
		wd = phid( dists ./ obj.radii(i) ).*distsd./obj.radii(i);
	end
	if nargout>3
		wdd = phidd( dists ./ obj.radii(i) ).* ...
			distsd(:,[1,1,1,2,2,3]).* ...
			distsd(:,[1,2,3,2,3,3])./obj.radii(i).^2 + ...
			phid( dists ./ obj.radii(i) ).*distsdd./obj.radii(i);
	end
	
	% Keep track of the weight sum so we can divide by it later
	wsum(J) = wsum(J) + w;
	if nargout>2
		wsumd(J,:) = wsumd(J,:) + wd;
	end
	if nargout>3
		wsumdd(J,:) = wsumdd(J,:) + wdd;
	end
	
	% Distance matrix
	Nc = size(obj.centres{i},1);
	if nargout>3
		[z,zd,zdd] = obj.distMatrixAct( x(J,:), obj.centres{i}, obj.weights{i}(1:Nc), rbf, rbfd, rbfdd);
	elseif nargout>2
		[z,zd] = obj.distMatrixAct( x(J,:), obj.centres{i}, obj.weights{i}(1:Nc), rbf, rbfd);
	else
		z = obj.distMatrixAct( x(J,:), obj.centres{i}, obj.weights{i}(1:Nc), rbf );
	end
	
	% Polynomial terms
	P = ones( length(J), size(powers,1) );
	if nargout>2
		Pd = ones( length(J), size(powers,1), 3 );
	end
	
	if nargout>3
		Pdd = ones( length(J), size(powers,1), 6 );
	end
	
	for k = 1:size(powers,1)
		Ptemp = x(J,:).^powers(k,:);
		P(:,k) = prod(Ptemp,2);
		if nargout>2
			Ptempd = powers(k,:).*x(J,:).^(powers(k,:)-1);
			Ptempd(:,powers(k,:)<=0) = 0;
			Pd(:,k,1) = Ptempd(:,1).*Ptemp(:,2).*Ptemp(:,3);
			Pd(:,k,2) = Ptemp(:,1).*Ptempd(:,2).*Ptemp(:,3);
			Pd(:,k,3) = Ptemp(:,1).*Ptemp(:,2).*Ptempd(:,3);
		end
		if nargout>3
			Ptempdd = powers(k,:).*(powers(k,:)-1).*x(J,:).^(powers(k,:)-2);
			Ptempdd(:,powers(k,:)<=1) = 0;
			Pdd(:,k,1) = Ptempdd(:,1).*Ptemp(:,2).*Ptemp(:,3);
			Pdd(:,k,2) = Ptempd(:,1).*Ptempd(:,2).*Ptemp(:,3);
			Pdd(:,k,3) = Ptempd(:,1).*Ptemp(:,2).*Ptempd(:,3);
			Pdd(:,k,4) = Ptemp(:,1).*Ptempdd(:,2).*Ptemp(:,3);
			Pdd(:,k,5) = Ptemp(:,1).*Ptempd(:,2).*Ptempd(:,3);
			Pdd(:,k,6) = Ptemp(:,1).*Ptemp(:,2).*Ptempdd(:,3);
		end
	end
	
	% Interpolate through the points in this subdomain
	F(J) = F(J) + ...
		w .* ( z + P*obj.weights{i}(Nc+1:end) );
	
	if nargout>2
		Fd(J,:) =  Fd(J,:) + ...
			wd .* ( z + P*obj.weights{i}(Nc+1:end) )+...
			w .* ( zd + [Pd(:,:,1)*obj.weights{i}(Nc+1:end),Pd(:,:,2)*obj.weights{i}(Nc+1:end),Pd(:,:,3)*obj.weights{i}(Nc+1:end)]);
	end
	
	if nargout>3
		doC1 = [1,1,1,2,2,3];
		doC2 = [1,2,3,2,3,3];
		Fdd(J,:) = Fdd(J,:) + wdd .* ( z + P*obj.weights{i}(Nc+1:end) )+w.*zdd +...
			wd(:,doC1).*zd(:,doC2)+wd(:,doC2).*zd(:,doC1);
		
		for j=1:6
			Fdd(J,j) = Fdd(J,j)+...
				wd(:,doC1(j)).*Pd(:,:,doC2(j))*obj.weights{i}(Nc+1:end)+...
				wd(:,doC2(j)).*Pd(:,:,doC1(j))*obj.weights{i}(Nc+1:end)+...
				w .* (Pdd(:,:,j)*obj.weights{i}(Nc+1:end));
		end
	end
	
	if exist('waitbarToggle','var') && waitbarToggle
		waitbar(i/M,wb,sprintf('%d out of %d patches done.',i,M))
	end
	
end

if nargout>2
	norms = (Fd.*wsum-F.*wsumd)./wsum.^2;
	normMag = sqrt(sum(norms.^2,2));
	unorms = norms./normMag;
end

if nargout>3
	doC1 = [1,1,1,2,2,3];
	doC2 = [1,2,3,2,3,3];
	Hdd = ((Fdd.*wsum+Fd(:,doC1).*wsumd(:,doC2)-Fd(:,doC2).*wsumd(:,doC1)-F.*wsumdd).*wsum.^2-...
		2*wsum.*wsumd(:,doC2).*(Fd(:,doC1).*wsum-F.*wsumd(:,doC1)))./wsum.^4;
	
	curv = 0.5*(-(sum(Hdd(:,[1,4,6]),2))./normMag+(norms(:,1).*sum(norms.*Hdd(:,[1,2,3]),2)+...
		norms(:,2).*sum(norms.*Hdd(:,[2,4,5]),2)+norms(:,3).*sum(norms.*Hdd(:,[3,5,6]),2))./normMag.^3);
end

F = F./wsum;
orphans = find(wsum==0);
F(orphans) = NaN;

if nargout>2
	unorms(orphans,:) = NaN;
end

if nargout>3
	curv(orphans) = NaN;
end

if exist('waitbarToggle','var') && waitbarToggle
	close(wb)
end

end

function out = phi(r)
out = zeros(size(r));
out(r<=1) = ( 1 - r(r<=1) ).^4 .* ( 4 .* r(r<=1) + 1 );
end

function out = phid(r)
out = zeros(size(r));
out(r<=1) = -4*( 1 - r(r<=1) ).^3 .* ( 4 .* r(r<=1) + 1 )+4*( 1 - r(r<=1) ).^4 ;
end

function out = phidd(r)
out = zeros(size(r));
out(r<=1) = 12*( 1 - r(r<=1) ).^2 .* ( 4 .* r(r<=1) + 1 )-32*( 1 - r(r<=1) ).^3;
end