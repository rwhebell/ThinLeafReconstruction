function [z,zd,zdd] = distMatrixAct(x,y,v,rbf,rbfd,rbfdd)
%distMatrixAct actions a distance matrix of x and y onto v
%	z = distMatrixAct(x,y,v)

wantNumEl = 1e8;

Nx = size(x,1);
Ny = size(y,1);

assert(Ny==length(v))

Nb = ceil( wantNumEl / Ny ); % how many in one block

M = ceil(Nx/Nb); % how many blocks

z = zeros( Nx, 1, 'like', x );
if nargout>1
    zd = zeros( Nx, 3, 'like', x );
end
if nargout>2
    zdd = zeros( Nx, 6, 'like', x );
end
for b = 0:M-1
	
	B = b*Nb+1:min((b+1)*Nb,Nx);
	
    if nargout==1
        D = surfFromPC.distMatrix( x( B, : ), y );
    elseif nargout==2
        [D,Dd] = surfFromPC.distMatrix( x( B, : ), y );
        Dd(isnan(Dd))=0;
    elseif nargout==3
        [D,Dd,Ddd] = surfFromPC.distMatrix( x( B, : ), y );
    end
	
	z(B) = rbf(D)*v;
    if nargout>1
        for i=1:3
            zd(B,i) = (rbfd(D).*Dd(:,:,i))*v;
        end
    end
    if nargout>2
        doC1 = [1,1,1,2,2,3];
        doC2 = [1,2,3,2,3,3];
        for i=1:6
            zdd(B,i) = (rbfdd(D).*Dd(:,:,doC1(i)).*Dd(:,:,doC2(i))+rbfd(D).*Ddd(:,:,i))*v;
        end
    end
	
end