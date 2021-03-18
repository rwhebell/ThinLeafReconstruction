function [D,Dd,Ddd] = distMatrix(x,y)
%distMatrix
%   D(i,j) = || x_i - y_j ||_2 ^ 2
%   D(i,j) = ||x_i||^2 + ||y_j||^2 - 2(x_i dot y_j)

D = pdist2(x,y,"euclidean");

if nargout>1
    Dd = zeros(size(x,1),size(y,1),3);
    for i=1:3
        Dd(:,:,i) = (x(:,i)-y(:,i)')./D;
    end
end

if nargout>2
    diC = [1,0,0,1,0,1];
    doC1 = [1,1,1,2,2,3];
    doC2 = [1,2,3,2,3,3];
    D3 = D.^3;
    Ddd = zeros(size(x,1),size(y,1),6);
    for i=1:6
        Ddd(:,:,i) = diC(i)./D-(x(:,doC1(i))-y(:,doC1(i))').*(x(:,doC2(i))-y(:,doC2(i))')./D3;
    end
end