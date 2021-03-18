function [weights, condition] = fitLocal( X, f, rbf, regFn, rho, polydegree )
%fitrbf Fits an approximate RBF interpolant.
%	Arguments:
%	X - datasites as row vectors
%	f - scalar values associated with datasites
%	polydegree - max polynomial degree
%	rbf - anonymous function
%	regFn - anonymous function
%	rho - smoothing parm passed to regFn
%
%	Solves:
%	[ rbf(distMatrix), P;] [ lambda ] = [ f ]
%	[       P',        0 ] [   a    ]   [ 0 ]

arguments
	X (:,:)
	f (:,1)
	rbf function_handle
	regFn function_handle
	rho (1,1)
	polydegree (1,1)
end

[N,~] = size(X);

[ xpowers, ypowers, zpowers ] = ...
	meshgrid( 0:polydegree, 0:polydegree, 0:polydegree );

powers = [ xpowers(:), ypowers(:), zpowers(:) ];

powers = powers( sum(powers,2) <= polydegree, : );
n = size(powers,1);

P = ones( N, n );
for k = 1:n
    for d = 1:3
        P(:,k) = P(:,k) .* (X(:,d)).^powers(k,d);
    end
end

A = rbf( surfFromPC.distMatrix(X,X) );

[Q,R] = qr(P);
Q1 = Q(:,1:n);
Q2 = Q(:,n+1:end);
R = R(1:n,:);

z = (Q2'*A*Q2 + regFn(N,rho)*eye(N-n)) \ (Q2'*f);
condition = condest(Q2'*A*Q2 + regFn(N,rho)*eye(N-n));
lambda = Q2*z;
a = R \ (Q1'*(f - A*lambda));

weights = [ lambda; a ];

end