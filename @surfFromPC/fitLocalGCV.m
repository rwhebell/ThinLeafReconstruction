function [weights, rho, condition] = fitLocalGCV( X, f, rbf, regFn, rhoRange, polydegree )
%fitrbf Fits an approximate RBF interpolant.
%	Arguments:
%	X - datasites as row vectors
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
	rhoRange (1,:)
	polydegree (1,1)
end

[N,~] = size(X);

[ xpowers, ypowers, zpowers ] = ...
	meshgrid( 0:polydegree, 0:polydegree, 0:polydegree );

powers = [ xpowers(:), ypowers(:), zpowers(:) ];

powers = powers( sum(powers,2) <= polydegree, : );
n = size(powers,1);

P = ones( N, n );
for k = 1:size(powers,1)
	for d = 1:3
		P(:,k) = P(:,k) .* (X(:,d)).^powers(k,d);
	end
end

A = rbf( surfFromPC.distMatrix(X,X) );

[Q,R] = qr(P);
Q1 = Q(:,1:n);
Q2 = Q(:,n+1:end);
R = R(1:n,:);

options = optimset('TolX',0.1);

% GCV
fun = @gcvFn;
[logrho_opt, V, flag, output] = ...
	fminbnd(fun, log10(rhoRange(1)), log10(rhoRange(2)), options);

if flag == 1
	fprintf('Opt rho = %.2e, with V(rho) = %.2e, in %d func evals.\n',...
		10^logrho_opt, V, output.funcCount )
elseif flag == 0
	error('fminbnd exceeded max iters or max function count.')
end

rho = 10^(logrho_opt);  % optimal rho value

z = (Q2'*A*Q2 + regFn(N,rho)*eye(N-n)) \ (Q2'*f);
condition = condest(Q2'*A*Q2 + regFn(N,rho)*eye(N-n));
lambda = Q2*z;
a = R \ (Q1'*(f - A*lambda));

weights = [ lambda; a ];

	function V = gcvFn(logrho)
		V = 0*logrho;
		for i = 1:numel(logrho)
			IminusB = regFn(N,10^logrho(i))*...
				Q2*inv(Q2'*(A+regFn(N,10^logrho(i))*eye(N))*Q2)*Q2';
			V(i) = N*norm(IminusB*f)^2 / trace(IminusB)^2;
		end
	end


end