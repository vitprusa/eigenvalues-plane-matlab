% Ldrum.m - eigenvalues of Laplacian on L-shaped region
%
% The first 3 eigenvalues are computed by the method of
% particular solutions, using interior as well as boundary
% points. T. Betcke and L.N. Trefethen, 2003.

% Compute subspace angles for various values of lambda:

N = 36; % accuracy parameter
np = 2*N; % no. of boundary and interior pts
k = 1:N; % orders in Bessel expansion
t1 = 1.5*pi*(.5:np-.5)'/np; % angles of bndry pts
r1 = 1./max(abs(sin(t1)),abs(cos(t1))); % radii of bndry pts
t2 = 1.5*pi*rand(np,1); % angles of interior pts
r2 = rand(np,1)./max(abs(sin(t2)),abs(cos(t2))); % radii interior pts
t = [t1;t2]; r = [r1;r2]; % bndry and interior combined
lamvec = .2:.2:25; % trial values of lambda

S = [];
for lam = lamvec
A = sin(2*t*k/3).*bsxfun(@besselj, 2*k/3,sqrt(lam)*r); % CORRECTION TO NEWER MATLAB VERSION
[Q,R] = qr(A,0);
s = min(svd(Q(1:np,:))); % subspace angle for this lam
S = [S s];
end

% Convert to signed subspace angles:
I = 1:length(lamvec); % all lam points
J = I(2:end-1); % interior points
J = J( S(J)<S(J-1) & S(J)<S(J+1) ); % local minima
J = J + (S(J-1)>S(J+1)); % points where sign changes
K = 0*I; 
K(J) = 1;
S = S.*(-1).^cumsum(K); % introduce sign flips

hold off, plot(lamvec,S), hold on % plot signed angle function
plot([0 max(lamvec)],[0 0],'-k') % plot lambda axis


% Find eigenvalues via local 9th-order interpolation:
for j = 1:length(J)
I = J(j)-5:J(j)+4;
lambda = polyval(polyfit(S(I)/norm(S(I)),lamvec(I),9),0);
disp(lambda) % display eigenvalue
plot(lambda*[1 1],.8*[-1 1],'r') % plot eigenvalue
end