% Ldrum_modified.m - first K eigenvalues of Laplacian on L-shaped region
%
% This is a modified version of the "Ldrum" code by T. Betcke and
% L.N. Trefethen. The original computes the first few eigenvalues of
% the Laplacian on the L-shaped region by the method of particular
% solutions, using interior as well as boundary points.
%
% Original source:
%   T. Betcke and L. N. Trefethen, "Reviving the Method of Particular
%   Solutions," SIAM Review, vol. 47, no. 3, pp. 469-491, 2005.
%   DOI: 10.1137/S0036144503437336
%
% Modifications (relative to the original):
%   - parametrized by K, the number of eigenvalues wanted;
%   - the trial range lammax is sized automatically from Weyl's law
%     so that it is guaranteed to contain at least K eigenvalues;
%   - the accuracy parameter N grows with lammax so higher modes
%     stay resolved;
%   - the computed eigenvalues are collected, sorted, and the first K
%     are returned;
%   - the random number generator is seeded so that the random interior
%     collocation points, and hence the output, are reproducible.

% --- Parameters ---
K = 8; % number of eigenvalues wanted

rng(0); % seed the RNG for reproducible interior collocation points

% Upper bound on lambda guaranteed to contain at least K eigenvalues.
% Weyl's law for the L-shape (area = 3): N(lambda) ~ area*lambda/(4*pi),
% so lambda_K ~ 4*pi*K/area.  A safety margin is added.
area = 3;
lammax = (4*pi/area)*K*1.5 + 20; % trial range upper limit

% Accuracy parameter, grown with lammax so high modes stay resolved.
N = max(36, ceil(6*sqrt(lammax))); % no. of terms in Bessel expansion
np = 2*N; % no. of boundary and interior pts
k = 1:N; % orders in Bessel expansion
t1 = 1.5*pi*(.5:np-.5)'/np; % angles of bndry pts
r1 = 1./max(abs(sin(t1)),abs(cos(t1))); % radii of bndry pts
t2 = 1.5*pi*rand(np,1); % angles of interior pts
r2 = rand(np,1)./max(abs(sin(t2)),abs(cos(t2))); % radii interior pts
t = [t1;t2]; r = [r1;r2]; % bndry and interior combined
lamvec = .2:.2:lammax; % trial values of lambda

% Compute subspace angles for various values of lambda:
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
K0 = 0*I;
K0(J) = 1;
S = S.*(-1).^cumsum(K0); % introduce sign flips

hold off, plot(lamvec,S), hold on % plot signed angle function
plot([0 max(lamvec)],[0 0],'-k') % plot lambda axis

% Find eigenvalues via local 9th-order interpolation:
evals = zeros(1,length(J));
for j = 1:length(J)
I = J(j)-5:J(j)+4;
evals(j) = polyval(polyfit(S(I)/norm(S(I)),lamvec(I),9),0);
end
evals = sort(evals); % ascending order

if length(evals) < K
warning('Only %d eigenvalues found below lambda = %.1f; increase K-range margin.', ...
        length(evals), lammax);
K = length(evals);
end
evals = evals(1:K); % keep the first K

disp(evals') % display the first K eigenvalues
for j = 1:K
plot(evals(j)*[1 1],.8*[-1 1],'r') % plot eigenvalue
end
