function laplace_valsvec = dst_laplace(valsvec, ndgrid)

% Now we reshape the input vector to an array that has the same shape as
% the x-y grid
valsgrid = reshape(valsvec, size(X));

% These are grid values of the second derivative with respect to x, scaling
% factor is included

[m, n] = size(X);

valsgridD2XX = zeros(m, n);   
parfor j = 1:n
    valsgridD2XX(:, j) = dst_derivative_slice(valsgrid(:, j), a, b);
end

valsgridD2YY = zeros(m, n);   
parfor i = 1:m
    valsgridD2YY(i, :) = dst_derivative_slice(valsgrid(i, :)', c, d)';
end


% These are grid values of the second derivative with respect to y, scaling
% factor is included
%valsgridD2YY = @(x) dst_derivative_slice(x', c, d)'; 

laplace_valsgrid = valsgridD2XX  + valsgridD2YY;
laplace_valsvec = laplace_valsgrid(:);

end