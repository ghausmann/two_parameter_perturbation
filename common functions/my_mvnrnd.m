function X = my_mvnrnd(mu, Sigma, n)
    % my_mvnrnd generates n samples from a multivariate normal distribution
    % with mean vector mu and covariance matrix Sigma.
    
    % Ensure mu is a column vector
    mu = mu(:);
    d = length(mu);
    
    % Compute the lower triangular Cholesky factor
    A = chol(Sigma, 'lower');
    
    % Generate n samples of independent standard normal variables
    Z = randn(n, d);
    
    % Transform the samples to have the desired covariance structure
    X = repmat(mu, 1, n) + A*Z';
    X = X';
end