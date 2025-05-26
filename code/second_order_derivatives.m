function [Dxv,Dyv,Dxxv,Dyyv,Dxyv] = second_order_derivatives(N,M,tx,v)
% Define grid dimensions and step sizes
% N = 100; % Number of grid points in x-direction
% M = 100; % Number of grid points in y-direction
% tx = 1; % Step size in x-direction
% ty = 1; % Step size in y-direction
ty = tx;

% Initialize the function v (example: random values)
%v = rand(N, M);

% Preallocate the differential matrices
Dxv = zeros(N, M);
Dyv = zeros(N, M);
Dxxv = zeros(N, M);
Dyyv = zeros(N, M);
Dxyv = zeros(N, M);

% Compute first order derivatives
for i = 1:N-1
    for j = 1:M-1
        Dxv(i,j) = (v(i+1,j) - v(i,j)) / tx;
        Dyv(i,j) = (v(i,j+1) - v(i,j)) / ty;
    end
end

% Compute second order derivatives
for i = 2:N-1
    for j = 2:M-1
        Dxxv(i,j) = (v(i+1,j) - 2*v(i,j) + v(i-1,j)) / tx^2;
        Dyyv(i,j) = (v(i,j+1) - 2*v(i,j) + v(i,j-1)) / ty^2;
        Dxyv(i,j) = (1/ty) * ((v(i+1,j+1) - v(i,j+1))/tx - (v(i+1,j) - v(i,j))/tx);
    end
end