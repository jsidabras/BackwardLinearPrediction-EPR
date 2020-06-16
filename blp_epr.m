function [blp_out, blp_only] = blp_epr(y,L,n)
% BLP-EPR backward linear prediction algorithm for 3Pulse DEER data
% function [blp_out, blp_only] = blp_epr(y,L,n) 
%
% arguments: 
%   y: real vector; data input currently does not support complex data.
%      Remember to phase and background correct *before* running blp_epr
%   L: real integer; number of points for backwards prediction     
%   n: real integer; number of coefficients for linear prediction 
%      (default: 25)
%   
% outputs: 
%   blp_out: full vector of input data y with M concatinated predicted points
%   blp_only: only M concatinated predicted points
%
% Author: Jason W. Sidabras (jason.sidabras@gmail.com)
%   Initial writing: 04/05/2020 JWS
%     New Algorithm: 16/06/2020 JWS
%   GPLv3 License.

% Default for n is 25. This gives a reasonable set of statistics during testing.
if nargin<3
  n = 25;
end

% currently the toeplitz matrix is setup for forward linear prediction
% just flip the input data and make sure it is real
% TODO: What happends if it is complex?
y = flip(real(y));

% Create a Toeplitz matrix
% this uses the entire input vector 
% TODO: it may be adventageous to limit the input length to filter out
% known spurrious signals (2+1 artifacts).
N = length(y);
H = toeplitz(y(n:N-1), y(n:-1:1));    
b = y(n+1:N)';

% Using the lsqminnorm really cleaned up the statistics. Much better solver
% for this set of problems.
a = lsqminnorm(H,b);
g = [y'; zeros(L, 1)];

% Use the solved vector to predict L points. 
for i = N+1:N+L
    g(i) = sum(a.*flip(g(i-n:i-1)));
% the above equation creates a dynamic version of the equation below. 
% this way one can choose the number of coefficients
%     g(i) = a(1) * g(i-1) + a(2) * g(i-2) + a(3) * g(i-3) + a(4) * g(i-4) + a(5) * g(i-5) ...
%         + a(6) * g(i-6) + a(7) * g(i-7) + a(8) * g(i-8)+ a(9) * g(i-9) + a(10) * g(i-10) ...
%         + a(11) * g(i-11) + a(12) * g(i-12) + a(13) * g(i-13) + a(14) * g(i-14) ...
%         + a(15) * g(i-15) + a(16) * g(i-16) + a(17) * g(i-17) + a(18) * g(i-18) ...
%         + a(19) * g(i-19) + a(20) * g(i-20) + a(21) * g(i-21) + a(22) * g(i-22) ...
%         + a(23) * g(i-23) + a(24) * g(i-24) + a(25) * g(i-25);
end

% output full data
g = flip(g);
blp_out = g;

% get rid of original input data and output only the BLP values
g(L+1:end) = [];
blp_only = g;
end