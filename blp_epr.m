function [blp_out, blp_only] = blp_epr(y,L,n)
% BLP-EPR backward linear prediction algorithm for 3Pulse DEER data
% function [blp_out, blp_only] = blp_epr(y,L,n) 
%
% arguments: 
%   y: complex vector; data input currently does not support complex data.
%      Remember to phase and background correct *before* running blp_epr
%   L: real integer; number of points for backwards prediction     
%   n: real integer; number of coefficients for linear prediction 
%      (default: 50)
%   
% outputs: 
%   blp_out: full vector of input data y with M concatinated predicted points
%   blp_only: only M concatinated predicted points
%
% Author: Jason W. Sidabras (jason.sidabras@gmail.com)
%   Initial writing: 04/05/2020 JWS
%     New Algorithm: 16/06/2020 JWS
%   GPLv3 License.

    % When n is 50: This gives a reasonable set of statistics 
    % during testing without background.
    % Default is maximum number of points: This works really well 
    % with data with real-world backgrounds and matches noise level.
    % 
    if nargin<3
      n = length(y)-50 ;
    end

    % currently the toeplitz matrix is setup for forward linear prediction
    % just flip the input data 
    y = flip(y')';

    % Create a Toeplitz matrix
    % this uses the entire input vector 
    % TODO: it may be adventageous to limit the input length to filter out
    % known spurrious signals (2+1 artifacts).
    N = length(y);
    H = toeplitz(y(n:N-1), y(n:-1:1));   
    b = y(n+1:N)';

    % Using the lsqminnorm really cleaned up the statistics. Much better solver
    % for this set of problems.
    % Use the solved vector to predict L points.
    if isreal(y)
        a = lsqminnorm(H,b);
        g = [y'; zeros(L, 1)];
        for i = N+1:N+L
            g(i) = sum(a.*flip(g(i-n:i-1)));
        end
        g = flip(g);
    else
        % With this algorithm it is better to have real/imag separate
        % Real Channel
        y_r = real(y);
        H_r = toeplitz(y_r(n:N-1), y_r(n:-1:1));    
        b_r = y_r(n+1:N)';
        a_r = lsqminnorm(H_r,b_r);       
        g_r = [y_r'; zeros(L, 1)];
        
        % Imaginary Channel
        y_i = imag(y);
        H_i = toeplitz(y_i(n:N-1), y_i(n:-1:1));    
        b_i = y_i(n+1:N)';
        a_i = lsqminnorm(H_i,b_i);       
        g_i = [y_i'; zeros(L, 1)];
        % complex answer
        for i = N+1:N+L
            g_r(i) = sum(a_r.*flip(g_r(i-n:i-1)));
            g_i(i) = sum(a_i.*flip(g_i(i-n:i-1)));
        end
        g = complex(flip(g_r),flip(g_i));
    end
    % NOTE:  
    %    g(i) = sum(a.*flip(g(i-n:i-1)));
    % the above equation creates a dynamic version of the equation below. 
    % this way one can choose the number of coefficients
    %     g(i) = a(1) * g(i-1) + a(2) * g(i-2) + a(3) * g(i-3) + a(4) * g(i-4) + a(5) * g(i-5) ...
    %         + a(6) * g(i-6) + a(7) * g(i-7) + a(8) * g(i-8)+ a(9) * g(i-9) + a(10) * g(i-10) ...
    %         + a(11) * g(i-11) + a(12) * g(i-12) + a(13) * g(i-13) + a(14) * g(i-14) ...
    %         + a(15) * g(i-15) + a(16) * g(i-16) + a(17) * g(i-17) + a(18) * g(i-18) ...
    %         + a(19) * g(i-19) + a(20) * g(i-20) + a(21) * g(i-21) + a(22) * g(i-22) ...
    %         + a(23) * g(i-23) + a(24) * g(i-24) + a(25) * g(i-25);


    % output full data
    blp_out = g;

    % get rid of original input data and output only the BLP values
    g(L+1:end) = [];
    blp_only = g;
end