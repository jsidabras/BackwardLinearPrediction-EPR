function blp_out = blp_epr(y,n,q,M)
% BLP-EPR backward linear prediction algorithm for 3Pulse DEER data
% function blp_out = blp_epr(y,n,q,M) 
%
% arguments: 
%   y: real vector; data input 
%   n: number of rows in n x q matrix for linear solution
%   q: number of columns in n x q matrix for linear solution
%   M: real scalar, number of points for backwards prediction
%   blp_out: output of the M number of points 
%
% Author: Jason W. Sidabras (jason.sidabras@gmail.com)
%   Initial writing: 04/05/2020 JWS
%   GPLv3 License.

% Initialize TMatrix and pVector
TMatrix = [];
pVector = [];
% initialize the answer vector
blp_out = [];
for j = 1:n
    % set up the TMatrix and pVector
    % the TMatrix is a n x q matrix that steps through the data
    % in order to predict one value right before the data starts
    TMatrixtmp = y(j+1:j+q);
    TMatrix = [TMatrix; TMatrixtmp];
    % the pVector is a vector of size n that seeds the predicition
    pVectortmp = y(j);
    pVector = [pVector; pVectortmp];
end

for p = 1:M
    % clear and initialize the next TVector and predicted p value
    NextTVector = [];
    NextTMatrix = [];

    % solve M.x = b, where x describes the system of knowns
    xSolved = lsqminnorm(TMatrix, pVector);
 
    % create an TVector that is missing its pPoint
    NextTVector = TMatrix(1:length(TMatrix(1,:))-1);
    NextTVector = [pVector(1), NextTVector];
    
    % find the backward predicted point
    pPoint = sum(NextTVector.*xSolved');
    
    % save predicted point in a solution vector
    blp_out = [pPoint, blp_out];
    
    % prepend new pPoint to pVector and delete last point keeping size n
    pVector = [pPoint; pVector];
    pVector(end,:) = [];
    
    % prepend new TVector to TMatrix and delete last vector 
    % keeping size n x q
    TMatrix(end,:) = [];
    TMatrix = [NextTVector; TMatrix];
    
end