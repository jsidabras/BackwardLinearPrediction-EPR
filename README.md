# BackwardLinearPrediction-EPR
A backward linear prediction algorithm implemented for 3Pulse DEER data

function blp_out = blp_epr(y,n,q,M) 

arguments: 
  y: real vector; data input 
  n: number of rows in n x q matrix for linear solution
  q: number of columns in n x q matrix for linear solution
  M: real scalar, number of points for backwards prediction
  blp_out: output of the M number of points 
