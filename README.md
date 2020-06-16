# BackwardLinearPrediction-EPR
BLP-EPR backward linear prediction algorithm for 3Pulse DEER data
function [blp_out, blp_only] = blp_epr(y,L,n) 

### arguments:

   **y**: real vector; data input currently does not support complex data.
      Remember to phase and background correct *before* running blp_epr

   **L**: real integer; number of points for backwards prediction     

   **n**: real integer; number of coefficients for linear prediction 
      (default: 25)
   
### outputs: 
   **blp_out**: full vector of input data y with M concatinated predicted points

   **blp_only**: only M concatinated predicted points

### Credits
 Author: Jason W. Sidabras (jason.sidabras@gmail.com)

   Initial writing: 04/05/2020 JWS

     New Algorithm: 16/06/2020 JWS

   GPLv3 License.
