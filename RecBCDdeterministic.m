%% RecBCDdeterministic(p,data,chi_sequence)..

% given the set of parameters (input: p) and
% a sequence of locations of Chi sites (input: chi_sequence)
% computes the PMF of a genomic location,
% positioned x nucleotides away from DSB,
% to belong to ssDNA (output: PMF)
% It also calculates log-liklihood (output: L) of the data (input: data),
% given the model parameter set (input: p)

function [L] = RecBCDdeterministic(p,data,ChiPositions)

% initialise random variable x
% which represents genomic location 
% as distance from DSB 

TOTALL = length(data);
x = 1:TOTALL;
PMF = zeros(1,TOTALL);
   
for j = 1 :length(ChiPositions)
    
    PMF = PMF + p(2)*(1-p(2)).^(j-1).*(1-p(1))...
       .^(max(0, x - (ChiPositions(j))/p(3))).*heaviside(x-ChiPositions(j));
            
end

% PMF normalised 
PMFL = (PMF(ChiPositions(1):TOTALL))/(sum(PMF(ChiPositions(1):TOTALL)));

% Data
OL = data(ChiPositions(1):TOTALL);

% Log-Likelihood
L = sum(log(PMFL(find(OL > 0))).*OL(find(OL > 0)));

end
