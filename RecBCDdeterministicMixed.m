%% RecBCDdeterministic(p,data,chi_sequence)..

% given the set of parameters (input: p) and
% a sequence of locations of Chi sites (input: chi_sequence)
% computes the PMF of a genomic location,
% positioned x nucleotides away from DSB,
% to belong to ssDNA (output: PMF)
% It also calculates log-liklihood (output: L) of the data (input: data),
% given the model parameter set (input: p)

function [L] = RecBCDdeterministicMixed(p,data,ChiPositions)%#codegen
c = @cmu.colors;

TOTALL = length(data);

x = 1:TOTALL;
PMF = zeros(1,TOTALL);
   ChiPositionsL = double(ChiPositions);
for j = 1 :length(ChiPositions)
    
    PMF = PMF + (1-p(6))*p(4)*(1-p(4)).^(j-1).*(1-p(1))...
       .^(max(0, x - (ChiPositions(j))/p(5))).*heaviside(x-ChiPositions(j)) + p(6)*p(2)*(1-p(2)).^(j-1).*(1-p(1))...
       .^(max(0, x - (ChiPositions(j))/p(3))).*heaviside(x-ChiPositions(j)) ;
%+ 
            
end

% PMF normalised 
PMF = (PMF(ChiPositions(1):TOTALL) )/sum(PMF(ChiPositions(1):TOTALL));

% correction for the case when backgriound is subtracted prior to
% optimisation
OL = data(ChiPositions(1):TOTALL);

L = sum(log(PMF(find(OL > 0))).*OL(find(OL > 0)));


end
