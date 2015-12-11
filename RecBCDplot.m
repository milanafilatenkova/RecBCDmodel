%% RecBCDdeterministic(p,data,chi_sequence)..

% given the set of parameters (input: p) and
% a sequence of locations of Chi sites (input: chi_sequence)
% computes the PMF of a genomic location,
% positioned x nucleotides away from DSB,
% to belong to ssDNA (output: PMF)
% It also calculates log-liklihood (output: L) of the data (input: data),
% given the model parameter set (input: p)

function [] = RecBCDplot(p, data,ChiPositions,w)
c = @cmu.colors;


TOTALL = length(data);
x = 1:TOTALL;
PMF = zeros(1,TOTALL);
   
for j = 1 :length(ChiPositions)
    
    PMF = PMF + p(2)*(1-p(2)).^(j-1).*(1-p(1))...
       .^(max(0, x - (ChiPositions(j))/p(3))).*heaviside(x-ChiPositions(j)) ;
            
end

% PMF normalised 
PMF = (PMF)/sum(PMF);


% Log-likelihood
OL = data(1:TOTALL);
x = [1: w: TOTALL*w];


figure, hold on,
plot(x,data/sum(data),'Color', c('silver'))
plot(x,PMF,'Color', c('navy blue'),'LineWidth',1.2)
plot(x, smooth(data/sum(data),0.057, 'loess'), 'Color', c( 'raspberry' ))
plot(ChiPositions*w,0,'o','MarkerSize', 10,'MarkerFaceColor',	c('aquamarine'),'MarkerEdgeColor',	'k' )


end
