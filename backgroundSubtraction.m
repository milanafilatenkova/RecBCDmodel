 function [data_B] = backgroundSubtraction(data, background_region, w)
 
 data_B = data - mean(data(1:background_region/w));
 
 %
 