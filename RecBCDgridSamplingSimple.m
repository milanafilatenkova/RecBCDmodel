function [L,T,J,K, ps_max,pchi_max,pm_max] = ...
     RecBCDgridSamplingSimple(ps_range,pchi_range,pm_range,data,ChiPositions);
 
L = [];


%compute log-likelohood on the Grid 
for t = 1:length(ps_range) 
        for j = 1:length(pchi_range)
             for k = 1:length(pm_range)
                 
                  [L(t,j,k)] = RecBCDdeterministic...
                          ([ps_range(t),pchi_range(j),pm_range(k)],data,ChiPositions);
                
             end
        end
end

%identify the point on the Grid where log-likelihood reaches its
%maximum

ind = find(L == max(max(max(max(L)))));
ind
[T,J,K] = ind2sub(size(L),ind);

ps_max = ps_range(T);
pchi_max = pchi_range(J);
pm_max = pm_range(K);



