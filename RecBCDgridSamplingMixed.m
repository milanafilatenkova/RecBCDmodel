
%% RecBCD_grid_sampling(n_chi,ps_range,pchi_range,pm_range)

% computes log-likelohood (output: L)in each point of the input Grid;
% estimates the optimal set of parameters coresponding
% .. to the maximum of log-likelehood function on the course grained
% .. parameter space (output:  ps_max,pchi_max,pm_max, pchi_max2,pm_max2,v_max);

% the function takes as input the number of Chi sites in the Chi-array
% .. (input: n_chi = 1..6);
% .. and the Grid definition (input: ps_range,pchi_range,pm_range);

function [L,T,J,K,JJ,KK,N, ps_max,pchi_max1,pm_max1,pchi_max2,pm_max2,a_max]...
    = RecBCDgridSamplingMixed(ps_range,pchi_range1,pm_range1,pchi_range2,pm_range2,a_range, data, ChiPositions)


L=zeros(length(ps_range),length(pchi_range1),length(pm_range1),length(pchi_range2),length(pm_range2), length(a_range));

%compute log-likelohood on the Grid 
for t = 1:length(ps_range) 
        for j = 1:length(pchi_range1)
             for k = 1:length(pm_range1)
                   for jj = 1:length(pchi_range2)
                      for kk = 1:length(pm_range2)
                        for  n =1:length(a_range)
                       
                           [L(t,j,k,jj,kk,n)] =  RecBCDdeterministicMixed...
                          ([ps_range(t),pchi_range1(j),pm_range1(k),pchi_range2(jj),pm_range2(kk), a_range(n)],data,ChiPositions);
                        end
                      end
                   end
             end
    end   
end

%identify the point on the Grid where log-likelihood reaches its
%maximum
ind = find(L == max(max(max(max(max(max(max(L))))))));
[T,J,K,JJ,KK,N] = ind2sub(size(L),ind);

ps_max = ps_range(T);
pchi_max1 = pchi_range1(J);
pm_max1 = pm_range1(K);
pchi_max2 = pchi_range2(JJ);
pm_max2 = pm_range2(KK);
a_max = a_range(N);
end


