function[data_binned] = data_bin(data, w)

for i = 1: floor(length(data)/w)
data_binned(i) = sum(data((i-1)*w +1 :i*w))/w;
end
