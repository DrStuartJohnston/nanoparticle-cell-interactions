function [x,y,z] = reverse_threeD_Mapping(index,M,N)

x = ceil((mod(index-1,M*N)+1)/N);
y = ceil(index/(M*N));
z = mod(index-1,N)+1;