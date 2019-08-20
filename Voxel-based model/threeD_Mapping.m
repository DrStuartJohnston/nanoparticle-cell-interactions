%% Calculate index in vector corresponding to (i,j,k) co-ordinates

function index = threeD_Mapping(row,column,depth,M,N,P)
% N is matrix height (Parameters.grid_Height), M is matrix width
% (Parameters.cell_Spaces_x), P is matrix depth (Parameters.cell_Spaces_y)

column(column==0) = M;
column(column==M+1) = 1;
depth(depth==0) = P;
depth(depth==P+1) = 1;
row(row==0) = N;
row(row==N+1) = 1;

index = (depth-1)*N*M+(column-1)*N+row;