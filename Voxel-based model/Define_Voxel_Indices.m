%% Generate the indices corresponding to transitions in down, left (x), right (x), left (y), right (y) and up directions

for i = 1:IBM_Parameters.cell_Spaces_y
    for j = 1:IBM_Parameters.cell_Spaces_x
        for k = 1:IBM_Parameters.grid_Height
            count = (i-1)*IBM_Parameters.cell_Spaces_x*IBM_Parameters.grid_Height + (j-1)*IBM_Parameters.grid_Height + k;
            index_Down(count) = threeD_Mapping(k+1,j,i,IBM_Parameters.cell_Spaces_x,IBM_Parameters.grid_Height,IBM_Parameters.cell_Spaces_y);
            index_Left_x(count) = threeD_Mapping(k,j-1,i,IBM_Parameters.cell_Spaces_x,IBM_Parameters.grid_Height,IBM_Parameters.cell_Spaces_y);
            index_Right_x(count) = threeD_Mapping(k,j+1,i,IBM_Parameters.cell_Spaces_x,IBM_Parameters.grid_Height,IBM_Parameters.cell_Spaces_y);
            index_Left_y(count) = threeD_Mapping(k,j,i-1,IBM_Parameters.cell_Spaces_x,IBM_Parameters.grid_Height,IBM_Parameters.cell_Spaces_y);
            index_Right_y(count) = threeD_Mapping(k,j,i+1,IBM_Parameters.cell_Spaces_x,IBM_Parameters.grid_Height,IBM_Parameters.cell_Spaces_y);
            index_Up(count) = threeD_Mapping(k-1,j,i,IBM_Parameters.cell_Spaces_x,IBM_Parameters.grid_Height,IBM_Parameters.cell_Spaces_y);
        end
    end
end