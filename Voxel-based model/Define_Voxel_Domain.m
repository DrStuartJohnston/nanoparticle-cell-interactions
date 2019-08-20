%% Script to define transitions and [ropensities in the voxel domain

top_Indices = 1:IBM_Parameters.grid_Height:IBM_Parameters.grid_Height*IBM_Parameters.total_Cell_Spaces;                                 %Vector indices corresponding to top boundary
bottom_Indices = IBM_Parameters.grid_Height:IBM_Parameters.grid_Height:IBM_Parameters.grid_Height*IBM_Parameters.total_Cell_Spaces;     %Vector indices corresponding to bottom boundary 
middle_Indices = setdiff(1:IBM_Parameters.grid_Height*IBM_Parameters.total_Cell_Spaces,[top_Indices,bottom_Indices]);                   %Vector indices corresponding to remainder of domain    

propensity_Down_Matrix = (IBM_Parameters.D/IBM_Parameters.grid_Size^2+IBM_Parameters.V/(2*IBM_Parameters.grid_Size))*ones(size(NP_Domain)); %Propensity for transitions to voxel below
propensity_Down_Matrix(bottom_Indices) = 0;                                                                                                 %Propensity for transitions to voxel below on bottom boundary              
propensity_Left_x_Matrix = (IBM_Parameters.D/IBM_Parameters.grid_Size^2)*ones(size(NP_Domain));                                             %Propensity for transitions to voxel left (x)
propensity_Right_x_Matrix = (IBM_Parameters.D/IBM_Parameters.grid_Size^2)*ones(size(NP_Domain));                                            %Propensity for transitions to voxel right (x)
propensity_Left_y_Matrix = (IBM_Parameters.D/IBM_Parameters.grid_Size^2)*ones(size(NP_Domain));                                             %Propensity for transitions to voxel left (y)
propensity_Right_y_Matrix = (IBM_Parameters.D/IBM_Parameters.grid_Size^2)*ones(size(NP_Domain));                                            %Propensity for transitions to voxel right (y)
propensity_Up_Matrix = (IBM_Parameters.D/IBM_Parameters.grid_Size^2-IBM_Parameters.V/(2*IBM_Parameters.grid_Size))*ones(size(NP_Domain));   %Propensity for transitions to voxel up
propensity_Up_Matrix(top_Indices) = 0;                                                                                                      %Propensity for transitions to voxel up on top boundary
propensity_Bind_Matrix = (IBM_Parameters.D/IBM_Parameters.grid_Size^2+IBM_Parameters.V/(2*IBM_Parameters.grid_Size))*ones(size(NP_Domain)); %Propensity for binding to cell layer
propensity_Bind_Matrix([top_Indices,middle_Indices]) = 0;                                                                                   %Propensity for binding for non-bottom voxels

index = 1:IBM_Parameters.total_Cell_Spaces*IBM_Parameters.grid_Height;                      %Indices for voxels
index_Down = zeros(IBM_Parameters.total_Cell_Spaces*IBM_Parameters.grid_Height,1);          %Corresponding voxels down
index_Left_x = zeros(IBM_Parameters.total_Cell_Spaces*IBM_Parameters.grid_Height,1);        %Corresponding voxels left (x)
index_Right_x = zeros(IBM_Parameters.total_Cell_Spaces*IBM_Parameters.grid_Height,1);       %Corresponding voxels right (x)
index_Left_y = zeros(IBM_Parameters.total_Cell_Spaces*IBM_Parameters.grid_Height,1);        %Corresponding voxels left (y)
index_Right_y = zeros(IBM_Parameters.total_Cell_Spaces*IBM_Parameters.grid_Height,1);       %Corresponding voxels right (y)
index_Up = zeros(IBM_Parameters.total_Cell_Spaces*IBM_Parameters.grid_Height,1);            %Corresponding voxels up