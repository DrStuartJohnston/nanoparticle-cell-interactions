number_Down_Events = poissrnd(dt*propensity_Down_Matrix.*NP_Domain);                                %Number of transitions down for each voxel
number_Left_x_Events = poissrnd(dt*propensity_Left_x_Matrix.*NP_Domain);                            %Number of transitions left (x) for each voxel
number_Right_x_Events = poissrnd(dt*propensity_Right_x_Matrix.*NP_Domain);                          %Number of transitions right (x) for each voxel
number_Left_y_Events = poissrnd(dt*propensity_Left_y_Matrix.*NP_Domain);                            %Number of transitions left (y) for each voxel
number_Right_y_Events = poissrnd(dt*propensity_Right_y_Matrix.*NP_Domain);                          %Number of transitions right (y) for each voxel
number_Up_Events = poissrnd(dt*propensity_Up_Matrix.*NP_Domain);                                    %Number of transitions up for each voxel
tmprng = rng;                                                                                       %Save rng seed
number_Bind_Events_Scaled = poissrnd(dt*propensity_Bind_Matrix.*NP_Domain/IBM_Parameters.scaling);  %Number of bind events (scaled) for each voxel
rng(tmprng);                                                                                        %Reset rng seed
number_Bind_Events = poissrnd(dt*propensity_Bind_Matrix.*NP_Domain);                                %Number of bind events (unscaled - rngset ensures same as above)