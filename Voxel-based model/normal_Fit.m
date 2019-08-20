[~,normal_Parameters] = fit_Dist(NPs_per_Cell(cell_Locations),numel(x_Range{counter}),'normal');
normal_PDF{counter} = normpdf(x_Range{counter},normal_Parameters.mu,normal_Parameters.sigma);
normal_Error(counter) = sqrt(sum((normal_PDF{counter}-hist_Temp{counter}).^2));
normal_Parameters_Stored(counter,1) = normal_Parameters.mu;
normal_Parameters_Stored(counter,2) = normal_Parameters.sigma;