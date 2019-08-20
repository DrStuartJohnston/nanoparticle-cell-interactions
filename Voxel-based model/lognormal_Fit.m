[~,lognormal_Parameters] = fit_Dist(NPs_per_Cell(cell_Locations),numel(x_Range{counter}(2:end)),'lognormal');
lognormal_PDF{counter} = lognpdf(1:max(x_Range{counter}),lognormal_Parameters.mu,lognormal_Parameters.sigma);
lognormal_Error(counter) = sqrt(sum((lognormal_PDF{counter}-hist_Temp{counter}(2:end)).^2));
lognormal_Parameters_Stored(counter,1) = lognormal_Parameters.mu;
lognormal_Parameters_Stored(counter,2) = lognormal_Parameters.sigma;