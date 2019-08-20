[poisson_Parameters] = fitdist(NPs_per_Cell(cell_Locations),'poisson');
poisson_PDF{counter} = poisspdf(x_Range{counter},poisson_Parameters.lambda);
poisson_Error(counter) = sqrt(sum((poisson_PDF{counter}-hist_Temp{counter}).^2));
poisson_Parameters_Stored(counter) = poisson_Parameters.lambda;