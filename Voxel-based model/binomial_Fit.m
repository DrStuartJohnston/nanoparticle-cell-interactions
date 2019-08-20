[binomial_Parameters] = fitdist(NPs_per_Cell(cell_Locations),'binomial','NTrials',round(t));
binomial_PDF{counter} = binopdf(x_Range{counter},t,binomial_Parameters.p);
binomial_Error(counter) = sqrt(sum((binomial_PDF{counter}-hist_Temp{counter}).^2));
binomial_Parameters_Stored(counter) = binomial_Parameters.p;