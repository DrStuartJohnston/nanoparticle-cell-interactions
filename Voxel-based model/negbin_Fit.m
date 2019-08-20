[negbin_Parameters] = fitdist(NPs_per_Cell(cell_Locations),'nbin');
negbin_PDF{counter} = nbinpdf(x_Range{counter},negbin_Parameters.R,negbin_Parameters.P);
negbin_Error(counter) = sqrt(sum((negbin_PDF{counter}-hist_Temp{counter}).^2));
negbin_Parameters_Stored(counter,1) = negbin_Parameters.R;
negbin_Parameters_Stored(counter,2) = negbin_Parameters.P;