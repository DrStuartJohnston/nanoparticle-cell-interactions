%% Fit the Poisson-lognormal distribution to the number of nanoparticles per cell

nBins = 1000;       %Number of terms in integral w.r.t. lognormal component

xInts = x_Range{end};   %Range of nanoparticles per cell

hist_Data = hist(NPs_per_Cell(cell_Locations),0:max(xInts));    %Histogram of nanoparticles per cell

sigma_Best = lsqnonlin(@(sigma) compound_Distribution(NPs_per_Cell(cell_Locations),nBins,0:max(xInts),sigma)-hist_Data,0.6,0,10,opts); %Calculate apparent heterogeneity via least squares

%% Plot fit between histogram and best-fit Poisson-lognormal distribution
figure; hist(NPs_per_Cell(cell_Locations),0:max(xInts))
hold on
plot(xInts,compound_Distribution(NPs_per_Cell(cell_Locations),nBins,xInts,sigma_Best))
