%% Plot relevant results for the dosing framework

%% Plot concentration profile
figure(1)
plot(linspace(0,PDE_Parameters.width,numel(x)),scaled_Final_Profile,'linewidth',2)
hold on
plot(linspace(0,PDE_Parameters.width,IBM_Parameters.grid_Height),vertical_Slice/IBM_Parameters.NP_Initial_Density,'linewidth',2)

%% Plot uptake curve
figure(2)
plot(time/max(time)*PDE_Parameters.final_Time/3600,uptake_Curve,'linewidth',2)
hold on
plot([100:100:IBM_Parameters.t_End]/3600,sum(NPs_per_Cell_Evolution(:,cell_Locations)./cell_Evolution,2),'linewidth',2)

%% Plot particle per cell histogram
figure(3)
hist(NPs_per_Cell(cell_Locations),x_Range{counter})

%% Plot error with each distribution fit
figure(4)
hold on
for i = 1:numel(IBM_Parameters.distribution_Fits)
    plot_Data = eval([sprintf('%s_Error',IBM_Parameters.distribution_Fits(i))]);
    plot([100:100:IBM_Parameters.t_End]/3600,plot_Data,'linewidth',2)
end

%% Plot compound distribution fit
figure(5)
hist(NPs_per_Cell(cell_Locations),0:max(xInts))
hold on
plot(xInts,compound_Distribution(NPs_per_Cell(cell_Locations),nBins,xInts,sigma_Best),'linewidth',2)

%% Plot each distribution fit

figure(6)
hist(NPs_per_Cell(cell_Locations),x_Range{end})
hold on
for i = 1:numel(IBM_Parameters.distribution_Fits)
    plot_Data = eval([sprintf('%s_PDF{end}',IBM_Parameters.distribution_Fits(i))]);
    if strcmpi(IBM_Parameters.distribution_Fits(i),'lognormal')
        plot(1:max(x_Range{end}),plot_Data*numel(cell_Locations),'linewidth',2)
    else
        plot(x_Range{end},plot_Data*numel(cell_Locations),'linewidth',2)
    end
end

figure(6)
fit_Times = time_Fit_Vec/100;
for i = 1:numel(fit_Times)
   figure(6+i)
   hist(NPs_per_Cell_Evolution(fit_Times(i),cell_Locations_Stored{i}),x_Range{fit_Times(i)});
   hold on
   plot(x_Range{fit_Times(i)},compound_Distribution(NPs_per_Cell_Evolution(fit_Times(i),cell_Locations_Stored{i}),nBins,x_Range{fit_Times(i)},sigma_Best_Vec(i)))
   plot(x_Range{fit_Times(i)},compound_Distribution(NPs_per_Cell_Evolution(fit_Times(i),cell_Locations_Stored{i}),nBins,x_Range{fit_Times(i)},IBM_Parameters.affinity_StDev))
end

figure;
plot(time_Fit_Vec/3600,sigma_Best_Vec)