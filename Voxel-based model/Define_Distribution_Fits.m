%% Define variables to store distribution fits as defined in Define_IBM_Parameters

for i = 1:numel(IBM_Parameters.distribution_Fits)
    evalc([sprintf('%s_Error',IBM_Parameters.distribution_Fits(i)) ' = zeros(IBM_Parameters.t_End/100,1) ']);
    evalc([sprintf('%s_PDF',IBM_Parameters.distribution_Fits(i)) ' = cell(IBM_Parameters.t_End/100,1) ']);
    if max(strcmpi(IBM_Parameters.distribution_Fits(i),["lognormal";"normal";"negbin"]))
        evalc([sprintf('%s_Stored',IBM_Parameters.distribution_Fits(i)) ' = zeros(IBM_Parameters.t_End/100,2) ']);
    else
        evalc([sprintf('%s_Stored',IBM_Parameters.distribution_Fits(i)) ' = zeros(IBM_Parameters.t_End/100,1) ']); 
    end
end