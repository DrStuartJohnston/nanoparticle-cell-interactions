%% Top level script to run both the voxel-based model and traditional dosage
%% model in Johnston, Faria and Crampin, 2019.

close all
clear all

PDE_Model_Type = 'hybrid';

Define_PDE_Parameters;      %Define parameters needed for traditional or hybrid dosage models

Define_IBM_Parameters;      %Define parameters needed for voxel-based model

if strcmpi(PDE_Model_Type,'traditional')
    Run_TraditionalDosageModel_Solver;  %Run traditional dosage model
elseif strcmpi(PDE_Model_Type,'hybrid')
    Run_HybridModel_Solver;             %Run hybrid model
end
    
Run_IBM;                                %Run voxel-based model

Compound_Distribution_Fit;              %Fit the Poisson-lognormal distribution to the voxel-based model output

Plot_Results;                           %Plot relevant results
