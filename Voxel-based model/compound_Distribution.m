%% Poisson-lognormal distribution function

function mean_Dist = compound_Distribution(input_Data,nBins,xInts,sigma)

compound_Dist = zeros(nBins,numel(xInts));

zz = 1:nBins;
logValue = logninv((zz-0.5)/nBins,log(mean(input_Data))-sigma^2/2,sigma);   %Calculate value corresponding to proportions of lognormal CDF

for i = 1:nBins
    compound_Dist(i,:) = 1/nBins*poisspdf(xInts,logValue(i)); %Generate Poisson PDF with mean equal to logValue
end

mean_Dist = sum(compound_Dist)*numel(input_Data); %Average the Poisson PDFs