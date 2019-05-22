function xp= resampleRegularize(xp,logW)
weights = exp(logW);
weights = weights/sum(weights);

I = resample(length(weights),weights,'stratified');
xp = xp(:,I);
end