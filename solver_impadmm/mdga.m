function [md,ga] = mdga(UniqueTensorCoefficients,sizex,sizey)
md = zeros(sizex,sizey);
ga = zeros(sizex,sizey);
for ir = 1:sizex
    for jc = 1:sizey
        t=UniqueTensorCoefficients(:,ir,jc);
        md(ir,jc)=(t(1)+t(5)+t(15)+(t(3)+t(10)+t(12))/3)/5;
        var=GeneralizedVariance(UniqueTensorCoefficients(:,ir,jc));
        ga(ir,jc)=1-1/(1+(250*var)^(1+1/(1+5000*var)));
    end
end
end