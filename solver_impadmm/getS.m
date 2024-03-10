

function Ssolution = getS(dsolution,S0,GradientOrientations,b_value,sizex,sizey)
% dsolution: 4th-order tensor  15*sizex*sizey
l = size(GradientOrientations,1);
Gs = zeros(15,l);
Ssolution = zeros(sizex,sizey,1,l);
for i = 1:l
    Gs(:,i) = v3tov15(GradientOrientations(i,:));
end

for ir = 1:sizex
    for jc = 1:sizey
        for kk = 1:l
          %  Ssolution(ir,jc,1,kk) = S0(ir,jc)*exp(-b_value(kk)*(dsolution(:,ir,jc)'*Gs(:,kk)));
            Ssolution(ir,jc,1,kk) = S0*exp(-b_value*(dsolution(:,ir,jc)'*Gs(:,kk)));
        end
    end
end


end
