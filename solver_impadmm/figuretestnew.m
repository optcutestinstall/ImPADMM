% firuges and average fitting errors for results obtained from DTR
%function [md_DTR,ga_DTR,avfe_dtr,minfe_dtr,maxfe_dtr,fe] = figuretestnew(out,UniqueTensorCoefficients0,GradientOrientations,sizex,sizey)
%function [md_DTR,ga_DTR,avfe_dtr,stdfe_dtr,UniqueUr] = figuretestnew(out,UniqueTensorCoefficients0,GradientOrientations,sizex,sizey)
function [md_DTR,ga_DTR,UniqueUr] = figuretestnew(out,sizex,sizey)
Usout = out.Us;
Xr = cell(sizex,sizey);
UniqueUr = zeros(15,sizex,sizey);
for ir = 1:sizex
    for jc = 1:sizey
        temp = Usout{ir,jc};
        temp = (temp + temp')/2;
        [ut,dt] = eig(temp);
        dt = diag(dt);
        dd = [zeros(3,1);max(0,dt(4));dt(5:6)];
        temp2 = ut*(diag(dd)*ut');
        T = GmtoT(temp2);
        Xr{ir,jc} = T;
        UniqueUr(:,ir,jc) = Ttouc(T);
    end
end
[md_DTR,ga_DTR] = mdga(UniqueUr,sizex,sizey);
%[avfe_dtr,minfe_dtr,maxfe_dtr,fe] = fittingerror(UniqueTensorCoefficients0,GradientOrientations,UniqueUr,sizex,sizey);
%[avfe_dtr,stdfe_dtr] = fittingerror(UniqueTensorCoefficients0,GradientOrientations,UniqueUr,sizex,sizey);
% figure;
% plotTensors(UniqueUr,1,[321 1]); title('DTR');
%  
end