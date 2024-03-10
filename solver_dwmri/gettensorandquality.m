function [UniqueTensorCoefficients,TensorODF,md,GA] = gettensorandquality(GradientOrientations,S,b_value)
S0 = 1;
%% get 4th order Cartesian Tensor ODF (ave in TensorODF)
order=4;
G=constructMatrixOfMonomials(GradientOrientations, order);
C=constructSetOf321Polynomials(order)';
BG=constructMatrixOfIntegrals(GradientOrientations, order, 100);
P2=BG*C;
P=G*C;
P=-diag(b_value)*P;
sizex = size(S,1);
sizey = size(S,2);
UniqueTensorCoefficients=zeros(15,sizex,sizey);
TensorODF = zeros(15,sizex,sizey);
for i=1:size(S,1)
    for j=1:size(S,2)
         y=log(squeeze(S(i,j,1,:))/S0);
         x=lsqnonneg(P, y);   % non-negative least squares
         UniqueTensorCoefficients(:,i,j) = C * x;
         x2 = lsqnonneg(P2,squeeze(S(i,j,1,:))/S0);
         TensorODF(:,i,j) = C*x2; 
   end
end

md = size(sizex,sizey);
GA = size(sizex,sizey);

% the mean diffusivity (only at voxel(1,1))
for ir = 1:sizex
    for jc = 1:sizey
        t=UniqueTensorCoefficients(:,ir,jc);
        md(ir,jc)=(t(1)+t(5)+t(15)+(t(3)+t(10)+t(12))/3)/5;
        % the generalized anisotropy GA (only at voxel(1,1))
        var=GeneralizedVariance(UniqueTensorCoefficients(:,ir,jc));
        GA(ir,jc)=1-1/(1+(250*var)^(1+1/(1+5000*var)));
    end
end


