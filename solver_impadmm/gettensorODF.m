function TensorODF = gettensorODF(GradientOrientations,S)
S0 = 1;
%% get 4th order Cartesian Tensor ODF (ave in TensorODF)
order=4;
C=constructSetOf321Polynomials(order)';
BG=constructMatrixOfIntegrals(GradientOrientations, order, 100);
P2=BG*C;
sizex = size(S,1);
sizey = size(S,2);
TensorODF = zeros(15,sizex,sizey);
for i=1:size(S,1)
    for j=1:size(S,2)
         x2 = lsqnonneg(P2,squeeze(S(i,j,1,:))/S0);
         TensorODF(:,i,j) = C*x2; 
   end
end
end
