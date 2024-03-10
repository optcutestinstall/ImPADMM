%-DESCRIPTION------------------------------------------------------------------------
% this Matlab script provides the demo for the second test in manuscript
%
% "An inexact majorized proximal proximal alternating direction methods 
%  of multipliers for diffusion tensor constraints regularization"
%  
%-AUTHOR-----------------------------------------------------------------------------
% Hong Zhu   
% School of Mathematical Sciences, Jiangsu University, 
% Xuefu Road,Zhenjiang, 212013, Jiangsu, China 
% zhuhongmath@126.com
% Michael K. NG
% Department of Mathematices, HongKong Baptist University, 
% Kowloon Tong, Hong Kong, 999077, China 
% michael-ng@hkbu.edu.hk

%% generate DW-MRI data
% 1 voxel only 
% 21 gradient directions
clc;clear;

addpath(genpath(pwd));
casepha = 21;   % number of gradients
rng(2000)
sizex = 6;    
sizey = 6;     % totle number of voxels: sizex - by - sizey  (5*5)
S0 = 1;
% generate casepha*3 matrix GradientOrientations (i.e., 21 3-dimensional directions)
GradientOrientations=[0.1639 0.5115 0.8435; 0.1176 -0.5388 0.8342; 0.5554 0.8278 -0.0797; -0.4804 0.8719 0.0948; 0.9251 -0.0442 0.3772; 0.7512 -0.0273 -0.6596; 0.1655 -0.0161 0.9861; 0.6129 -0.3427 0.7120; 0.6401 0.2747 0.7175; -0.3724 -0.3007 0.8780; -0.3451 0.3167 0.8835; 0.4228 0.7872 0.4489; 0.0441 0.9990 0.0089; -0.1860 0.8131 0.5515; 0.8702 0.4606 0.1748; -0.7239 0.5285 0.4434; -0.2574 -0.8032 0.5372; 0.3515 -0.8292 0.4346; -0.7680 -0.4705 0.4346; 0.8261 -0.5384 0.1660; 0.9852 -0.0420 -0.1660];
b_value = 1500;
% angles of fibers
fiber_orientation1=[cos(20*pi/180) sin(20*pi/180) 0];
fiber_orientation2=[0 cos(60*pi/180) sin(60*pi/180)];
fiber_orientation3=[cos(100*pi/180) sin(100*pi/180) 0];

% generate S
S=zeros(sizex,sizey,1,size(GradientOrientations,1),1);

for i=1:size(GradientOrientations,1)
    for x = 1:sizex/2
        for y = 1:sizey/2
            S(x,y,1,i)=S0* (SimulateDWMRI(fiber_orientation1,GradientOrientations(i,:))+ SimulateDWMRI(fiber_orientation2,GradientOrientations(i,:)))/2;
        end
    end
end
for i=1:size(GradientOrientations,1)
    for x = 1:sizex/2
        for y = sizey/2+1:sizey
            S(x,y,1,i)=S0* (SimulateDWMRI(fiber_orientation1,GradientOrientations(i,:))+ SimulateDWMRI(fiber_orientation3,GradientOrientations(i,:)))/2;
        end
    end
end

for i=1:size(GradientOrientations,1)
    for x = sizex/2+1:sizex
        for y = 1:sizey/2
            S(x,y,1,i)=S0*SimulateDWMRI(fiber_orientation2,GradientOrientations(i,:));
        end
    end
end

for i=1:size(GradientOrientations,1)
    for x = sizex/2+1:sizex
        for y = sizey/2+1:sizey
            S(x,y,1,i)=S0* (SimulateDWMRI(fiber_orientation2,GradientOrientations(i,:))+ SimulateDWMRI(fiber_orientation3,GradientOrientations(i,:)) + SimulateDWMRI(fiber_orientation1,GradientOrientations(i,:)))/3;
        end
    end
end
%% --- Estimation of Diffusion Tensors using MATLAB (DTI)
% GA0:  generalized anisotropy GA
[UniqueTensorCoefficients0,TensorODF0,md0,GA0] = gettensorandquality(GradientOrientations,S,b_value);
figure;
plotTensors(UniqueTensorCoefficients0,1,[321 1]); title('Original tensor field');
figure;
plotTensors(TensorODF0,1,[321 1]); title('ODF: Original tensor field');

%%

ae0 = size(6,6);
alp = 20; beta =60;
ae0(1:3,1:3) = angleerror3(TensorODF0(:,1:3,1:3),alp,beta,3,3);
%
alp = 20; beta =100;
ae0(1:3,4:6) = angleerror(TensorODF0(:,1:3,4:6),alp,beta,3,3);
%
beta =60;
ae0(4:6,1:3) = angleerror4(TensorODF0(:,4:6,1:3),beta,3,3);
%
%
alp = 20; beta =60; gamma =100; 
ae0(4:6,4:6) = angleerror2(TensorODF0(:,4:6,4:6),alp,beta,gamma,3,3);


%% --- add gaussian noise

sigma = 0.05;
UniqueTensorCoefficients = UniqueTensorCoefficients0;
sindex2 = [];  
for ir = 1:sizey
    sindex2 = [sindex2;[(1:sizex)',ir*ones(sizex,1)]];
end
for iii = 1:size(sindex2,1)
     k1 = sindex2(iii,1);k2 = sindex2(iii,2);
     for ii = 1:size(GradientOrientations,1)
         S(k1,k2,1,ii) = sqrt((S(k1,k2,1,ii)+sigma*randn(1))^2+(sigma*randn(1))^2);
     end
end



%% test dissusion tensor MRI using Matlab
tic;
[HODTI,TensorODF,md_matlab,ga_matlab]= gettensorandquality(GradientOrientations,S,b_value);
tmatlab = toc;
figure;
plotTensors(HODTI,1,[321 1]); title('HODT by NNLS');

figure;
plotTensors(TensorODF,1,[321 1]); title('ODF by NNLS');


aennls2 = size(sizex,sizey);

alp = 20; beta =60;  %  
aennls2(1:3,1:3) = angleerror3(TensorODF(:,1:3,1:3),alp,beta,3,3);
%
alp = 20; beta =100;
aennls2(1:3,4:6) = angleerror(TensorODF(:,1:3,4:6),alp,beta,3,3);
%
 beta =60;
aennls2(4:6,1:3) = angleerror4(TensorODF(:,4:6,1:3),beta,3,3);
%
%
alp = 20; beta =60; gamma =100; 
aennls2(4:6,4:6) = angleerror2(TensorODF(:,4:6,4:6),alp,beta,gamma,3,3);
%

%% get A
cAuc = getA(S,GradientOrientations,b_value,sizex,sizey);

figure;
plotTensors(cAuc,1,[321 1]); title('noisy');

%% test ImPADMM
opts = [];
opts.sizex = sizex;
opts.sizey = sizey;
opts.tol = 1e-5;
opts.epsy = 1e-4;
opts.maxiter = 200000;
opts.gam1 = .1;
opts.gam2 = .1;
opts.gam3 = .1;
opts.gam4 = .1;
opts.gam5 = .1;
opts.alp = 1e-5;  % alp should be small

%% initial point
% generate Gs sizex*sizey cell, each entry is a symmetric sdp matrix
% generate Us sizex*sizey cell, each entry is a matrix with rank lack of 3
% generate X: sizex*sizey cell, each entry is a 4th order tensor
% generate Zl: sizex*sizey cell, each entry is a 4th order tensor
Gs = cell(sizex,sizey);
PGs = cell(sizex,sizey);
Us = cell(sizex,sizey);
X = cell(sizex,sizey);
Y1 = cell(sizex,sizey);
Y4 = cell(sizex,sizey);
Y5 = cell(sizex,sizey);
Lam1 = cell(sizex,sizey);
Lam5 = cell(sizex,sizey);
A = cell(sizex,sizey);
for ir = 1:sizex
    for ic = 1:sizey
        temp = max(cAuc(:,ir,ic),0);
        tempx = uctoT(temp);
        A{ir,ic} = uctoT(UniqueTensorCoefficients(:,ir,ic));   
        gij = TtoM(tempx);
        [ug,dg] = eig(gij);
        dg = diag(dg);
        ndg = max(dg,0);
        matg = ug*(diag(ndg)*ug');
        Gs{ir,ic} = matg;
        ndgr3 = [ndg(1:3);zeros(3,1)];
        matr3 = ug*(diag(ndgr3)*ug');
        Us{ir,ic} = matr3;
        Xij = GmtoT(matr3); 
        X{ir,ic} = Xij;
        Y1{ir,ic} = zeros(3,3,3,3);
        tmatg = GmtoT(matg);
        PGs{ir,ic} = tmatg;
        Y4{ir,ic} = Xij - tmatg;
        Y5{ir,ic} = matg - matr3;
        Lam1{ir,ic} = zeros(3,3,3,3);
        Lam5{ir,ic} = zeros(6,6);
    end
end
Z1 = X; Z2 = X; Z3 = X;
Y2 = Y1; Y3 = Y1;
Lam2 = Lam1; Lam3 = Lam1; Lam4 = Lam1; 

tic;
[out,iter] = ImPADMM(A,X,Z1,Z2,Z3,Gs,Us,Y1,Y2,Y3,Y4,Y5,Lam1,Lam2,Lam3,Lam4,Lam5,PGs,opts);
tdtr = toc;

[md_DTR,ga_DTR,UniqueUr] = figuretestnew(out,sizex,sizey);
figure;
plotTensors(UniqueUr,1,[321 1]); title('HODT by ImPADMM');



SDTR = getS(UniqueUr,S0,GradientOrientations,b_value,sizex,sizey);
TODFDTR = gettensorODF(GradientOrientations,SDTR);
figure;
plotTensors(TODFDTR,1,[321 1]); title('ODF by ImPADMM');


aedtr2 = size(sizex,sizey);

alp = 20;  beta =60;
aedtr2(1:3,1:3) = angleerror3(TODFDTR(:,1:3,1:3),alp,beta,3,3);
%
alp = 20; beta =100;
aedtr2(1:3,4:6) = angleerror(TODFDTR(:,1:3,4:6),alp,beta,3,3);
%
 beta =60;
aedtr2(4:6,1:3) = angleerror4(TODFDTR(:,4:6,1:3),beta,3,3);
%
%
alp = 20; beta =60; gamma =100; 
aedtr2(4:6,4:6) = angleerror2(TODFDTR(:,4:6,4:6),alp,beta,gamma,3,3);
 
[mean(mean(aennls2(1:3,4:6))),mean(mean(aennls2(1:3,1:3))),mean(mean(aennls2(4:6,1:3))),mean(mean(aennls2(4:6,4:6)))]
[mean(mean(aedtr2(1:3,4:6))),mean(mean(aedtr2(1:3,1:3))),mean(mean(aedtr2(4:6,1:3))),mean(mean(aedtr2(4:6,4:6)))]
%%