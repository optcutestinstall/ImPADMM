%-DESCRIPTION------------------------------------------------------------------------
% this Matlab script implements the ImPADMM method proposed in manuscript
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

function [out,iter] = ImPADMM(A,X,Z1,Z2,Z3,Gs,Us,Y1,Y2,Y3,Y4,Y5,Lam1,Lam2,Lam3,Lam4,Lam5,PGs,opts)

sizex = opts.sizex;
sizey = opts.sizey;
tol = opts.tol;
epsy = opts.epsy;
gam1 = opts.gam1;
gam2 = opts.gam2;
gam3 = opts.gam3;
gam4 = opts.gam4;
gam5 = opts.gam5;
alp = opts.alp;
maxiter = opts.maxiter;

pg = 1/(1 + gam4*epsy);
pu = 1/(1 + gam5*epsy);
px = alp*epsy/(epsy+4*alp);
pz1 = 1/(1 + gam1*epsy);
pz2 = 1/(1 + gam2*epsy);
pz3 = 1/(1 + gam3*epsy);
%
for iter = 1:maxiter
    %% update Zl
    Z1hat = cell(sizex,sizey); Z2hat = cell(sizex,sizey); Z3hat = cell(sizex,sizey);
    for ir = 1:sizex
        for jc = 1:sizey
            Xij = X{ir,jc};
            Z1hat{ir,jc} = pz1*(Xij - Y1{ir,jc} - epsy*(Lam1{ir,jc} - gam1*Z1{ir,jc}));
            Z2hat{ir,jc} = pz2*(Xij - Y2{ir,jc} - epsy*(Lam2{ir,jc} - gam2*Z2{ir,jc}));
            Z3hat{ir,jc} = pz3*(Xij - Y3{ir,jc} - epsy*(Lam3{ir,jc} - gam3*Z3{ir,jc}));
        end
    end
    %
    [Dset,K1,K2,K3] = generateDset(sizex,sizey);
    Z1 = updateZ1new(Z1hat,Dset,K1,gam1,epsy);
    Z2 = updateZ1new(Z2hat,Dset,K2,gam2,epsy);
    Z3 = updateZ1new(Z3hat,Dset,K3,gam3,epsy);
    %% update Gs & Us
    for ir = 1:sizex
        for jc = 1:sizey
            Xij = X{ir,jc};
            Y4ij = Y4{ir,jc};
            Lam4ij = Lam4{ir,jc};
            Gij = Gs{ir,jc};
            PGij = PGs{ir,jc};
            Uij = Us{ir,jc};
            Lam5ij = Lam5{ir,jc};
            Y5ij = Y5{ir,jc};
            temp1 = Xij - Y4ij - epsy*Lam4ij - PGij;
            temp2 = TtoM(temp1);
            temp3 = (temp2 + Uij + Y5ij + epsy*(Lam5ij + gam4*Gij))*pg;
            temp3 = (temp3+temp3')/2;
            [ug,dg] = eig(temp3);
            ndg = max(diag(dg),0);
            Gijnew = ug*(diag(ndg)*ug');
            Gs{ir,jc} = Gijnew;
            PGs{ir,jc} = GmtoT(Gijnew);
            temp4 = Gijnew - Y5ij - epsy*(Lam5ij - gam5*Uij);
            temp4 = (temp4 + temp4')/2;
            [uu,du,vu] = svd(temp4);
            du = diag(du);
            temp5 = uu*(diag([du(1:3);zeros(3,1)])*vu');
            Us{ir,jc} = temp5*pu;
        end
    end   
    %% update X
    for ir = 1:sizex
        for jc = 1:sizey
            X{ir,jc} = px*((1/alp)*A{ir,jc} + Lam1{ir,jc} + Lam2{ir,jc} + Lam3{ir,jc} + Lam4{ir,jc} + (1/epsy)*(Y1{ir,jc} + Y2{ir,jc} + Y3{ir,jc} + Y4{ir,jc} + PGs{ir,jc} + Z1{ir,jc} + Z2{ir,jc} + Z3{ir,jc}));
        end
    end
    %% update Yl
    maxys = 0;
    for ir = 1:sizex
        for jc = 1:sizey
             Xij = X{ir,jc};
             Y1ij = 2/3*(Xij - Z1{ir,jc} - epsy*Lam1{ir,jc});
             Y1{ir,jc} = Y1ij;
             Y2ij = 2/3*(Xij - Z2{ir,jc} - epsy*Lam2{ir,jc});
             Y2{ir,jc} = Y2ij;
             Y3ij = 2/3*(Xij - Z3{ir,jc} - epsy*Lam3{ir,jc});
             Y3{ir,jc} = Y3ij;
             Y4ij = 2/3*(Xij - PGs{ir,jc} - epsy*Lam4{ir,jc});
             Y4{ir,jc} = Y4ij;
            Y5ij = 2/3*(Gs{ir,jc} - Us{ir,jc} - epsy*Lam5{ir,jc});
            Y5{ir,jc} = Y5ij;
            maxys = max(maxys,max(abs(Y5ij(:))));
        end
    end
    %% test stop
    if maxys < tol 
        break;
    end
    %% update Lams
    for ir = 1:sizex
        for jc = 1:sizey
            Xij = X{ir,jc};
            Lam1{ir,jc} = Lam1{ir,jc} + (1/epsy)*(Z1{ir,jc} - Xij + Y1{ir,jc});
            Lam2{ir,jc} = Lam2{ir,jc} + (1/epsy)*(Z2{ir,jc} - Xij + Y2{ir,jc});
            Lam3{ir,jc} = Lam3{ir,jc} + (1/epsy)*(Z3{ir,jc} - Xij + Y3{ir,jc});
            Lam4{ir,jc} = Lam4{ir,jc} + (1/epsy)*(PGs{ir,jc} - Xij + Y4{ir,jc});
            Lam5{ir,jc} = Lam5{ir,jc} + (1/epsy)*(Us{ir,jc} - Gs{ir,jc} + Y5{ir,jc});
        end
    end
    
end
out.X = X;
out.Z1 = Z1;
out.Z2 = Z2;
out.Z3 = Z3;
out.Gs = Gs;
out.Us = Us;
out.Y1 = Y1;
out.Y2 = Y2;
out.Y3 = Y3;
out.Y4 = Y4;
out.Y5 = Y5;
end