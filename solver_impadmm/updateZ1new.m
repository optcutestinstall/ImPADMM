% K1 = {0}; D_0 = {(1,1),(2,2),(3,3),...,(sizex,sizex)}
function Z1 = updateZ1new(Z1hat,Dset,K1,gam1,epsy)
Z1 = Z1hat;
sizex = size(Z1,1);
sizey = size(Z1,2);
% 
ln = length(K1);
for jj = 1:ln
    dindex = K1(jj);
    if dindex == 0
        dind = 1;
    elseif dindex < 0
        dind = 2*abs(dindex);
    else
        dind = 2*dindex+1;
    end
    ijdset = Dset{dind,2};
    lnij = size(ijdset,1);
    for kk = 1:lnij
        i = ijdset(kk,1); j = ijdset(kk,2);
        z10hat11 = Z1hat{i,j}(:);
        z10hat12 = Z1hat{min(i,sizex),min(j+1,sizey)}(:);
        z10hat21 = Z1hat{min(i+1,sizex),min(j,sizey)}(:);
        Z10hat = [z10hat11,z10hat12,z10hat21];
        z10test1 = z10hat11 - 2*z10hat12 + z10hat21;
        z10test2 = z10hat11 + z10hat12 - 2*z10hat21;
        if norm(z10test1,2) <= (3*epsy)/(1 + gam1*epsy) && norm(z10test2,2) <= (3*epsy)/(1 + gam1*epsy)
            zt1 = z10hat11 + z10hat12 + z10hat21;
            zt2 = 3*z10hat11 + z10hat12 - z10hat21;
            zbar = [zt1,zt2,zt1]/3;
        else
        zbar =  bcgpai2(Z10hat,gam1,epsy);
        end
        Z1{i,j} = reshape(zbar(:,1),3,3,3,3);
        if j+1 <= sizey
           Z1{i,j+1} = reshape(zbar(:,2),3,3,3,3);
        end
        if i+1 <= sizex
           Z1{i+1,j} = reshape(zbar(:,3),3,3,3,3); 
        end
    end
end
end