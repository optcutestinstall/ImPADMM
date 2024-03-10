function zbar =  bcgpai2(Z10hat,gam1,epsy)
% Z10hat, 81*3 matrix
% U1, U2: 81*1 matrices
gam = 0.55;
imax = 1000;
[rzbar,czbar] = size(Z10hat);
U1 = zeros(rzbar,1);
U2 = zeros(rzbar,1);
l = epsy/(1 + gam1*epsy);
l2 = 2*l;
c = gam*l2;

zhat1 = Z10hat(:,1);
zhat2 = Z10hat(:,2);
zhat3 = Z10hat(:,3);
g1 = zhat1 - zhat2;
g2 = zhat1 - zhat3;

for iter = 1:imax
    U1old = U1;
    U2old = U2;
    U1t = U1 - (l*(2*U1 + U2) + g1)/c;
    U1 = U1t;
    np1 = norm(U1,2);
    if np1 > 1
        U1 = U1t/np1;
    end
    
    U2t = U2 - (l*(U1 + 2*U2) + g2)/c;
    U2 = U2t;
    np2 = norm(U2,2);
    if np2 > 1
        U2 = U2t/np2;
    end
    if norm(U1 - U1old,2) < 1e-3 && norm(U2 - U2old,2) < 1e-3
        break;
    end
end
zbar = Z10hat + l*[U1 + U2, -U1, -U2];
end