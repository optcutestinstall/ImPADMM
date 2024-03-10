function fiber_reconstruction_error_degrees = angleerror4(TensorODF0,beta,sizex,sizey)
% beta = 60
% check one fiber
fiber_reconstruction_error_degrees = zeros(sizex,sizey);
for ii = 1:sizex
    for jj = 1:sizey
        [v,l]=eig_dt4(TensorODF0(:,ii,jj));
        d1=acos(v(1,:)*[0; cos(beta*pi/180);sin(beta*pi/180)])*180/pi;
        d1=min(abs(d1),abs(180-d1));
        fiber_reconstruction_error_degrees(ii,jj)=d1;
    end
end
