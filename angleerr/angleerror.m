function fiberr = angleerror(TensorODF0,alp,beta,sizex,sizey)
fiber_reconstruction_error_degrees = zeros(sizex,sizey);
% check 2 cross fibers
for ii = 1:sizex
    for jj = 1:sizey
        [v,l]=eig_dt4(TensorODF0(:,ii,jj));
        if isempty(v)
            fiber_reconstruction_error_degrees(ii,jj)=180;
        else
        d1=acos(v(1,:)*[cos(alp*pi/180);sin(alp*pi/180);0])*180/pi;
        d1=min(abs(d1),abs(180-d1));
        d2=acos(v(2,:)*[cos(beta*pi/180);sin(beta*pi/180);0])*180/pi;
        d2=min(abs(d2),abs(180-d2));
        error1=(d1+d2)/2;
        d1=acos(v(2,:)*[cos(alp*pi/180);sin(alp*pi/180);0])*180/pi;
        d1=min(abs(d1),abs(180-d1));
        d2=acos(v(1,:)*[cos(beta*pi/180);sin(beta*pi/180);0])*180/pi;
        d2=min(abs(d2),abs(180-d2));
        error2=(d1+d2)/2;
        fiber_reconstruction_error_degrees(ii,jj)=min(error1,error2);
        end
    end
end
fiberr = mean(fiber_reconstruction_error_degrees(:));
end
