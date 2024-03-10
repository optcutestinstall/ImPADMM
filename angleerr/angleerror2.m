function fiber_reconstruction_error_degrees = angleerror2(TensorODF0,alp,beta,gamma,sizex,sizey)
% beta = 60
% check 3 cross fibers
fiber_reconstruction_error_degrees = zeros(sizex,sizey);
for ii = 1:sizex
    for jj = 1:sizey
        [v,l]=eig_dt4(TensorODF0(:,ii,jj));
        % 1, 2, 3
        d1=acos(v(1,:)*[cos(alp*pi/180);sin(alp*pi/180);0])*180/pi;
        d1=min(abs(d1),abs(180-d1));
        d2=acos(v(2,:)*[0; cos(beta*pi/180);sin(beta*pi/180)])*180/pi;
        d2=min(abs(d2),abs(180-d2));
        d3=acos(v(3,:)*[cos(gamma*pi/180);sin(gamma*pi/180);0])*180/pi;
        d3=min(abs(d3),abs(180-d3));
        error1=(d1+d2+d3)/3;
        % 2, 1, 3,
        d1=acos(v(2,:)*[cos(alp*pi/180);sin(alp*pi/180);0])*180/pi;
        d1=min(abs(d1),abs(180-d1));
        d2=acos(v(1,:)*[0; cos(beta*pi/180);sin(beta*pi/180)])*180/pi;
        d2=min(abs(d2),abs(180-d2));
        d3=acos(v(3,:)*[cos(gamma*pi/180);sin(gamma*pi/180);0])*180/pi;
        d3=min(abs(d3),abs(180-d3));
        error2=(d1+d2+d3)/3;
        % 1 3,2
        d1=acos(v(1,:)*[cos(alp*pi/180);sin(alp*pi/180);0])*180/pi;
        d1=min(abs(d1),abs(180-d1));
        d2=acos(v(3,:)*[0;cos(beta*pi/180);sin(beta*pi/180)])*180/pi;
        d2=min(abs(d2),abs(180-d2));
        d3=acos(v(2,:)*[cos(gamma*pi/180);sin(gamma*pi/180);0])*180/pi;
        d3=min(abs(d3),abs(180-d3));
        error3=(d1+d2+d3)/3;
        % 3, 1, 2
        d1=acos(v(3,:)*[cos(alp*pi/180);sin(alp*pi/180);0])*180/pi;
        d1=min(abs(d1),abs(180-d1));
        d2=acos(v(1,:)*[0;cos(beta*pi/180);sin(beta*pi/180)])*180/pi;
        d2=min(abs(d2),abs(180-d2));
        d3=acos(v(2,:)*[cos(gamma*pi/180);sin(gamma*pi/180);0])*180/pi;
        d3=min(abs(d3),abs(180-d3));
        error4=(d1+d2+d3)/3;
        % 2 3 1
        d1=acos(v(2,:)*[cos(alp*pi/180);sin(alp*pi/180);0])*180/pi;
        d1=min(abs(d1),abs(180-d1));
        d2=acos(v(3,:)*[0;cos(beta*pi/180);sin(beta*pi/180)])*180/pi;
        d2=min(abs(d2),abs(180-d2));
        d3=acos(v(1,:)*[cos(gamma*pi/180);sin(gamma*pi/180);0])*180/pi;
        d3=min(abs(d3),abs(180-d3));
        error5=(d1+d2+d3)/3;
        % 3 2 1
        d1=acos(v(3,:)*[cos(alp*pi/180);sin(alp*pi/180);0])*180/pi;
        d1=min(abs(d1),abs(180-d1));
        d2=acos(v(2,:)*[0;cos(beta*pi/180);sin(beta*pi/180)])*180/pi;
        d2=min(abs(d2),abs(180-d2));
        d3=acos(v(1,:)*[cos(gamma*pi/180);sin(gamma*pi/180);0])*180/pi;
        d3=min(abs(d3),abs(180-d3));
        error6=(d1+d2+d3)/3;
        %
        error = [error1;error2;error3;error4;error5;error6];
        fiber_reconstruction_error_degrees(ii,jj)=min(error(:));
    end
end
