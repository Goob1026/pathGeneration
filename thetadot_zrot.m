function thetadot = thetadot_zrot(thetax,thetay,thetaz,k)

thetazdot = 0*thetaz;

thetaydot = -k*sin(thetay/2);
thetaxdot = -k*(sin(thetax/2)*cos(thetay)^2 - sin(thetax)^2*sin(thetay/2));


thetadot = [thetaxdot;thetaydot;thetazdot];
end