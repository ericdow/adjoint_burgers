function f = godunov_flux( u )

% returns -(F_j+1/2 - F_j-1/2)

fp = zeros(size(u));
fm = zeros(size(u));

v1 = sign(u(1:end-1));
v2 = sign(u(2:end));
v3 = sign(u(1:end-1) + u(2:end));

fp(1:end-1) = -0.5*u(1:end-1).^2*0.5.*(-v1+1)*0.5.*(-v2+1) +...
    -0.5*u(2:end).^2*0.5.*(v1+1)*0.5.*(v2+1) +...
    -0.5*u(2:end).^2*0.5.*(-v1+1)*0.5.*(v2+1)*0.5.*(v3+1) +...
    -0.5*u(1:end-1).^2*0.5.*(-v1+1)*0.5.*(v2+1)*0.5.*(-v3+1);

fm(2:end) = fp(1:end-1);

fm(2) = fp(1);
fm(1) = -0.5*u(1)^2;
fp(end) = -0.5*u(end)^2;

f = -(fp - fm);
