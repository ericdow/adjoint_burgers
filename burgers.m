function u = burgers(u0, dx, T, dt)
% Solve forward problem
% u_t - u u_x = 0, -L < x < 2L, 0 < t < T

u = zeros(length(0:dt:T),length(u0));
u(1,:) = u0;
i = 2;
fp = zeros(size(u0));
fm = zeros(size(u0));
for t = dt:dt:T
    v1 = sign(u(i-1,1:end-1));
    v2 = sign(u(i-1,2:end));
    v3 = sign(u(i-1,1:end-1) + u(i-1,2:end));
    
    fp(1:end-1) = -0.5*u(i-1,1:end-1).^2*0.5.*(-v1+1)*0.5.*(-v2+1) +...
        -0.5*u(i-1,2:end).^2*0.5.*(v1+1)*0.5.*(v2+1) +...
        -0.5*u(i-1,2:end).^2*0.5.*(-v1+1)*0.5.*(v2+1)*0.5.*(v3+1) +...
        -0.5*u(i-1,1:end-1).^2*0.5.*(-v1+1)*0.5.*(v2+1)*0.5.*(-v3+1);
    
    fm(2:end) = fp(1:end-1);
    
    fm(2) = fp(1);
    fm(1) = -0.5*u(i-1,1)^2;
    fp(end) = -0.5*u(i-1,end)^2;
    
    f = -(fp - fm);
    u(i,:) = u(i-1,:) + dt/dx*f;
    
    i = i+1;
end

end

