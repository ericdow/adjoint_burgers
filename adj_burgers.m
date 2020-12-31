function psi = adj_burgers( psi0, u, dx, T, dt )

psi = psi0;
i = 0;
for t = dt:dt:T;
    % grab the solution
    ut = u(end-i,:);
   
    % flux
    fp = zeros(size(psi));
    fm = zeros(size(psi));
    fp(1:end-1) = 0.5*(ut(2:end).*psi(2:end)+ut(1:end-1).*psi(1:end-1)) -...
        0.5*(abs(ut(2:end)).*psi(2:end)-abs(ut(1:end-1)).*psi(1:end-1));
    fm(2:end) = fp(1:end-1);
    fm(1) = ut(1)*psi(1);
    fp(end) = ut(end)*psi(end);
    f = -(fp - fm);
    
    psi = psi + dt/dx*f;
    
    i = i+1;
end

