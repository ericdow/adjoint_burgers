from numpy import *
import scipy.sparse
import matplotlib.pyplot as plt

def burgers(u0, dx, t):
    
    # solve burgers equation using upwind FVM
    dt = t[1] - t[0]
    u = zeros((len(t),len(u0)))
    i = 0
    fp = zeros(len(u0))
    fm = zeros(len(u0))
    for ti in t:
        if i == 0:
            u[i,:] = u0
        else:
            v1 = sign(u[i-1,:-1])
            v2 = sign(u[i-1,1:])
            v3 = sign(u[i-1,:-1] + u[i-1,1:])
            
            fp[:-1] = -0.5*u[i-1,:-1]**2*0.5*(-v1+1)*0.5*(-v2+1) + \
                      -0.5*u[i-1,:-1]**2*0.5*(-v1+1)*0.5*(v2+1)*0.5*(-v3+1) + \
                      -0.5*u[i-1,1:]**2*0.5*(v1+1)*0.5*(v2+1) + \
                      -0.5*u[i-1,1:]**2*0.5*(-v1+1)*0.5*(v2+1)*0.5*(v3+1)
            
            fm[1:] = fp[:-1]
            
            fm[0]  = -0.5*u[i-1,0]**2
            fp[-1] = 0.0
            
            f = -(fp - fm);
            u[i,:] = u[i-1,:] + dt/dx*f;
        
        i += 1;

    return u

def burgers_complex(u0, dx, t):
    
    # solve burgers equation using upwind FVM
    dt = t[1] - t[0]
    u = zeros((len(t),len(u0)),dtype=complex)
    i = 0
    fp = zeros(len(u0),dtype=complex)
    fm = zeros(len(u0),dtype=complex)
    for ti in t:
        if i == 0:
            u[i,:] = u0
        else:
            v1 = sign(u[i-1,:-1])
            v2 = sign(u[i-1,1:])
            v3 = sign(u[i-1,:-1] + u[i-1,1:])
            
            fp[:-1] = -0.5*u[i-1,:-1]**2*0.5*(-v1+1)*0.5*(-v2+1) + \
                      -0.5*u[i-1,:-1]**2*0.5*(-v1+1)*0.5*(v2+1)*0.5*(-v3+1) + \
                      -0.5*u[i-1,1:]**2*0.5*(v1+1)*0.5*(v2+1) + \
                      -0.5*u[i-1,1:]**2*0.5*(-v1+1)*0.5*(v2+1)*0.5*(v3+1)
            
            fm[1:] = fp[:-1]
            
            fm[0]  = -0.5*u[i-1,0]**2
            fp[-1] = 0.0
            
            f = -(fp - fm);
            u[i,:] = u[i-1,:] + dt/dx*f;
        
        i += 1;

    return u

def flux(u):
    fp = zeros(len(u),dtype=complex)
    fm = zeros(len(u),dtype=complex)

    v1 = sign(real(u[:-1]))
    v2 = sign(real(u[1:]))
    v3 = sign(real(u[:-1] + u[1:]))
    
    fp[:-1] = -0.5*u[:-1]**2*0.5*(-v1+1)*0.5*(-v2+1) + \
              -0.5*u[:-1]**2*0.5*(-v1+1)*0.5*(v2+1)*0.5*(-v3+1) + \
              -0.5*u[1:]**2*0.5*(v1+1)*0.5*(v2+1) + \
              -0.5*u[1:]**2*0.5*(-v1+1)*0.5*(v2+1)*0.5*(v3+1)
    
    fm[1:] = fp[:-1]
    
    fm[0]  = -0.5*u[0]**2
    fp[-1] = 0.0
    
    f = -(fp - fm);

    return f

def adj_burgers(u, dx, t):
   
    # solve the discrete adjoint problem
    N = len(u[0,:]) 
    dt = t[1] - t[0]
    v = zeros((len(t),N))
    i = 0
    for ti in t:

        # check dF_du with complex step
        '''
        h = 1.0e-20
        dF_du_fd = zeros((N,N))
        for j in range(N):
            uh = zeros(N,dtype=complex)
            uh += ut
            uh[j] += complex(0,h)
            fh = flux(uh)/dx
            dF_du_fd[:,j] = imag(fh)/h

        print dF_du_fd
        print dF_du.todense()
        print " "
        '''

        A2 =  eye(N)/dt

        if i == 0:
            G = dx*u[-1,:]

            v[i,:] = linalg.solve(A2,G)

        else:
            ut = u[-i,:]

            v1 = sign(ut[:-1])
            v2 = sign(ut[1:])
            v3 = sign(ut[:-1] + ut[1:])

            data = zeros((2,N))
            data[0,:-1] = -ut[:-1]*0.5*(-v1+1)*0.5*(-v2+1) + \
                          -ut[:-1]*0.5*(-v1+1)*0.5*( v2+1)*0.5*(-v3+1)
            data[1,1: ] = -ut[1:] *0.5*( v1+1)*0.5*( v2+1) + \
                          -ut[1:] *0.5*(-v1+1)*0.5*( v2+1)*0.5*( v3+1)
            dfp_du = scipy.sparse.dia_matrix((data, [0, 1]), shape=(N, N))
            
            data = zeros((2,N))
            data[0,:-1] = -ut[:-1]*0.5*(-v1+1)*0.5*(-v2+1) + \
                          -ut[:-1]*0.5*(-v1+1)*0.5*( v2+1)*0.5*(-v3+1)
            data[1,1:]  = -ut[1:] *0.5*( v1+1)*0.5*( v2+1) + \
                          -ut[1:] *0.5*(-v1+1)*0.5*( v2+1)*0.5*( v3+1)
            data[1,0] = -ut[0]
            dfm_du = scipy.sparse.dia_matrix((data, [-1, 0]), shape=(N, N))

            dF_du = 1./dx*(dfm_du - dfp_du)
            A1 = -eye(N)/dt - dF_du

            G = zeros(N)

            if i == len(t)-1:
                v[i,:] = (G - dot(A1.T,v[i-1,:])).squeeze()
            else:
                v[i,:] = linalg.solve(A2,(G - dot(A1.T,v[i-1,:])).T).squeeze()

        i += 1

    return v

def J(dx, u):

    J = 0.5*dx*sum(u**2)

    return J

# problem parameters
N = 100
x = linspace(0.0,1.0,N+1)
dx = x[1] - x[0]
T = 5
t = linspace(0.0,T,500)

xcc = 0.5*(x[1:] + x[:-1])

# solve primal problem
u0 = 1.0 - xcc
u0[(xcc < 0.8)] = 0.0
u = burgers(u0, dx, t)

# solve the adjoint problem
v = adj_burgers(u, dx, t)

# compute the sensitivity to the initial condition with complex step
h = 1.0e-30
dJ_du0_fd = zeros(N)
for i in range(N):
    u0p = zeros(N,dtype=complex)
    u0p += u0
    u0p[i] += complex(0,h)
    up = burgers_complex(u0p, dx, t)
    Jp = J(dx, up[-1,:])
    dJ_du0_fd[i] = imag(Jp)/h

# plot the results    
i = 0
for ti in t:
    if mod(i,2) == 0:
        plt.ion()
        plt.show()
        plt.figure(1)
        ax1 = plt.gca()
        ax2 = ax1.twinx()
        ax1.plot(flipud(xcc),u[i,:],'b',linewidth=3)
        ax1.set_xlabel('x',fontsize=20)
        ax1.set_ylabel('u', color='b',fontsize=20)
        for t1 in ax1.get_yticklabels():
            t1.set_color('b')
        ax1.set_ylim([0,0.2])

        ax2.plot(flipud(xcc),v[i,:]/dx,'r',linewidth=3,label="Adjoint")
        for t2 in ax2.get_yticklabels():
            t2.set_color('r')
        ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        ax2.set_ylabel('v', color='r',fontsize=20)
        ax2.set_ylim([0,8.0e-4])

        if (ti == t[-1]) or (ti == t[-2]):
            ax2.plot(flipud(xcc),dJ_du0_fd,'k--',linewidth=3,label="Complex Step FD")
            plt.legend()

        plt.draw()
        plt.savefig('burgers_'+str(i/2)+'.pdf',format='pdf')
        plt.clf()
    
    i += 1
