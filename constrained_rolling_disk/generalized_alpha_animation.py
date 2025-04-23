## Importing packages
import numpy as np
import os
import argparse

# # optaining certain parametric arguments from an input file
# def parse_args():
#     parser = argparse.ArgumentParser(description='Process some integers and floats.')
#     parser.add_argument('weight', type=float, help='disk weight in lb (float)')
#     parser.add_argument('radius', type=float, help='disk radius in ft (float)')
#     parser.add_argument('theta', type=float, help='angle that disk make with vertical in rad (float)')
#     parser.add_argument('xbar0_1', type=float, help='x-coordinate of initial location of CG of disk (float)')
#     parser.add_argument('xbar0_2', type=float, help='y-coordinate of initial location of CG of disk (float)')
#     parser.add_argument('xbar0_3', type=float, help='z-coordinate of initial location of CG of disk (float)')
#     parser.add_argument('--output_path', '-o', type=str, help='Output path (optional)')
    
#     args = parser.parse_args()
    
#     return args.weight, args.radius, args.height, args.theta, args.I_xbar0_1, args.I_xbar0_2, args.I_xbar0_3, args.II_xbar0_1, args.II_xbar0_2, args.II_xbar0_3, args.output_path

# weight, radius, height, theta, I_xbar0_1, I_xbar0_2, I_xbar0_3, II_xbar0_1, II_xbar0_2, II_xbar0_3, output_path = parse_args()

weight = 1      # disk weight in lb (float)
radius = 1      # disk radius in ft (float)
theta = np.pi/3 # angle that disk make with vertical in rad (float)
xbar0_1 = 0     # x-coordinate of initial location of CG of disk (float)
xbar0_2 = 0     # x-coordinate of initial location of CG of disk (float)
xbar0_3 = radius*np.sin(np.pi/3)    # x-coordinate of initial location of CG of disk (float)
output_path = "/Users/theresahonein/Desktop/three_wheel_vehicle/constrained_rolling_disk/outputs"

directory_containing_output = os.path.dirname(output_path)

directory_containing_output = os.path.dirname(output_path)
f = open(f"{directory_containing_output}/saved_values_all.txt",'a')
g = open(f"{directory_containing_output}/saved_values_plots.txt",'a')

## Problem Constants
# number of degrees of freedom
ndof = 6

# cylinder dimensions (considering two identical cylinders)
gravity = 32.2          # m/sec^2, gravitational acceleration
mass = weight/gravity   # slugs, cylinder mass

# restitution coefficients
e_N = 0                 # dimensionless, normal impact restitution coefficient
e_F = 0                 # dimensionaless, tangential impact restitution coefficient

# friction coefficients
# source: https://www.engineeringtoolbox.com/friction-coefficients-d_778.html
mu_s = np.inf           # static friction coefficient
mu_k = 0.4              # kinetic friction coefficient

## Nondimensionalization
# Nondimensionalization Parameters
mass_nd_param = mass
length_nd_param = 1
time_nd_param = np.sqrt(length_nd_param/gravity)
velocity_nd_param = np.sqrt(gravity*length_nd_param)
acceleration_nd_param = gravity
force_nd_param = mass_nd_param*acceleration_nd_param

# Nondimensionalized Values
m = mass/mass_nd_param
R = radius/length_nd_param
h = 0
gr = 1                  # dimensionless, gravity

# Rotational Inertia Terms
lambdat = 1/4*m*R**2    # [mass].[length]**2, tranverse moment of inertia
lambdaa = 1/2*m*R**2    # [mass].[length]**2, axial moment of inertia
I = np.array([[lambdat,0,0],[0,lambdat,0],[0,0,lambdaa]])   # moment of intertia matrix

## Simulation Parameters
ti = 0                      # [time], initial time
ntime = 2000                # number of iterations
dtime = 2e-3/time_nd_param  # [time], time step duration

## Parameters of generalized alpha scheme
# Differentiation parameters
rho_inf = 0.5
alpha_m = (2*rho_inf-1)/(rho_inf+1)
alpha_f = rho_inf/(rho_inf+1)
gama = 0.5+alpha_f-alpha_m
beta = 0.25*(0.5+gama)**2

# Set approximation parameter
r = 0.3

# Loop parameters
maxiter_n = 100
tol_n = 1.0e-8

## Kinetic Quantities
ng = 1         # num position level constraints, constant
nN = 1        # num normal no penetration constraints, constant
ngamma = 1     # num velocity level constraints, contant
nF = 2    # no slip constraints/friction force, constant

## Initialize arrays to save results
a = np.zeros((ntime,ndof))
U = np.zeros((ntime,ndof))
Q = np.zeros((ntime,ndof))
Kappa_g = np.zeros((ntime,ng))
Lambda_g = np.zeros((ntime,ng))
lambda_g = np.zeros((ntime,ng))
Lambda_gamma = np.zeros((ntime,ngamma))
lambda_gamma = np.zeros((ntime,ngamma))
Kappa_N = np.zeros((ntime,nN))
Lambda_N = np.zeros((ntime,nN))
lambda_N = np.zeros((ntime,nN))
Lambda_F = np.zeros((ntime,nF))
lambda_F = np.zeros((ntime,nF))

R_array = np.zeros(ntime)

gamma_F = np.zeros((ntime,nF))
gdot_N = np.zeros((ntime,nN))
q = np.zeros((ntime,ndof))
u = np.zeros((ntime,ndof))

# obtained from > valid_touching_config
psi0 = 0
theta0 = theta
phi0 = 0

psidot0 = 1
thetadot0 = 0 

phidot = -6*np.pi*time_nd_param*np.ones(ntime)
phiddot = np.zeros(ntime)

xbar0 = np.array([xbar0_1, xbar0_2, xbar0_3])/length_nd_param

vbar0 = np.array([0,0,0])

q0 = np.concatenate((xbar0,psi0,theta0,phi0),axis=None)
u0 = np.concatenate((vbar0,0,0,0),axis=None)

# better obtain these values from static case
lambda_N[0,:] = m*gr 
lambda_g[0,:] = 0.2*m*gr

x0 = np.concatenate((a[0,:],U[0,:],Q[0,:],
                     Kappa_g[0,:],Lambda_g[0,:],lambda_g[0,:],
                     Lambda_gamma[0,:],lambda_gamma[0,:],
                     Kappa_N[0,:],Lambda_N[0,:],lambda_N[0,:],
                     Lambda_F[0,:],lambda_F[0,:]),axis=None)

# initial auxiliary variables
a_bar0 = np.zeros(ndof)
lambda_N_bar0 = lambda_N[0,:]
lambda_F_bar0 = np.zeros(nF)

prev_AV = np.concatenate((a_bar0,lambda_N_bar0,lambda_F_bar0),axis=None)

gamma_F0 = np.zeros((nF))
gdot_N0 = np.zeros((nN))

def get_x_components(x):
    a = x[0:ndof]
    U = x[ndof:2*ndof]
    Q = x[2*ndof:3*ndof]
    Kappa_g = x[3*ndof:3*ndof+ng]
    Lambda_g = x[3*ndof+ng:3*ndof+2*ng]
    lambda_g = x[3*ndof+2*ng:3*ndof+3*ng]
    Lambda_gamma = x[3*ndof+3*ng:3*ndof+3*ng+ngamma]
    lambda_gamma = x[3*ndof+3*ng+ngamma:3*ndof+3*ng+2*ngamma]
    Kappa_N = x[3*ndof+3*ng+2*ngamma:3*ndof+3*ng+2*ngamma+nN]
    Lambda_N = x[3*ndof+3*ng+2*ngamma+nN:3*ndof+3*ng+2*ngamma+2*nN]
    lambda_N = x[3*ndof+3*ng+2*ngamma+2*nN:3*ndof+3*ng+2*ngamma+3*nN]
    Lambda_F = x[3*ndof+3*ng+2*ngamma+3*nN:3*ndof+3*ng+2*ngamma+3*nN+nF]
    lambda_F = x[3*ndof+3*ng+2*ngamma+3*nN+nF:3*ndof+3*ng+2*ngamma+3*nN+2*nF]
    return a, U, Q, Kappa_g, Lambda_g, lambda_g, Lambda_gamma, lambda_gamma, Kappa_N, Lambda_N, lambda_N, Lambda_F, lambda_F


def get_disk_ground_contact(q,u,a):

    # Rotation Matricies
    R1 = np.array([[np.cos(q[3]), np.sin(q[3]), 0],
                   [-np.sin(q[3]), np.cos(q[3]), 0],
                   [0, 0, 1]])
    
    R2 = np.array([[1, 0, 0],
                   [0, np.cos(q[4]), np.sin(q[4])],
                   [0, -np.sin(q[4]), np.cos(q[4])]])
    
    R3 = np.array([[np.cos(q[5]), np.sin(q[5]), 0],
                   [-np.sin(q[5]), np.cos(q[5]), 0],
                   [0, 0, 1]])

    # Fixed basis vectors
    E1 = np.array([1,0,0])
    E2 = np.array([0,1,0])
    E3 = np.array([0,0,1])

    # Needed moving basis vectors
    e1pp = np.transpose(R2@R1)@E1
    e2pp = np.transpose(R2@R1)@E2
    e3 = np.transpose(R3@R2@R1)@E3

    # Relevant angular velocities
    omegapp = u[3]*E3+u[4]*e1pp
    omega = u[3]*E3+u[4]*e1pp+u[5]*e3

    # Relevant basis vectors first time derivatives
    e1ppdot = np.cross(omegapp,e1pp)
    e2ppdot = np.cross(omegapp,e2pp)
    e3dot = np.cross(omega,e3)
    
    # Relevant angular accelerations
    alphapp = a[3]*E3+a[4]*e1pp+u[4]*e1ppdot
    alpha = a[3]*E3+a[4]*e1pp+u[4]*e1ppdot+a[5]*e3+u[5]*e3dot

    # Relevant basis vectors second time derivatives 
    e2ppddot = np.cross(alphapp,e2pp)+np.cross(omegapp,e2ppdot)
    e3ddot = np.cross(alpha,e3)+np.cross(omega,e3dot)

    pi = -R*e2pp
    pidot = np.cross(omega,pi) #-R*e2ppdot
    piddot = np.cross(alpha,pi)+np.cross(omega,pidot) #+R*e2ppddot

    

    vC = np.array([-(u[3] + u[5]*np.cos(q[4]))*(2*R*np.cos(q[4]) - h*np.sin(q[4]))*np.cos(q[3])/2 + (2*R*np.sin(q[4]) + h*np.cos(q[4]))*(u[4]*np.sin(q[3]) - u[5]*np.sin(q[4])*np.cos(q[3]))/2,
                   -(u[3] + u[5]*np.cos(q[4]))*(2*R*np.cos(q[4]) - h*np.sin(q[4]))*np.sin(q[3])/2 - (2*R*np.sin(q[4]) + h*np.cos(q[4]))*(u[4]*np.cos(q[3]) + u[5]*np.sin(q[3])*np.sin(q[4]))/2,
                   u[4]*(R*np.cos(q[4]) - h*np.sin(q[4])/2)])

    aC = np.array([-R*a[3]*np.cos(q[3])*np.cos(q[4]) + R*a[4]*np.sin(q[3])*np.sin(q[4]) - R*a[5]*np.cos(q[3]) + R*u[3]**2*np.sin(q[3])*np.cos(q[4]) + 2*R*u[3]*u[4]*np.sin(q[4])*np.cos(q[3]) + R*u[3]*u[5]*np.sin(q[3]) + R*u[4]**2*np.sin(q[3])*np.cos(q[4]) + a[3]*h*np.sin(q[4])*np.cos(q[3])/2 + a[4]*h*np.sin(q[3])*np.cos(q[4])/2 - h*u[3]**2*np.sin(q[3])*np.sin(q[4])/2 + h*u[3]*u[4]*np.cos(q[3])*np.cos(q[4]) - h*u[4]**2*np.sin(q[3])*np.sin(q[4])/2,
                   -R*a[3]*np.sin(q[3])*np.cos(q[4]) - R*a[4]*np.sin(q[4])*np.cos(q[3]) - R*a[5]*np.sin(q[3]) - R*u[3]**2*np.cos(q[3])*np.cos(q[4]) + 2*R*u[3]*u[4]*np.sin(q[3])*np.sin(q[4]) - R*u[3]*u[5]*np.cos(q[3]) - R*u[4]**2*np.cos(q[3])*np.cos(q[4]) + a[3]*h*np.sin(q[3])*np.sin(q[4])/2 - a[4]*h*np.cos(q[3])*np.cos(q[4])/2 + h*u[3]**2*np.sin(q[4])*np.cos(q[3])/2 + h*u[3]*u[4]*np.sin(q[3])*np.cos(q[4]) + h*u[4]**2*np.sin(q[4])*np.cos(q[3])/2,
                   a[4]*(2*R*np.cos(q[4]) - h*np.sin(q[4]))/2 - u[4]**2*(2*R*np.sin(q[4]) + h*np.cos(q[4]))/2])


    g_N = q[2]+pi[2]
    gdot_N = u[2]+pidot[2]
    gddot_N = a[2]-aC[2]
    W_N = np.zeros(6)
    W_N[2] = 1
    W_N[4] = h/2*(np.sin(q[4]))-R*(np.cos(q[4]))

    gamma_F1 = u[0]-vC[0]
    gamma_F2 = u[1]-vC[1]

    gammadot_F1 = a[0]-aC[0]
    gammadot_F2 = a[1]-aC[1]

    W_F1 = np.zeros(6)
    W_F1[0] = 1
    W_F1[3] = (R*np.cos(q[4]) - h*np.sin(q[4])/2)*np.cos(q[3])
    W_F1[4] = (-R*np.sin(q[4]) - h*np.cos(q[4])/2)*np.sin(q[3])
    W_F1[5] = -(-R*np.sin(q[4]) - h*np.cos(q[4])/2)*np.sin(q[4])*np.cos(q[3]) + (2*R*np.cos(q[4]) - h*np.sin(q[4]))*np.cos(q[3])*np.cos(q[4])/2

    W_F2 = np.zeros(6)
    W_F2[1] = 1
    W_F2[3] = (R*np.cos(q[4]) - h*np.sin(q[4])/2)*np.sin(q[3])
    W_F2[4] = (R*np.sin(q[4]) + h*np.cos(q[4])/2)*np.cos(q[3])
    W_F2[5] = (R*np.sin(q[4]) + h*np.cos(q[4])/2)*np.sin(q[3])*np.sin(q[4]) + (2*R*np.cos(q[4]) - h*np.sin(q[4]))*np.sin(q[3])*np.cos(q[4])/2

    return g_N,gdot_N,gddot_N,W_N,gamma_F1,gammadot_F1,W_F1,\
        gamma_F2,gammadot_F2,W_F2


def get_R(x, prev_x, prev_AV, prev_gamma_F, prev_gdot_N, prev_q, prev_u,
           J_calc_bool):
    global iter
    global A, B, C, D_st, E_st
    global phiIdot

    # Data extraction
    prev_a, _, _, _, _, _, _, _, _, _, prev_lambda_N, _, prev_lambda_F = \
        get_x_components(prev_x)
    a, U, Q, Kappa_g, Lambda_g, lambda_g, Lambda_gamma, lambda_gamma,\
        Kappa_N, Lambda_N, lambda_N, Lambda_F, lambda_F = get_x_components(x)

    # Getting previous auxiliary variables
    prev_a_bar = prev_AV[0 : ndof]
    prev_lambda_N_bar = prev_AV[ndof : ndof+nN]
    prev_lambda_F_bar = prev_AV[ndof+nN : ndof+nN+nF]

    # Calculating new auxiliary variables
    a_bar = (alpha_f*prev_a+(1-alpha_f)*a -
             alpha_m*prev_a_bar)/(1-alpha_m)  # (71)
    lambda_N_bar = (alpha_f*prev_lambda_N+(1-alpha_f)*lambda_N
                    - alpha_m*prev_lambda_N_bar)/(1-alpha_m)  # (96)
    lambda_F_bar = (alpha_f*prev_lambda_F+(1-alpha_f)*lambda_F
                    - alpha_m*prev_lambda_F_bar)/(1-alpha_m)  # (114)
    AV = np.concatenate((a_bar, lambda_N_bar, lambda_F_bar), axis=None)

    # Calculate q and u (73)
    u = prev_u+dtime*((1-gama)*prev_a_bar+gama*a_bar)+U
    q = prev_q+dtime*prev_u+dtime**2/2*((1-2*beta)*prev_a_bar+2*beta*a_bar)+Q

    # Mass matrix
    Mdiag = np.array([m, m, m, lambdat, lambdat, lambdaa])
    M = np.diag(Mdiag)

    # Vector of applied forces and moments
    # applied vertical weight force
    f = np.array([0, 0, -m*gr, 0, 0, 0])

    # Constraint and constraint gradients
    # initializing constraints
    g = np.zeros(ng)
    gamma = np.zeros(ngamma)
    g_N = np.zeros(nN)
    gamma_F = np.zeros(nF)
    # initialize constraint derivatives
    gdot = np.zeros(ng)
    gddot = np.zeros(ng)
    gammadot = np.zeros(ngamma)
    gdot_N = np.zeros(nN)
    gddot_N = np.zeros(nN)
    gammadot_F = np.zeros(nF)
    # initializing constraint gradients
    W_g = np.zeros((ng, ndof))
    W_gamma = np.zeros((ngamma, ndof))
    W_N = np.zeros((nN, ndof))
    W_F = np.zeros((nF, ndof))

    # Constraints on the disk
    g_N[0], gdot_N[0], gddot_N[0], W_N[0, :],\
        gamma_F[0], gammadot_F[0], W_F[0, :],\
        gamma_F[1], gammadot_F[1], W_F[1, :]\
        = get_disk_ground_contact(q, u, a)

    # Constraints prescribing the motion of the disk

    # # angle I_psi
    # g[0] = q[3]-psi0
    # gdot[0] = u[3]
    # gddot[0] = a[3]
    # W_g[0,3] = 1

    # coordinate x3
    c1 = xbar0[2]/length_nd_param+h/2*np.cos(theta)
    g[0] = q[2]+h/2*np.cos(q[4])-c1
    gdot[0] = u[2]-h/2*u[4]*np.sin(q[4])
    gddot[0] = a[2]-h/2*(a[4]*np.sin(q[4])+u[4]**2*np.cos(q[4]))
    W_g[0,2] = 1
    W_g[0,4] = -h/2*np.sin(q[4])

    # The driving motion
    # angle I_phi
    gamma[0] = u[5]-phidot[iter]
    gammadot[0] = a[5]-phiddot[iter]
    W_gamma[0,5] = 1

    # Kinetic quantities
    # normal
    ksi_N = gdot_N+e_N*prev_gdot_N # (86)
    P_N = Lambda_N+dtime*((1-gama)*prev_lambda_N_bar+gama*lambda_N_bar) # (95)
    Kappa_hat_N = Kappa_N+dtime**2/2*((1-2*beta)*prev_lambda_N_bar+2*beta*lambda_N_bar) # (102)
    # frictional
    ksi_F = gamma_F+e_F*prev_gamma_F
    P_F = Lambda_F+dtime*((1-gama)*prev_lambda_F_bar+gama*lambda_F_bar) # (113)

    # Smooth residual Rs
    temp1 = M@a-f-np.transpose(W_g)@lambda_g-np.transpose(W_gamma)@lambda_gamma-np.transpose(W_N)@lambda_N-np.transpose(W_F)@lambda_F
    temp2 = M@U-np.transpose(W_g)@Lambda_g-np.transpose(W_gamma)@Lambda_gamma-np.transpose(W_N)@Lambda_N-np.transpose(W_F)@Lambda_F
    temp3 = M@Q-np.transpose(W_N)@Kappa_N-np.transpose(W_g)@Kappa_g-dtime/2*np.transpose(W_F)@Lambda_F-dtime/2*np.transpose(W_gamma)@Lambda_gamma
    Rs = np.concatenate((temp1,temp2,temp3,g,gdot,gddot,gamma,gammadot),axis=None)

    # Contact residual Rc
    R_Kappa_N = np.zeros(nN) # (129)
    R_Lambda_N = np.zeros(nN)
    R_lambda_N = np.zeros(nN)
    R_Lambda_F = np.zeros(nF) # (138)
    R_lambda_F = np.zeros(nF) # (142)

    gammaF_lim = np.array([[0,1],[2,3],[4,5]])

    for k in range(nN):
        if A[k]:
            R_Kappa_N[k] = g_N[k]
            if D_st[k]:
                R_Lambda_F[gammaF_lim[k,:]] = ksi_F[gammaF_lim[k,:]]
                if E_st[k]:
                    R_lambda_F[gammaF_lim[k,:]] = gammadot_F[gammaF_lim[k,:]]
                else:
                    R_lambda_F[gammaF_lim[k,:]] = lambda_F[gammaF_lim[k,:]]+mu_k*lambda_N[k]*np.sign(gammadot_F[gammaF_lim[k,:]])                    
            else:
                R_Lambda_F[gammaF_lim[k,:]] = P_F[gammaF_lim[k,:]]+mu_k*P_N[k]*np.sign(ksi_F[gammaF_lim[k,:]])
                R_lambda_F[gammaF_lim[k,:]] = lambda_F[gammaF_lim[k,:]]+mu_k*lambda_N[k]*np.sign(gamma_F[gammaF_lim[k,:]])
        else:
            R_Kappa_N[k] = Kappa_hat_N[k]
            R_Lambda_F[gammaF_lim[k,:]] = P_F[gammaF_lim[k,:]]
            R_lambda_F[gammaF_lim[k,:]] = lambda_F[gammaF_lim[k,:]]
        # (132)
        if B[k]:
            R_Lambda_N[k] = ksi_N[k]
        else:
            R_Lambda_N[k] = P_N[k]
        # (135)
        if C[k]:
            R_lambda_N[k] = gddot_N[k]
        else:
            R_lambda_N[k] = lambda_N[k]

    Rc = np.concatenate((R_Kappa_N,R_Lambda_N,R_lambda_N,R_Lambda_F,R_lambda_F),axis=None)

    # Assembling residual array
    Res = np.concatenate((Rs,Rc),axis=None)

    if J_calc_bool == False:
        norm_R = np.linalg.norm(Res,np.inf)
        print(f'norm_R = {norm_R}')

    return Res, AV, gdot_N, gamma_F, q, u

def get_R_J(x,prev_x,prev_AV,prev_gamma_F,prev_gdot_N,prev_q,prev_u):
    
    epsilon = 1e-6
    R_x, AV, gdot_N, gamma_F, q, u = get_R(x,prev_x,prev_AV,prev_gamma_F,
                                           prev_gdot_N,prev_q,prev_u,False)
    n = np.size(R_x) # Jacobian dimension
    # Initializing the Jacobian
    J = np.zeros((n,n))
    I = np.identity(n)
    # Constructing the Jacobian column by column
    for i in range(n):
        # print(i)
        R_x_plus_epsilon,_,_,_,_,_ = get_R(x+epsilon*I[:,i],prev_x,prev_AV,
                                           prev_gamma_F,prev_gdot_N,prev_q,prev_u,True)
        J[:,i] = (R_x_plus_epsilon-R_x)/epsilon

    return R_x,AV,gdot_N,gamma_F,q,u,J


## Solution
x_temp = x0
prev_x = x0

q[0,:] = q0
u[0,:] = u0
gamma_F[0,:] = gamma_F0
gdot_N[0,:] = gdot_N0

n_Rs = 3*ndof+3*ng+2*ngamma
n = np.size(x_temp)

y_indices = range(n_Rs)
z_indices = range(n_Rs,n)

A=[1]
B=[1]
C=[1]
D_st = [1]
E_st = [1]

for iter in range(1,ntime):
    print(f'i={iter}')

    # First semismooth Newton calculation
    t = dtime*iter
    nu = 0

    Res,AV_temp,gdot_N_temp,gamma_F_temp,q_temp,u_temp,J\
         = get_R_J(x_temp,prev_x,prev_AV,gamma_F[iter-1,:],gdot_N[iter-1,:],
                   q[iter-1,:],u[iter-1,:])
    
    norm_R = np.linalg.norm(Res,np.inf)
    print(f'lambda_N = {x_temp[3*ndof+3*ng+2*ngamma+2*nN:3*ndof+3*ng+2*ngamma+3*nN]}')
    print(f'lambda_g = {x_temp[3*ndof+2*ng:3*ndof+3*ng]}')

    # Semismooth Newton iterations
    while norm_R>tol_n and nu<maxiter_n:
        # Newton update
        x_temp = x_temp-np.linalg.solve(J,Res)
        # Calculate new EOM and residual
        nu = nu+1
        print(f'nu = {nu}')
        Res,AV_temp,gdot_N_temp,gamma_F_temp,q_temp,u_temp,J = \
            get_R_J(x_temp,prev_x,prev_AV,gamma_F[iter-1,:],gdot_N[iter-1,:],\
                    q[iter-1,:],u[iter-1,:])
        norm_R = np.linalg.norm(Res,np.inf)
        # print(f'norm_R = {norm_R}')
        print(f'lambda_N = {x_temp[3*ndof+3*ng+2*ngamma+2*nN:3*ndof+3*ng+2*ngamma+3*nN]}')
        print(f'lambda_g = {x_temp[3*ndof+2*ng:3*ndof+3*ng]}')
    
    R_array[iter] = norm_R

    # Saving output results
    a[iter,:],U[iter,:],Q[iter,:],Kappa_g[iter,:],Lambda_g[iter,:],lambda_g[iter,:],\
        Lambda_gamma[iter,:],lambda_gamma[iter,:],Kappa_N[iter,:],\
            Lambda_N[iter,:],lambda_N[iter,:],Lambda_F[iter,:],lambda_F[iter,:]\
                  = get_x_components(x_temp)
        
    # Updating reusable results
    gamma_F[iter,:] = gamma_F_temp
    gdot_N[iter,:] = gdot_N_temp
    q[iter,:] = q_temp
    u[iter,:] = u_temp
    prev_x = x_temp
    prev_AV = AV_temp 


f.write(f'{theta}\n {lambda_N}\n {lambda_g}\n')
g.write(f'{theta} {lambda_N[iter,0]}\n')

import scipy.io
file_name_q = str(f'{output_path}/q.mat')
file_name_u = str(f'{output_path}/u.mat')
file_name_a = str(f'{output_path}/a.mat')
file_name_g = str(f'{output_path}/g.mat')
file_name_gamma = str(f'{output_path}/gamma.mat')
file_name_F = str(f'{output_path}/F.mat')
file_name_N = str(f'{output_path}/N.mat')
scipy.io.savemat(file_name_q,dict(q=q))
scipy.io.savemat(file_name_u,dict(u=u))
scipy.io.savemat(file_name_a,dict(a=a))
scipy.io.savemat(file_name_g,dict(lambda_g=lambda_g))
scipy.io.savemat(file_name_gamma,dict(lambda_gamma=lambda_gamma))
scipy.io.savemat(file_name_N,dict(lambda_N=lambda_N))
scipy.io.savemat(file_name_F,dict(lambda_F=lambda_F))


print('done')