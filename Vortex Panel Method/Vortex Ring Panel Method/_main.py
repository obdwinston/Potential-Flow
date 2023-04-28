import math as m
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# airfoil parameters
Nx = 10  # number of panels in x-direction (even or odd number)
M = 0.06  # maximum camber [% chord]
P = 0.40  # maximum camber position [% chord]
c = 1.0  # chord length [m]

# wing parameters
Ny = 50  # number of panels in y-direction (even number)
a_le = m.radians(35)  # leading edge angle [rad]
a_te = m.radians(25)  # trailing edge angle [rad]
a_tt = m.radians(-5)  # angle of twist [rad]
b = 5.0  # wing span [m]
d = m.radians(10)  # dihedral angle [rad]

# flow parameters
V_inf = 10.0  # freestream velocity [m/s]
a_inf = m.radians(5)  # freestream angle of attack [rad]
b_inf = m.radians(0)  # freestream sideslip angle [rad]
rho_inf = 1.225  # freestream density [kg/m3]

# program options
show_wing = 1  # shows wing
show_control = 0  # shows control points on wing
show_ring = 0  # shows vortex rings on wing
show_normal = 0  # shows panel normals on wing

elevation = 40  # elevation of saved figures [deg]
azimuth = 200  # azimuth of saved figures [deg]

# program start

# freestream speeds
u_inf = V_inf*m.cos(a_inf)*m.cos(b_inf)  # freestream speed in x-direction
v_inf = V_inf*m.cos(a_inf)*m.sin(b_inf)  # freestream speed in y-direction
w_inf = V_inf*m.sin(a_inf)  # freestream speed in z-direction

# airfoil coordinates
x = (1 - np.cos(np.linspace(0, m.pi, Nx + 1)))/2  # x-coordinates of mean camber line
z = np.zeros(Nx + 1)  # z-coordinates of mean camber line
for i in range(Nx + 1):
    if x[i] < P:
        z[i] = (M/P**2)*(2*P*x[i] - x[i]**2)
    else:
        z[i] = (M/(1 - P)**2)*(1 - 2*P + 2*P*x[i] - x[i]**2)

x = c*x  # scale by chord length
z = c*z  # scale by chord length

# wing coordinates
y = np.linspace(-b/2, b/2, Ny + 1)  # y-sections

xx = np.zeros((Nx + 1, Ny + 1))  # x-coordinates of wing mean camber surface
yy = np.zeros((Nx + 1, Ny + 1))  # y-coordinates of wing mean camber surface
zz = np.zeros((Nx + 1, Ny + 1))  # z-coordinates of wing mean camber surface

for j in range(Ny + 1):
    # scale local
    xle = abs(y[j])*m.tan(a_le)  # local leading edge
    xte = abs(y[j])*m.tan(a_te) + c  # local trailing edge
    sc = (xte - xle)/c  # local scale
    xx[:, j] = sc*x + xle
    yy[:, j] = y[j]*np.ones(Nx + 1)
    zz[:, j] = sc*z + abs(y[j])*m.tan(d)
    # rotate local
    tt = abs(y[j])/(b/2)*a_tt  # local angle of twist
    dx = np.zeros(Nx)
    dz = np.zeros(Nx)
    xt = np.zeros(Nx)
    zt = np.zeros(Nx)
    dx[:] = xx[1:, j] - xx[0, j]
    dz[:] = zz[1:, j] - zz[0, j]
    xt[:] = xx[0, j] + dx*m.cos(tt) + dz*m.sin(tt)
    zt[:] = zz[0, j] - dx*m.sin(tt) + dz*m.cos(tt)
    xx[1:, j] = xt
    zz[1:, j] = zt

# control points
xc = np.zeros((Nx, Ny))
yc = np.zeros((Nx, Ny))
zc = np.zeros((Nx, Ny))
xcc = np.zeros((Nx + 1, Ny))
ycc = np.zeros((Nx + 1, Ny))
zcc = np.zeros((Nx + 1, Ny))

xcc[:, :] = 0.50*(xx[:, :-1] + xx[:, 1:])
ycc[:, :] = 0.50*(yy[:, :-1] + yy[:, 1:])
zcc[:, :] = 0.50*(zz[:, :-1] + zz[:, 1:])

xc[:, :] = xcc[:-1, :] + 0.75*(xcc[1:, :] - xcc[:-1, :])
yc[:, :] = ycc[:-1, :] + 0.50*(ycc[1:, :] - ycc[:-1, :])
zc[:, :] = zcc[:-1, :] + 0.75*(zcc[1:, :] - zcc[:-1, :])

# vortex rings
# vortex ring point 1
xv1 = np.zeros((Nx, Ny))
yv1 = np.zeros((Nx, Ny))
zv1 = np.zeros((Nx, Ny))

xv1[:-1, :] = xx[1:-1, :-1] + 0.25*(xx[2:, :-1] - xx[1:-1, :-1])
yv1[:-1, :] = yy[1:-1, :-1] + 0.25*(yy[2:, :-1] - yy[1:-1, :-1])
zv1[:-1, :] = zz[1:-1, :-1] + 0.25*(zz[2:, :-1] - zz[1:-1, :-1])

xv1[-1, :] = xx[-1, :-1] + 0.25*(xx[-1, :-1] - xx[-2, :-1])
yv1[-1, :] = yy[-1, :-1] + 0.25*(yy[-1, :-1] - yy[-2, :-1])
zv1[-1, :] = zz[-1, :-1] + 0.25*(zz[-1, :-1] - zz[-2, :-1])
# vortex ring point 2
xv2 = np.zeros((Nx, Ny))
yv2 = np.zeros((Nx, Ny))
zv2 = np.zeros((Nx, Ny))

xv2[:, :] = xx[:-1, :-1] + 0.25*(xx[1:, :-1] - xx[:-1, :-1])
yv2[:, :] = yy[:-1, :-1] + 0.25*(yy[1:, :-1] - yy[:-1, :-1])
zv2[:, :] = zz[:-1, :-1] + 0.25*(zz[1:, :-1] - zz[:-1, :-1])
# vortex ring point 3
xv3 = np.zeros((Nx, Ny))
yv3 = np.zeros((Nx, Ny))
zv3 = np.zeros((Nx, Ny))

xv3[:, :] = xx[:-1, 1:] + 0.25*(xx[1:, 1:] - xx[:-1, 1:])
yv3[:, :] = yy[:-1, 1:] + 0.25*(yy[1:, 1:] - yy[:-1, 1:])
zv3[:, :] = zz[:-1, 1:] + 0.25*(zz[1:, 1:] - zz[:-1, 1:])
# vortex ring point 4
xv4 = np.zeros((Nx, Ny))
yv4 = np.zeros((Nx, Ny))
zv4 = np.zeros((Nx, Ny))

xv4[:-1, :] = xx[1:-1, 1:] + 0.25*(xx[2:, 1:] - xx[1:-1, 1:])
yv4[:-1, :] = yy[1:-1, 1:] + 0.25*(yy[2:, 1:] - yy[1:-1, 1:])
zv4[:-1, :] = zz[1:-1, 1:] + 0.25*(zz[2:, 1:] - zz[1:-1, 1:])

xv4[-1, :] = xx[-1, 1:] + 0.25*(xx[-1, 1:] - xx[-2, 1:])
yv4[-1, :] = yy[-1, 1:] + 0.25*(yy[-1, 1:] - yy[-2, 1:])
zv4[-1, :] = zz[-1, 1:] + 0.25*(zz[-1, 1:] - zz[-2, 1:])

# panel normals
nx = np.zeros((Nx, Ny))
ny = np.zeros((Nx, Ny))
nz = np.zeros((Nx, Ny))
for j in range(Ny):
    for i in range(Nx):
        n1 = np.array([xx[i + 1, j] - xx[i, j + 1], yy[i + 1, j] - yy[i, j + 1], zz[i + 1, j] - zz[i, j + 1]])
        n2 = np.array([xx[i + 1, j + 1] - xx[i, j], yy[i + 1, j + 1] - yy[i, j], zz[i + 1, j + 1] - zz[i, j]])
        n3 = np.cross(n1, n2)
        nx[i, j] = n3[0]/np.linalg.norm(n3)
        ny[i, j] = n3[1]/np.linalg.norm(n3)
        nz[i, j] = n3[2]/np.linalg.norm(n3)

# wing area
dS = np.zeros((Nx, Ny))
for j in range(Ny):
    for i in range(Nx):
        v1 = np.array([xx[i + 1, j] - xx[i, j + 1], yy[i + 1, j] - yy[i, j + 1], zz[i + 1, j] - zz[i, j + 1]])
        v2 = np.array([xx[i + 1, j + 1] - xx[i, j], yy[i + 1, j + 1] - yy[i, j], zz[i + 1, j + 1] - zz[i, j]])
        v3 = np.cross(v1, v2)
        dS[i, j] = 0.50*np.linalg.norm(v3)
S = sum(sum(dS))

# show wing
if show_wing:
    # create plot
    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')
    ax.view_init(elev=elevation, azim=azimuth)
    # axis equal
    X, Y, Z = [-0.5, c + 0.5], [-b/2 - 0.5, b/2 + 0.5], [-0.5, 0.5]
    max_range = np.amax(np.array([np.amax(X) - np.amin(X), np.amax(Y) - np.amin(Y), np.amax(Z) - np.amin(Z)]))
    Xb = 0.5*max_range*np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + 0.5*(np.amax(X) + np.amin(X))
    Yb = 0.5*max_range*np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + 0.5*(np.amax(Y) + np.amin(Y))
    Zb = 0.5*max_range*np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + 0.5*(np.amax(Z) + np.amin(Z))
    for X_value, Y_value, Z_value in zip(Xb, Yb, Zb):
        ax.plot([X_value], [Y_value], [Z_value])
    # plot wing
    for i in range(Nx):
        for j in range(Ny):
            # mean camber surface
            x_fill = [xx[i, j], xx[i, j + 1], xx[i + 1, j + 1], xx[i + 1, j], xx[i, j]]
            y_fill = [yy[i, j], yy[i, j + 1], yy[i + 1, j + 1], yy[i + 1, j], yy[i, j]]
            z_fill = [zz[i, j], zz[i, j + 1], zz[i + 1, j + 1], zz[i + 1, j], zz[i, j]]
            vertices = [list(zip(x_fill, y_fill, z_fill))]
            panel = Poly3DCollection(vertices, linewidths=0.5, edgecolors='k', alpha=0.5)
            ax.add_collection3d(panel)
            # control points
            if show_control:
                ax.scatter3D(xc[i, j], yc[i, j], zc[i, j], c='red', marker='.')
            # vortex rings
            if show_ring:
                xv = [xv1[i, j], xv2[i, j], xv3[i, j], xv4[i, j], xv1[i, j]]
                yv = [yv1[i, j], yv2[i, j], yv3[i, j], yv4[i, j], yv1[i, j]]
                zv = [zv1[i, j], zv2[i, j], zv3[i, j], zv4[i, j], zv1[i, j]]
                ax.plot3D(xv, yv, zv, c='red', linewidth=0.5)
            # panel normals
            if show_normal:
                xn = [xc[i, j], xc[i, j] + 0.1*nx[i, j]]
                yn = [yc[i, j], yc[i, j] + 0.1*ny[i, j]]
                zn = [zc[i, j], zc[i, j] + 0.1*nz[i, j]]
                ax.plot3D(xn, yn, zn, c='green', linewidth=0.5)
    # show plot
    plt.grid('on')
    plt.show()
    plt.close()


# vortex functions
def vortex_line(p, pa, pb):
    xp, yp, zp = p[0], p[1], p[2]
    xa, ya, za = pa[0], pa[1], pa[2]
    xb, yb, zb = pb[0], pb[1], pb[2]

    r0 = np.array([xb - xa, yb - ya, zb - za])  # a pointing to b
    r1 = np.array([xp - xa, yp - ya, zp - za])  # a pointing to p
    r2 = np.array([xp - xb, yp - yb, zp - zb])  # b pointing to p
    r3 = np.cross(r1, r2)

    limit = 1e-10
    if np.linalg.norm(r1) < limit or np.linalg.norm(r2) < limit or np.linalg.norm(r3)**2 < limit:
        u = 0
        v = 0
        w = 0
    else:
        K1 = 1/(4*m.pi*np.linalg.norm(r3)**2)
        K2 = np.dot(r0, r1)/np.linalg.norm(r1)
        K3 = np.dot(r0, r2)/np.linalg.norm(r2)
        K = K1*(K2 - K3)

        u = K*r3[0]
        v = K*r3[1]
        w = K*r3[2]

    return np.array([u, v, w])


def vortex_ring(p, pq, pr, ps, pt):
    V1 = vortex_line(p, pq, pr)
    V2 = vortex_line(p, pr, ps)
    V3 = vortex_line(p, ps, pt)
    V4 = vortex_line(p, pt, pq)
    V = V1 + V2 + V3 + V4

    return V


# initialise variables
pc = np.zeros((Ny*Nx, 3))
nc = np.zeros((Ny*Nx, 3))
pv1 = np.zeros((Ny*Nx, 3))
pv2 = np.zeros((Ny*Nx, 3))
pv3 = np.zeros((Ny*Nx, 3))
pv4 = np.zeros((Ny*Nx, 3))

k = 0
for i in range(Nx):
    for j in range(Ny):
        pc[k] = np.array([xc[i, j], yc[i, j], zc[i, j]])
        nc[k] = np.array([nx[i, j], ny[i, j], nz[i, j]])
        pv1[k] = np.array([xv1[i, j], yv1[i, j], zv1[i, j]])
        pv2[k] = np.array([xv2[i, j], yv2[i, j], zv2[i, j]])
        pv3[k] = np.array([xv3[i, j], yv3[i, j], zv3[i, j]])
        pv4[k] = np.array([xv4[i, j], yv4[i, j], zv4[i, j]])
        k += 1

# steady wake
pwl = np.zeros((Ny + 1, 3))  # wake leading edge
pwl[:-1] = pv1[-Ny:]
pwl[-1] = pv4[-1]

pwt = np.zeros((Ny + 1, 3))  # wake trailing edge
pwt[:-1] = pv1[-Ny:]
pwt[-1] = pv4[-1]

dt = 20*b/u_inf
for i in range(len(pwt)):
    pwt[i] += dt*np.array([u_inf, v_inf, w_inf])

pwc = np.concatenate((pwl, pwt), axis=0)  # wake corner points

# wake vortex rings
Nwx = 1
Nwy = Ny
pw1 = np.zeros((Nwx*Nwy, 3))
pw2 = np.zeros((Nwx*Nwy, 3))
pw3 = np.zeros((Nwx*Nwy, 3))
pw4 = np.zeros((Nwx*Nwy, 3))
# wake vortex ring point 1
i = 0
for k in range((Nwx + 1)*(Nwy + 1)):
    if k > Nwy and (k + 1) % (Nwy + 1) != 0:
        pw1[i] = pwc[k]
        i += 1
# wake vortex ring point 2
i = 0
for k in range((Nwx + 1)*(Nwy + 1)):
    if k < (Nwx + 1)*(Nwy + 1) - (Nwy + 1) and (k + 1) % (Nwy + 1) != 0:
        pw2[i] = pwc[k]
        i += 1
# wake vortex ring point 3
i = 0
for k in range((Nwx + 1)*(Nwy + 1)):
    if k < (Nwx + 1)*(Nwy + 1) - (Nwy + 1) and k % (Nwy + 1) != 0:
        pw3[i] = pwc[k]
        i += 1
# wake vortex ring point 4
i = 0
for k in range((Nwx + 1)*(Nwy + 1)):
    if k > Nwy and k % (Nwy + 1) != 0:
        pw4[i] = pwc[k]
        i += 1

# influence coefficient matrix
A = np.zeros((Nx*Ny, Nx*Ny))
for i in range(Nx*Ny):
    k = 0
    for j in range(Nx*Ny):
        A[i, j] = np.dot(vortex_ring(pc[i], pv1[j], pv2[j], pv3[j], pv4[j]), nc[i])
        if j >= (Nx*Ny - Ny):
            A[i, j] += np.dot(vortex_ring(pc[i], pw1[k], pw2[k], pw3[k], pw4[k]), nc[i])
            k += 1

# constants vector
B = np.zeros(Nx*Ny)
for i in range(Nx*Ny):
    B[i] = -np.dot(np.array([u_inf, v_inf, w_inf]), nc[i])

# solution vector
gam = np.linalg.solve(A, B)
gam = gam.reshape(Nx, Ny)

# lift computation
dL = np.zeros((Nx, Ny))
dP = np.zeros((Nx, Ny))
for j in range(Ny):
    for i in range(Nx):
        if i == 0:
            dL[i, j] = rho_inf*V_inf*gam[i, j]*(yv3[i, j] - yv2[i, j])
            dP[i, j] = dL[i, j]/dS[i, j]
        else:
            dL[i, j] = rho_inf*V_inf*(gam[i, j] - gam[i - 1, j])*(yv3[i, j] - yv2[i, j])
            dP[i, j] = dL[i, j]/dS[i, j]

# plot pressure distribution
fig = plt.figure(figsize=(12, 10))
fig.suptitle('Pressure Distribution', fontweight='bold')
ax = plt.axes(projection='3d')
ax.view_init(elev=elevation, azim=azimuth)

X, Y, Z = [-0.5, c + 0.5], [-b/2 - 0.5, b/2 + 0.5], [-0.5, 0.5]
max_range = np.amax(np.array([np.amax(X) - np.amin(X), np.amax(Y) - np.amin(Y), np.amax(Z) - np.amin(Z)]))
Xb = 0.5*max_range*np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + 0.5*(np.amax(X) + np.amin(X))
Yb = 0.5*max_range*np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + 0.5*(np.amax(Y) + np.amin(Y))
Zb = 0.5*max_range*np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + 0.5*(np.amax(Z) + np.amin(Z))
for X_value, Y_value, Z_value in zip(Xb, Yb, Zb):
    ax.plot([X_value], [Y_value], [Z_value])

cmap = plt.get_cmap('jet')
norm = clr.Normalize(vmin=np.max(np.min(dP)), vmax=np.max(np.max(dP)))
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
plt.colorbar(sm, orientation='vertical', label='Pressure [Pa]', shrink=0.5)

for i in range(Nx):
    for j in range(Ny):
        x_fill = [xx[i, j], xx[i, j + 1], xx[i + 1, j + 1], xx[i + 1, j], xx[i, j]]
        y_fill = [yy[i, j], yy[i, j + 1], yy[i + 1, j + 1], yy[i + 1, j], yy[i, j]]
        z_fill = [zz[i, j], zz[i, j + 1], zz[i + 1, j + 1], zz[i + 1, j], zz[i, j]]
        vertices = [list(zip(x_fill, y_fill, z_fill))]
        panel = Poly3DCollection(vertices, edgecolors='None')
        panel.set_facecolor(cmap(norm(dP[i, j])))
        ax.add_collection3d(panel)

plt.grid('on')
plt.show()

# plot lift distribution
fig = plt.figure(figsize=(12, 6))
fig.suptitle('Lift Distribution', fontweight='bold')
plt.plot(yc[0, :]/b, sum(dL))
plt.title(r'AR = %.2f, L = %.2f N, $c_L$ = %.2f, $\alpha$ = %.2f$\degree$' % (b*b/S, sum(sum(dL)), (2*sum(sum(dL)))/(rho_inf*V_inf*V_inf*S), m.degrees(a_inf)))
plt.xlabel('y/b')
plt.ylabel('L/b')

plt.grid('on')
plt.show()
