import math as m
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

# airfoil parameters
M = 0.02  # maximum camber [% chord]
P = 0.40  # maximum camber position [% chord]
t = 0.12  # thickness [% chord]
n = 50  # half number of panels
C = 4.0  # chord length [m]

# flow parameters
V = 1.0  # freestream velocity [m/s]
alpha = m.radians(5.0)  # angle of attack [rad]


def solve_airfoil(M, P, t, alpha):
    # create airfoil
    x = (1 - np.cos(np.linspace(0, m.pi, n)))/2
    y = np.zeros(n)

    for i in range(len(x)):
        if x[i] < P:
            y[i] = (M/P**2)*(2*P*x[i] - x[i]**2)
        else:
            y[i] = (M/(1 - P)**2)*(1 - 2*P + 2*P*x[i] - x[i]**2)

    a0 = 0.2969
    a1 = -0.1260
    a2 = -0.3516
    a3 = 0.2843
    a4 = -0.1036  # closed trailing edge
    yt = (t/0.2)*(a0*x**0.5 + a1*x + a2*x**2 + a3*x**3 + a4*x**4)
    yu = y + yt
    yl = y - yt

    X = C*np.concatenate((np.flip(x[1:]), x))
    Y = C*np.concatenate((np.flip(yl[1:]), yu))

    # panel geometry
    xc = 0.5*(X[1:] + X[:-1])
    yc = 0.5*(Y[1:] + Y[:-1])
    dx = X[1:] - X[:-1]
    dy = Y[1:] - Y[:-1]
    ds = np.sqrt(dx*dx + dy*dy)
    theta = np.arctan2(dy, dx)

    # influence coefficients
    N = 2*(n - 1)  # number of control points
    A = np.zeros((N + 1, N + 1))
    B = np.zeros(N + 1)
    C1 = np.zeros((N, N))
    C2 = np.zeros((N, N))
    D1 = np.zeros((N, N))
    D2 = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            a = -(xc[i] - X[j])*m.cos(theta[j]) - (yc[i] - Y[j])*m.sin(theta[j])
            b = (xc[i] - X[j])*(xc[i] - X[j]) + (yc[i] - Y[j])*(yc[i] - Y[j])
            c = m.sin(theta[i] - theta[j])
            d = m.cos(theta[i] - theta[j])
            e = (xc[i] - X[j])*m.sin(theta[j]) - (yc[i] - Y[j])*m.cos(theta[j])
            f = m.log(1 + (ds[j]*ds[j] + 2*a*ds[j])/b)
            g = m.atan2(e*ds[j], b + a*ds[j])
            p = (xc[i] - X[j])*m.sin(theta[i] - 2*theta[j]) + (yc[i] - Y[j])*m.cos(theta[i] - 2*theta[j])
            q = (xc[i] - X[j])*m.cos(theta[i] - 2*theta[j]) - (yc[i] - Y[j])*m.sin(theta[i] - 2*theta[j])

            C2[i, j] = d + 0.5*q*f/ds[j] - (a*c + d*e)*g/ds[j]
            C1[i, j] = 0.5*d*f + c*g - C2[i, j]
            D2[i, j] = c + 0.5*p*f/ds[j] + (a*d - c*e)*g/ds[j]
            D1[i, j] = 0.5*c*f - d*g - D2[i, j]

    np.fill_diagonal(C1, -1.0)
    np.fill_diagonal(C2, 1.0)
    A[:-1, 0] = C1[:, 0]
    A[:-1, 1:-1] = C1[:, 1:] + C2[:, :-1]
    A[:-1, -1] = C2[:, -1]
    B[:-1] = 2*m.pi*V*np.sin(theta - alpha)

    # Kutta condition
    A[-1, 0] = 1
    A[-1, -1] = 1

    # circulation density
    gamma = np.linalg.solve(A, B)
    cl = 2/(V*C)*sum(0.5*(gamma[1:] + gamma[:-1])*ds)

    # pressure distribution
    np.fill_diagonal(D1, m.pi/2)
    np.fill_diagonal(D2, m.pi/2)
    Vt = np.cos(theta - alpha) + 1/(2*m.pi*V)*(D1@gamma[:-1] + D2@gamma[1:])
    cp = 1 - Vt*Vt

    return X, Y, xc, yc, cl, cp


# show airfoil
X, Y, xc, yc, cl, cp = solve_airfoil(M, P, t, alpha)
fig, ax = plt.subplots(figsize=(7, 7))
plt.suptitle('Linear Vortex Panel Method', fontweight='bold')
ax.set_title('NACA %d%d%02d Airfoil\n'r'$c_l$ = %.2f, $\alpha$ = %.2f$\degree$' % (round(100*M), round(10*P), round(100*t), cl, m.degrees(alpha)))
line1, = plt.plot(X, Y, color='black', label='airfoil')
line2, = plt.plot(xc[n:], cp[n:], color='red', label='upper')
line3, = plt.plot(xc[:(n + 1)], cp[:(n + 1)], color='blue', label='lower')
plt.axis('equal')
plt.grid('on')
plt.legend(loc='lower right')

# create slider
plt.subplots_adjust(bottom=0.3)
axM = plt.axes([0.25, 0.2, 0.5, 0.03])
axP = plt.axes([0.25, 0.15, 0.5, 0.03])
axt = plt.axes([0.25, 0.1, 0.5, 0.03])
axalpha = plt.axes([0.25, 0.05, 0.5, 0.03])
sM = Slider(ax=axM, label='Camber', valmin=0, valmax=0.06, valinit=M)
sP = Slider(ax=axP, label='Position', valmin=0.3, valmax=0.7, valinit=P)
st = Slider(ax=axt, label='Thickness', valmin=0, valmax=0.15, valinit=t)
salpha = Slider(ax=axalpha, label='Alpha', valmin=0, valmax=0.2, valinit=alpha)


# update slider
def update(val):
    M = sM.val
    P = sP.val
    t = st.val
    alpha = salpha.val

    X, Y, xc, yc, cl, cp = solve_airfoil(M, P, t, alpha)
    plt.suptitle('Linear Vortex Panel Method', fontweight='bold')
    ax.set_title('NACA %d%d%02d Airfoil\n'r'$c_l$ = %.2f, $\alpha$ = %.2f$\degree$' % (round(100*M), round(10*P), round(100*t), cl, m.degrees(alpha)))
    line1.set_xdata(X)
    line1.set_ydata(Y)
    line2.set_xdata(xc[n:])
    line2.set_ydata(cp[n:])
    line3.set_xdata(xc[:(n + 1)])
    line3.set_ydata(cp[:(n + 1)])

    fig.canvas.draw_idle()


sM.on_changed(update)
sP.on_changed(update)
st.on_changed(update)
salpha.on_changed(update)

# reset slider
resetax = plt.axes([0.85, 0.05, 0.1, 0.05])
button = Button(resetax, 'Reset')


def reset(event):
    sM.reset()
    sP.reset()
    st.reset()
    salpha.reset()


button.on_clicked(reset)

plt.show()
