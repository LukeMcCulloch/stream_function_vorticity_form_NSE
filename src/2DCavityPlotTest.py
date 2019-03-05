#
# Luke McCulloch
# Plot the 2D driven cavity flow
# NSE
# Stream Function, Vorticity Formulation
# Nov 2012

import numpy as np
from pylab import *
import matplotlib.pyplot as plt


BC=False

origin = 'lower'

stepsize=np.loadtxt('util.dat')

startx = 0.0
stopx  = 1.+stepsize[0]
stepx  = stepsize[0]

starty = 1.0
stopy  = 0.-stepsize[1]
stepy  = -stepsize[1]

X,Y = meshgrid( arange(startx,stopx,stepx),arange(starty,stopy,stepy) )


utlm   = np.loadtxt('u.dat')
vtlm   = np.loadtxt('v.dat')
uB   = np.loadtxt('uB.dat')
vB   = np.loadtxt('vB.dat')
vortlm= np.loadtxt('vorticity.dat')
streamtlm = np.loadtxt('strm-func.dat')

U = (utlm)
V = (vtlm)

# Bounding Box:
bx=[0.,1.,1.,0.,0.]
by=[0.,0.,1.,1.,0.]



M=sqrt(pow(U, 2) + pow(V, 2))
#1
# http://matplotlib.org/examples/pylab_examples/quiver_demo.html
figure()
if len(X)==11:
    Q = quiver( X,Y, U, V, M, units='x', pivot='tip',width=.005, scale=1./.15)
    if BC==True:
        qk = quiverkey(Q, 0.5, .94, 1, r'$u = 1 \frac{m}{s}$', labelpos='W',
               fontproperties={'weight': 'bold'})
        qk = quiverkey(Q, 0.75, .94, 0, r'$v=0\frac{m}{sec}$ ', labelpos='W',color='w',
               fontproperties={'weight': 'bold'})
        #qk = quiverkey(Q, .08, .5, 0, r'  $0 \frac{m}{s}$', labelpos='W',
        #       fontproperties={'weight': 'bold'})
        
        qk = quiverkey(Q, .15, .5, 0, r'  $u=v=0\frac{m}{sec}$', labelpos='W',color='w',
               fontproperties={'weight': 'bold'})
        qk = quiverkey(Q, 1.01, .5, 0, r'  $u=v=0\frac{m}{sec}$', labelpos='W',color='w',
               fontproperties={'weight': 'bold'})
        qk = quiverkey(Q, 0.6, .1, 0, r'  $u=v=0\frac{m}{sec}$', labelpos='W',color='w',
               fontproperties={'weight': 'bold'})
elif len(X)==51:
    Q = quiver( X,Y, U, V, M, units='x', pivot='tip',width=.005, scale=2./.15)
elif len(X)==101:
    Q = quiver( X,Y, U, V, M, units='x', pivot='tip',width=.005, scale=3./.15)
else:
    Q = quiver( X,Y, U, V, M, units='x', pivot='tip',width=.005, scale=3.3/.15)
      


l,r,b,t = axis()
dx, dy = r-l, t-b

#axis([-.2,1.2,-.2,1.2])
if BC==True:
    title('U & V, Vector Magnitude Initial Conditions')
    axis([l-0.2*dx, r+0.2*dx, b-0.2*dy, t+0.2*dy])
else:
    title('U & V, Vector Magnitudes')
    axis([l-0.05*dx, r+0.05*dx, b-0.05*dy, t+0.05*dy])
CB = plt.colorbar(Q, shrink=0.8, extend='both')
xlabel('X location')
ylabel('Y location')
plt.plot(bx,by)
plt.axis('equal')
show()



if BC==False:

    #levels=np.linspace(amin(streamtlm),amax(streamtlm),50,endpoint=True)
    levels=np.linspace(-.105,0.,50,endpoint=True)
    #c = plt.contour(X, Y,streamtlm,10)# good for Re 10
    c = plt.contour(X, Y,streamtlm,levels)
    plt.clabel(c, inline=1, fontsize=5)
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    mytitle = r"Stream Function"
    plt.title(mytitle)
    # we switch on a grid in the figure for orientation
    plt.grid()
    # colorbar
    CB = plt.colorbar(c, shrink=0.8, extend='both')
    plt.axis('equal')
    plt.plot(bx,by)
    plt.show()


    #levels=np.linspace(amin(vortlm),amax(vortlm),50,endpoint=True)
    levels=np.linspace(-144.,75.,50,endpoint=True)
    #levels=np.linspace(-150,75.
    #c = plt.contour(X, Y,vortlm,10)# good for Re 10
    c = plt.contour(X, Y,vortlm,30)
    plt.clabel(c, inline=1, fontsize=10)
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    mytitle = r"Vorticity"
    plt.title(mytitle)
    # we switch on a grid in the figure for orientation
    plt.grid()
    #colorbar
    CB = plt.colorbar(c, shrink=0.8, extend='both')
    plt.axis('equal')
    plt.plot(bx,by)
    plt.show()
