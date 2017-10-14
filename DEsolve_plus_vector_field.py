import numpy as np
import matplotlib.pyplot as plt

''' 
2D AB4 solver. Solve coupled Differential equations of the form
    x' = f(x,y)
    y' = f(x,y)

AB4 because it took me 2 seconds to set up and it seems to work well
'''

##################################################################################
'''
    QUESTION 2 SOLUTIONS: explicit
'''

def q2x(t,x0):
    c = x0*x0/(x0*x0 + 1.)
    a = c*np.exp(2.*t)
    b = 1. + a
    return np.sqrt(a/b)

def q2y(t,y0):
    return y0*np.exp(-t)



##################################################################################
# ENTER THE SYSTEM OF DEs HERE:
##################################################################################
def xprime(x,y):
    return  x*x + 3*y*y - 1

def yprime(x,y):
    return -2*x*y
        
def ab4_2d(x0,y0,f1,f2,timeseries):

    h = np.abs(timeseries[1] - timeseries[0])
    
    tn = timeseries[0]
    yn = y0
    xn = x0
    
    tn1 = tn + h
    yn1 = yn + h*f1(xn,yn)
    xn1 = xn + h*f2(xn,yn)    

    tn2 = tn1 + h
    yn2 = yn1 + h*f1(xn1,yn1)    
    xn2 = xn1 + h*f2(xn1,yn1)
    
    tn3 = tn2 + h
    yn3 = yn2 + h*f1(xn2,yn2)
    xn3 = xn2 + h*f2(xn2,yn2)
    
    # t is the time series, 'x' is the first de soln, 'y' is second de soln
    xvals,yvals,tvals = [],[],[]    
    for i in range(0,len(timeseries)):
        try:

            yn3 += (h/24.0)*( 55.0*(f2(xn3,yn3)) - 59.0*(f2(xn2,yn2)) + 37.0*(f2(xn1,yn1)) - 9.0*(f2(xn,yn)))
            xn3 += (h/24.0)*( 55.0*(f1(xn3,yn3)) - 59.0*(f1(xn2,yn2)) + 37.0*(f1(xn1,yn1)) - 9.0*(f1(xn,yn)))            
            
            yn = yn1
            tn = tn1
            xn = xn1
            
            yn1 = yn2
            tn1 = tn2
            xn1 = xn2

            yn2 = yn3
            tn2 = tn3
            xn2 = xn3
            
            tn3 += h  

            xvals.append(xn3)
            tvals.append(tn3)
            yvals.append(yn3)
            
            if np.abs(yn3) > 1e2:                
                # Dont want our solution to continue for ridiculous y-values and
                # blow up...
                #print 'large nums'
                break
            if np.abs(xn3) > 1e2:                
                            # Dont want our solution to continue for ridiculous y-values and
                            # blow up...
                            #print 'large nums'
                            break
        except OverflowError:
            print 'values blew up!!'
            break

    return tvals,xvals,yvals

def graph_arrows(x,y):
    """
        takes a list of x values and list of y values of the same size.
        plots arrows pointing in the direction of time of the solution.
    """
    if len(x) >5 and len(y) > 5:
        try:
            plt.annotate("",
                    xy=(x[len(x)/5 + 1], y[len(y)/5 + 1]), xycoords='data',
                    xytext=(x[len(x)/5], y[len(y)/5]), textcoords='data',
                    arrowprops=dict(arrowstyle="->",
                                    connectionstyle="arc3"),
                    )
            plt.annotate("",
                    xy=(x[len(x)/2 + 1], y[len(y)/2 + 1]), xycoords='data',
                    xytext=(x[len(x)/2], y[len(y)/2]), textcoords='data',
                    arrowprops=dict(arrowstyle="->",
                                    connectionstyle="arc3"),
                    )
        except None:
            print" ?"
    return True
    
    
def vector_field(xprime,yprime,gridmin=-3,gridmax=3,gridspace=0.1,scale=0.03):
    xpoints = np.arange(gridmin,gridmax,gridspace)
    ypoints = xpoints

    def arrow(x,y,f,g, scale = scale):
        xhead = f(x,y) 
        yhead = g(x,y) 
        norm = np.sqrt((xhead - x)**2 + (yhead - y)**2)
        ax1.arrow(x,y,scale*(xhead)/norm,scale*(yhead)/norm,head_starts_at_zero=True)
        normalvector = [xhead/norm, yhead/norm]
        normedlength = np.hypot(normalvector[0], normalvector[1])
        return norm, normalvector, normedlength
    

    for i in range(len(xpoints)):    
        for j in range(len(ypoints)):    
            arrow(xpoints[i],ypoints[j],xprime,yprime)


ax1 = plt.subplot2grid((4,4),(0,0), colspan=4,rowspan=2)
ax2 = plt.subplot2grid((4,4),(2,0), colspan=2,rowspan=2)
ax3 = plt.subplot2grid((4,4),(2,2), colspan=2,rowspan=2)
h= 0.0001 # this is the step size for the ab4 method. too large = unstable shit solution. too small = too long to calculate
time = np.arange(-5,5,h) # range of time values to draw the solution through. doesn't matter that its negative.
# these are the initial values to start at around the plane. all combinations of them.
vector_field(xprime,yprime,gridmin=-5,gridmax=5,gridspace=0.2)
y0 = [0.25,1./np.sqrt(3),1,0,0.5,0.9]
x0 = [0]
ax1.scatter(0,1./np.sqrt(3))
xval = np.arange(-2,2,0.1)
def f(x): return np.sqrt(1./3. - x*x/3.)
ax1.plot(xval,f(xval))
for i in x0:
    for j in y0:
        t,x,y = ab4_2d(i,j,xprime,yprime,time)
        
        #plot1 = plt.subplot(311)
        ax1.plot(x,y,c='r')
        ax1.set_xlim(-5,5)
        ax1.set_ylim(-5,5)
        ax1.grid(True)
        #plot1.xaxis.tick_top()
        #plt.title('Phase plot in 2-D\n')
        #plt.grid(True)
        #plt.arrow(x[100],y[100],x[101],y[101],head_width=0.5,head_length=0.5)
        #graph_arrows(x,y)

        #plot2 = plt.subplot(312)
        ax2.plot(t,x, c='b')
        ax2.set_xlim(-5,5)
        ax2.set_ylim(-5,5)
        #plt.setp(plot2.get_xticklabels(),visible=False)
        #plt.title("x(t)")
        #plt.grid(True)
   
        #plot3 = plt.subplot(313, sharex=plot2)
        ax3.plot(t,y,c='b')
        ax3.set_xlim(-5,5)
        ax3.set_ylim(-5,5)
        #plt.title("y(t)")
        #plt.grid(True)
