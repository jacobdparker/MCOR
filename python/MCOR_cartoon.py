import numpy as np
import matplotlib.pyplot as plt

#gaussian parameters for generating plots
sigma = 10
lag = -20 #negative lag would be a blue shifted object, positive red
range = 50
#Define a gaussian function internally that takes three arguments

def gaussian(x,lag,sigma):
    return np.exp(-(x - lag) ** 2 / (2 * sigma))


# Set up zero order object in black
x = np.linspace(-range,range,range*100)


plt.subplot(311)
y = gaussian(x,0,sigma)
plt.plot(x,y,color = 'black')



#plus order
y_plus = gaussian(x,lag,sigma)
plt.plot(x,y_plus, 'b')

#minus order
lag = -lag
y_minus = gaussian(x,lag,sigma)
plt.plot(x,y_minus, 'b--')


plt.title('Blue Shifted Object in Three MOSES Orders')
plt.legend(['Zero','Plus','Minus'])


#Sub-plot showing subtracted images

plt.subplot(312)
plt.title('Difference Images')
plt.plot(x,y_plus-y)
plt.plot(x,y_minus-y)
plt.legend(['Plus - Zero','Minus - Zero'])


#Sub-plot showing cross correlation of subtracted images
plt.subplot(313)
cc = np.correlate(y_plus-y,y_minus-y,mode='same')
plt.title('Cross Correlation of Plus - Zero and Minus - Zero')
plt.xlabel('Pixels')
plt.plot(x,cc,'r')

plt.tight_layout()
plt.savefig('MCOR_cartoon.pdf', bbox_inches='tight')
plt.show()

