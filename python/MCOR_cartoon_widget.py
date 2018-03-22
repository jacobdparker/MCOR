import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons


#gaussian parameters for generating initial plots
sigma = 1
shift = -1  #negative lag would be a blue shifted object, positive red
range = 50

#build x-axis
x = np.linspace(-range,range,range*100)

#Define a gaussian function internally that takes three arguments
def gaussian(x,shift,sigma):
    return np.exp(-(x - shift) ** 2 / (2 * sigma))

#begin creating a big figure
fig = plt.figure()

#make initial plots based on starting parameters

#Set up zero order object in black
plot1 = fig.add_subplot(411)
zero = gaussian(x,0,sigma)
[zero_plot] = plot1.plot(x,zero,color = 'black')

#plus order
plus = gaussian(x,shift,sigma)
[plus_plot]= plot1.plot(x,plus, 'b')

minus = gaussian(x,-shift,sigma)
[minus_plot]=plot1.plot(x,minus, 'b--')

plot1.legend(['Zero','Plus','Minus'])


#Sub-plot showing subtracted images

plot2 = fig.add_subplot(412)
#plot2.title('Difference Images')
[pmz] = plot2.plot(x,plus-zero)
[mmz] = plot2.plot(x,minus-zero)
plot2.legend(['Plus - Zero','Minus - Zero'])
plot2.set_ylim([-1.1,1.1])


#Sub-plot showing cross correlation of subtracted images
plot3 = fig.add_subplot(413)
cc = np.correlate(plus-zero,minus-zero,mode='same')
[cc_plot] = plot3.plot(x,cc,'r')
plot3.set_ylim([-50,50])

#plt.title('Cross Correlation of Plus - Zero and Minus - Zero')
#plot3.xlabel('Pixels')


# Add two sliders for tweaking the parameters

# Define an axes area and draw a slider in it
sigma_slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03])
sigma_slider = Slider(sigma_slider_ax, 'Sigma', 0.1, 100.0, valinit=sigma)

# Draw another slider
shift_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03])
shift_slider = Slider(shift_slider_ax, 'Shift', -30, 30, valinit=shift)

# Define an action for modifying the line when any slider's value changes
def sliders_on_changed(val):
    zero = gaussian(x, 0,sigma_slider.val)
    plus = gaussian(x,shift_slider.val,sigma_slider.val)
    minus = gaussian(x,-shift_slider.val,sigma_slider.val)
    zero_plot.set_ydata(zero)
    minus_plot.set_ydata(minus)
    plus_plot.set_ydata(plus)

    pmz.set_ydata(plus-zero)
    mmz.set_ydata(minus-zero)

    cc = np.correlate(plus - zero, minus - zero, mode='same')
    cc_plot.set_ydata(cc)
    fig.canvas.draw_idle()



#Call to function that controls slider change
sigma_slider.on_changed(sliders_on_changed)
shift_slider.on_changed(sliders_on_changed)


plt.tight_layout()
#plt.savefig('MCOR_cartoon.pdf', bbox_inches='tight')
plt.show()

