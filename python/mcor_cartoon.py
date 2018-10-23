import numpy as np
import matplotlib.pyplot as plt


#set plot size
plt.rcParams.update({'font.size': 14})


#build x-axis
rnge = 500
resolution = 1
x = np.linspace(-rnge,rnge,rnge*resolution+1)

#start by building and plotting three orders: plus, zero, minus
fig1 = plt.figure()
plot_zero = fig1.add_subplot(111)
fig2 = plt.figure()
plot_plus = fig2.add_subplot(111)
fig3 = plt.figure()
plot_minus = fig3.add_subplot(111)

#Define a gaussian function internally that takes three arguments
def gaussian(x,shift,sigma,int,center):
    return int*np.exp(-(x -center - shift) ** 2 / (2 * sigma))

#build red noise background for all orders
mu = 0
sigma = .1
white_noise = np.random.normal(mu,sigma,np.size(x))
white_noise[0] = 4  #set background noise value
red = np.cumsum(white_noise)

#start each order with a red noise background (boring, stationary He II)
#and add He II features to background
zero = red + gaussian(x, 0, 200, 10, 0) + gaussian(x, 0, 100, 10, 300)
plus = np.copy(zero)
minus = np.copy(zero)

#add He II to plots
plot_zero.plot(x, zero,color = 'b',label = 'He II')
plot_plus.plot(x, plus,color = 'b',label = 'He II')
plot_minus.plot(x,minus,color = 'b',label = 'He II')


#build arrays of parameters for spectral contam features
sigma = np.array([200., 200., 200.,200.])
shift = np.array([-17, 200.,150.,170.])
int = np.array([2,.8,.5,1])
center = np.array([0,0,0,0])
color = ['red','orange','green','purple']
labels = ['Si XI','Contam','Contam','Contam']

#initialize total signal for each channel
t_zero = np.copy(zero)
t_plus = np.copy(zero)
t_minus = np.copy(zero)

#add spectral contam to plots


for i in range(len(sigma)):
    zero = gaussian(x, 0, sigma[i], int[i],center[i])
    plus = gaussian(x, shift[i], sigma[i], int[i],center[i])
    minus = gaussian(x, -shift[i], sigma[i], int[i],center[i])
    plot_zero.plot(x, zero,color = color[i],label = labels[i])
    plot_plus.plot(x, plus,color = color[i],label = labels[i])
    plot_minus.plot(x,minus,color = color[i],label = labels[i])
    t_zero += zero
    t_plus += plus
    t_minus += minus

plot_zero.plot(x,t_zero,color = 'black',label='Total Signal')
plot_plus.plot(x,t_plus,color = 'black',label='Total Signal')
plot_minus.plot(x,t_minus,color='black',label='Total Signal')

#build cc x-axis for plotting
x_cc = np.linspace(-(np.size(x)*2-1),np.size(x)*2-1,np.size(x)*resolution*2-1)

#cross correlate each channel
cc_pz = np.correlate(t_plus,t_zero, mode='full')
cc_mz = np.correlate(t_minus,t_zero, mode='full')
fig4 = plt.figure()
cc_pz_plot = fig4.add_subplot(111)
fig5 = plt.figure()
cc_mz_plot = fig5.add_subplot(111)

#... and plot each cross_correlation
cc_pz_plot.plot(x_cc,cc_pz,color='black')
cc_mz_plot.plot(x_cc,cc_mz,color='black')

#plot subtracted images
fig6 = plt.figure()
pz_plot = fig6.add_subplot(111)
fig7 = plt.figure()
mz_plot = fig7.add_subplot(111)

pz = t_plus - t_zero
mz = t_minus - t_zero

pz_plot.plot(x,pz,'black')
mz_plot.plot(x,mz,'black')

#... and one last plot of this cross correlation

fig8 = plt.figure()
cc_pzmz_plot = fig8.add_subplot(111)

cc_pzmz = np.correlate(pz,mz, mode='full')
cc_pzmz_plot.plot(x_cc,cc_pzmz,'black')


#add a bunch of shit to the figures
plot_zero.set_title('M=0 order')
plot_plus.set_title('M=1 order')
plot_minus.set_title('M=-1 order')
pz_plot.set_title('M=1 minus M=0 (PZ)')
mz_plot.set_title('M=-1 minus M=0 (MZ)')
cc_pzmz_plot.set_title('PZ Cross-Correlated w/ MZ')
cc_pz_plot.set_title('Cross-Correlation of M=1 and M=0')
cc_mz_plot.set_title('Cross-Correlation of M=-1 and M=0')


plot_zero.set_xlabel('x (pixels)')
plot_plus.set_xlabel('x (pixels)')
plot_minus.set_xlabel('x (pixels)')

# plot_zero.set_ylabel('Intensity')
# plot_plus.set_ylabel('Intensity')
# plot_minus.set_ylabel('Signal')

plot_zero.legend()
plot_minus.legend()
plot_plus.legend()

pz_plot.set_xlabel('x (pixels)')
mz_plot.set_xlabel('x (pixels)')

cc_mz_plot.set_xlabel('Lag (pixels)')
cc_mz_plot.set_ylabel('Cross-Correlation')
cc_pz_plot.set_xlabel('Lag (pixels)')
cc_pz_plot.set_ylabel('Cross-Correlation')
cc_pzmz_plot.set_xlabel('Lag (pixels)')
cc_pzmz_plot.set_ylabel('Cross-Correlation')

#label all the figures
fig1.text(0,.9,'b)')
fig2.text(0,.9,'a)')
fig3.text(0,.9,'c)')
fig4.text(0,.9,'d)')
fig5.text(0,.9,'e)')
fig6.text(0,.9,'d)')
fig7.text(0,.9,'e)')
fig8.text(0,.9,'f)')

#save the figures!!!
plt.tight_layout()
fig1.savefig('fig1.pdf', bbox_inches='tight')
fig2.savefig('fig2.pdf', bbox_inches='tight')
fig3.savefig('fig3.pdf', bbox_inches='tight')
fig4.savefig('fig4.pdf', bbox_inches='tight')
fig5.savefig('fig5.pdf', bbox_inches='tight')
fig6.savefig('fig6.pdf', bbox_inches='tight')
fig7.savefig('fig7.pdf', bbox_inches='tight')
fig8.savefig('fig8.pdf', bbox_inches='tight')



