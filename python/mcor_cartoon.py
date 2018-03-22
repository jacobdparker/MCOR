import numpy as np
import matplotlib.pyplot as plt
import code

#build x-axis
rnge = 500
x = np.linspace(-rnge,rnge,rnge*100)

#Define a gaussian function internally that takes three arguments
def gaussian(x,shift,sigma,int,center):
    return int*np.exp(-(x -center - shift) ** 2 / (2 * sigma))

#build arrays of parameters for objects being plotted
sigma = np.array([200., 200., 200., 200.,75.,75., 50.,20,10,5,10])
shift = np.array([0, -17.,300.,320.,0,-17,0,0,0,0,0])
int = np.array([10,4,.8,.5,10,4,10,1,2,1,3])
center = np.array([0,0,0,0,120,120,-150,300,-300,200,170])
color = ['blue','red','orange','green','blue','red','blue','b','b','b','b']



#start by building and plotting three orders: plus, zero, minus
fig1 = plt.figure()
plot_zero = fig1.add_subplot(111)
fig2 = plt.figure()
plot_plus = fig2.add_subplot(111)
fig3 = plt.figure()
plot_minus = fig3.add_subplot(111)

#keep track of total symbols
t_zero = np.zeros_like(x)
t_plus = np.zeros_like(x)
t_minus = np.zeros_like(x)


for i in range(len(sigma)):
    zero = gaussian(x, shift[0], sigma[i], int[i],center[i])
    plus = gaussian(x, shift[i], sigma[i], int[i],center[i])
    minus = gaussian(x, -shift[i], sigma[i], int[i],center[i])
    plot_zero.plot(x, zero,color = color[i])
    plot_plus.plot(x, plus,color = color[i])
    plot_minus.plot(x,minus,color = color[i])
    t_zero += zero
    t_plus += plus
    t_minus += minus

plot_zero.plot(x,t_zero,color = 'black')
plot_plus.plot(x,t_plus,color = 'black')
plot_minus.plot(x,t_minus,color='black')


#cross correlate each channel
cc_pz = np.correlate(t_plus,t_zero, mode='same')
cc_mz = np.correlate(t_minus,t_zero, mode='same')

fig4 = plt.figure()
cc_pz_plot = fig4.add_subplot(111)
fig5 = plt.figure()
cc_mz_plot = fig5.add_subplot(111)

cc_pz_plot.plot(x,cc_pz,color='black')
cc_mz_plot.plot(x,cc_mz,color='black')

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

cc_pzmz = np.correlate(pz,mz, mode='same')
cc_pzmz_plot.plot(x,cc_pzmz,'black')

plt.tight_layout()
fig1.savefig('fig1.pdf', bbox_inches='tight')
fig2.savefig('fig2.pdf', bbox_inches='tight')
fig3.savefig('fig3.pdf', bbox_inches='tight')
fig4.savefig('fig4.pdf', bbox_inches='tight')
fig5.savefig('fig5.pdf', bbox_inches='tight')
fig6.savefig('fig6.pdf', bbox_inches='tight')
fig7.savefig('fig7.pdf', bbox_inches='tight')
fig8.savefig('fig8.pdf', bbox_inches='tight')




