import numpy as np
from matplotlib import pyplot as plt

e, t = np.loadtxt('h2o.txt').transpose()
fontsize = 12
plt.ylabel('relative permittivity $\\epsilon_r$',fontsize=fontsize)
plt.xlabel('temperature ($^o$C)',fontsize=fontsize)
plt.title('relative permittivity of water vs. temperature',fontsize=1.3*fontsize)
plt.plot(t,e,'.',color='k')
plt.show()
