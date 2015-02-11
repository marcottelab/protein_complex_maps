from scipy.stats import gaussian_kde
import numpy as np

#kdrew: copied from http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html
def violin_plot(ax,data,pos, bp=False):
	'''
	create violin plots on an axis
	'''
	dist = max(pos)-min(pos)
	w = min(0.15*max(dist,1.0),0.5)
	nPoints = sum([len(y) for y in data])
	for d,p in zip(data,pos):
		k = gaussian_kde(d) #calculates the kernel density
		m = k.dataset.min() #lower bound of violin
		M = k.dataset.max() #upper bound of violin
		x = np.arange(m,M,(M-m)/100.) # support for violin
		v = k.evaluate(x) #violin profile (density curve)
		#v = v/v.max()*w #scaling the violin to the available space
		v = v/v.max()*(1.0*len(d)/nPoints) #scaling the violin to be proportional to number of counts
		ax.fill_betweenx(x,p,v+p,facecolor='y',alpha=0.3)
		ax.fill_betweenx(x,p,-v+p,facecolor='y',alpha=0.3)
	if bp:
		ax.boxplot(data,notch=1,positions=pos,vert=1)

    #pos = range(5)
    #data = [normal(size=100) for i in pos]
    #fig=plt.figure()
    #ax = fig.add_subplot(111)
    #violin_plot(ax,data,pos,bp=1)
    #plt.show()

