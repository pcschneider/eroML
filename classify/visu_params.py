import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

Blues = plt.get_cmap('Blues')
Reds = plt.get_cmap('Reds')



dd = np.genfromtxt("params.txt", unpack=False)

#def draw_triangle():
print(dd[0], dd[0])
dy = np.log10(dd[1,1]) - np.log10(dd[1,0])
print("dy: ",dy)
x = np.sort(np.log10(np.unique(dd[::,0])))
dx = x[1]-x[0]

l4max = np.max(dd[::,4])
l5max = np.max(dd[::,5])

for line in dd:
    #print(line)
    x = np.log10(line[0])
    y = np.log10(line[1])
    print(x,y)
    c1 = Blues(line[4]/l4max) # stars as other
    c2 = Reds(line[5]/l5max) # others as stars
    
    t1 = plt.Polygon(np.transpose([[x-dx/2,x-dx/2, x+dx/2], [y-dy/2,y+dy/2,y+dy/2]]), color=c1)
    plt.gca().add_patch(t1)
    
    t2 = plt.Polygon(np.transpose([[x-dx/2,x+dx/2, x+dx/2], [y-dy/2,y-dy/2,y+dy/2]]), color=c2)
    plt.annotate("%i" % line[4], xy=(x-dx/4,y+dy/4))
    plt.annotate("%i" % line[5], xy=(x,y-dy/4))
    plt.gca().add_patch(t2)
    
plt.xlabel("C")
plt.ylabel("class=0 weight")
plt.show()
