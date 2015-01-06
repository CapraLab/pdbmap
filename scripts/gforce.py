import scipy.constants.G
def gforce(coord1,coord2,m1=1,m2=1):
  x1,y1,z1 = coord1
  x2,y2,z2 = coord2
  r    = numpy.array([x2-x1,y2-y1,z2-z1])
  dist = np.sum([d**2 for d in r])
  fx,fy,fz = r * ((G * M * m) / dist)