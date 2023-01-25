import math

mn = 0.94
outfile = open("pp_elastic_s.dat","w")
for str in open("pp_elastic_plab.dat"):
 cols = str.split()
 try:
  plab = float(cols[0])
  print >>outfile, 2.0*mn*(mn+math.sqrt(plab**2 + mn**2)), "   ", cols[1]
 except ValueError:
  print "Not a float:", cols[0]
