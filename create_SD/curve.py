from scipy import optimize
import sys
import math

x=float(sys.argv[1])
c_inf=float(sys.argv[2])
a=float(sys.argv[3])
beta_min=float(sys.argv[4])
print(beta_min - 1/a * math.log(1- x/c_inf))



