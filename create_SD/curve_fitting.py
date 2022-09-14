from scipy import optimize
import sys
import numpy as np
def test_func(x,c_inf,a,beta_min):
    return c_inf * (1 - np.exp(-a*(x-beta_min))) 

def obtain_params(f,x,y):
    params, params_covariance = optimize.curve_fit(f, x, y,maxfev = 50000)
    return params

x=sys.argv[1].split(",")
y=sys.argv[2].split(",")
x=list(map(float,x))
y=list(map(float,y))
print(obtain_params(test_func,x,y))
