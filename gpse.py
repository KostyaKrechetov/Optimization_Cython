from scipy.optimize.optimize import OptimizeResult as Result
import math

from matplotlib import pyplot
import pandas as pd
from hyperopt import fmin, rand, tpe, hp, STATUS_OK, Trials
import itertools
import random


import numpy as np
from sklearn.linear_model import LinearRegression
from math import log
import json


def gpse(fun, x0 , bounds, options, params):

    """
    Constrained optimizer for vector of N values. f --> min, c >= 0

    x0:           base control vector
    xl:           min values of control vector
    xu:           max values of control vector
    f_tolerance:     optimization stops if function change is less than <f_tolerance>. Must be positive.
    mode:            min | max

    return:          -
    """

    # Need to add maxiter, expected_parallel

    final_sten_size = options.get('final_sten_size',0.001)
    ftol            = options.get('ftol',           0.001)
    ftol_abs        = options.get('ftol_abs',       None)
    init_sten_size  = options.get('init_sten_size', 1.0)
   # use_constraints = options.get('use_constraints', False)
    Whether_resume  = options.get('Whether_resume',False)
    opt_object      = options.get('opt_object',None)
    maxiter         = options.get('maxiter',     100)
    
     
    cf = params.get('k' , 2)
    m = params.get('m')
    #sten_steps_num  = params.get('sten_steps_num',1)
    #sten_steps_step  = params.get('sten_steps_step',0)
    sten_steps = params.get('step_size_coeff',10)
    
    move_choose_size = int(m)+round(log(len(bounds)))

    #sten_steps = [1+x*sten_steps_step for x in range(1,sten_steps_num+1)]
    sten_steps = 0.1*sten_steps
    sten_steps = [sten_steps] * len(bounds)

    
    move_min_size = 2
    move_max_size = move_choose_size

    xl = [x[0] for x in bounds]
    xu = [x[1] for x in bounds]
    #x0 = list(x0)

    if ftol_abs is None and ftol is None:
        raise Exception('Either absolute or relative function tolerance must be specified!')

    opts = {
            'init_sten_size':           init_sten_size,
            'final_sten_size':          final_sten_size,
            'ftol':                     ftol,
            'ftol_abs':                 ftol_abs,
            'move_choose_size':         move_choose_size,
#            'use_constraints':          use_constraints,
    }

    



    
    
    
    
    if Whether_resume==False :
        f_best_list=[]
        u_calc = []
        f_calc = []
        got_base = False
        #c_calc = []
        x0_model = x0
        x0_math = norm(x0_model, xl, xu)
        x0 = x0_math
        u_cent = x0
        f_base = 0
        #c_base = 0
        u_calc.append(u_cent)
        
        #c_calc.append(c_cent)
    else:
        x0 = opt_object[0]
        
        x0_model = x0
        x0_math = norm(x0_model, xl, xu)
        x0 = x0_math
        u_cent = x0
        
        f_best_list=opt_object[4]
        u_calc = opt_object[2]
        f_calc = opt_object[3]
        init_sten_size = opt_object[5]
        f_base = opt_object[3][0]
        f_cent= opt_object[1]
        got_base = True
   
    

    #    print(' x0_model: %s' % (x0_model))
#    print(' x0_math: %s' % (x0_math))
#
#    print('Starting optimization, f --> min...')
#    print(' method: GPS enchanced')
#    print(' x0: %s' % (x0))
#    print(' xl: %s' % (xl))
#    print(' xu: %s' % (xu))
#
#    print('Options:')
#    print(json.dumps(opts, indent=2))


    sten_size = init_sten_size 
    q=0
    q1 = 0
    while sten_size >= final_sten_size and q1 < maxiter:
        q1+=1
        
        #u_cent = u_cent
        print('Stencil center: %s\nStencil size: %s' % (denorm(u_cent, xl, xu), sten_size))

        # Generate a pool of controls for one stencil
        xs = []

        for i in range(len(u_cent)):
            for direction in [1, -1]:
                x_new = [x for x in u_cent]
                x_new[i] += sten_size * sten_steps[i] * direction
                if x_new[i] > 1:
                    x_new[i] = 1
                if x_new[i] < 0:
                    x_new[i] = 0
                if x_new not in u_calc and x_new not in xs:
                    xs.append(x_new)

        if len(xs) == 0:
            sten_size /= cf
            continue

        if not got_base:
            xs.insert(0, x0_math)
            # xs.insert(0, [1] * len(x0_math))

        # Run exploration
        # ----------------------------------------------------------------------

#        print('x0_model: ', x0_model)
#        print('xl: ', xl)
#        print('xu: ', xu)
        
        results = math_fun(fun, xs, xl, xu)
        q+= len(xs)
        if not got_base:
            res = results[0][1]
            f_base = res
            #c_base = res.get('c', 0)

            if ftol is not None:
                if abs(f_base) < 1e-5:
                    print('*** Warning! You chose relative function tolerance, but function base value is too low: %s ***' % (f_base))
                if f_base == 0:
                    raise Exception('Can not use relative function tolerance with function base value = 0')

            f_cent = f_base
            #c_cent = c_base
            f_calc.append(f_cent)
            f_best_list.append(f_cent)
            q+= 1
            
            #c_calc.append(c_cent)
            got_base = True

        u_calc_sten = []
        f_calc_sten = []
        c_calc_sten = []
        for i, u in enumerate(xs):
            res = results[i][1]
            f = res
            #c = res.get('c', 0)

            u_calc_sten.append(u)
            f_calc_sten.append(f)
            if f >= f_best_list[len(f_best_list)-1] :
                f_best_list.append(f)
            else:
                f_best_list.append(f_best_list[len(f_best_list)-1])
            #c_calc_sten.append(c)

            u_calc.append(u)
            f_calc.append(f)
           
            #c_calc.append(c)

        # Find better variants
        betters = []
        for k in range(len(f_calc_sten)):
            fsten = f_calc_sten[k]
            #csten = c_calc_sten[k]

#            if use_constraints and csten < c_base:
#                continue

            if fsten > f_cent:
                betters.append((k, fsten))

        betters.sort(key=lambda b: b[1],reverse=True)
        
        if len(betters) == 0:
            print('No better point at this stencil. Reducing stencil size.')
            sten_size /= cf
            continue

        for_combs = betters[:move_choose_size]
        for_combs_is = [b[0] for b in for_combs]

        # Make combs
        combs = []
        for k in range(move_min_size, move_max_size+1):
            for subset in itertools.combinations(for_combs_is, k):
                combs.append(subset)

        # For each comb find new control
        xs_move = []
        for comb in combs:

            x_diff = [0] * len(x0)
            x_diff = np.array(x_diff)

            for k in comb:
                x_diff = x_diff + np.array(xs[k]) - np.array(u_cent)

            x_move = np.array(u_cent) + x_diff
            x_move = list(x_move)
            xs_move.append(x_move)

        # Run validation
        # ----------------------------------------------------------------------
     
        if len(xs_move) > 0:
            print('Validation of %s...' % len(xs_move))

        results = math_fun(fun, xs_move, xl, xu)
        q+= len(xs_move)
        for i, u in enumerate(xs_move):
            res = results[i][1]
            f = res
            #c = res.get('c', 0)

            u_calc_sten.append(u)
            f_calc_sten.append(f)
            if f >= f_best_list[len(f_best_list)-1] :
                f_best_list.append(f)
            else:
                f_best_list.append(f_best_list[len(f_best_list)-1])
            #c_calc_sten.append(c)

            u_calc.append(u)
            f_calc.append(f)
            #c_calc.append(c)

        for k in range(len(f_calc_sten)):
            fsten = f_calc_sten[k]
            #csten = c_calc_sten[k]

#            if use_constraints and csten < c_base:
#                continue

            if fsten > f_cent:
                betters.append((k, fsten))

        betters.sort(key=lambda b: b[1],reverse=True)

        # Choose best point
        f_best_sten = betters[0][1]
        if not f_best_sten > f_cent:
            print('No better point at this step.\n stencil center: %s\n current best: %s' % (f_best_sten, f_cent))
            sten_size /= cf

        i_best_sten = f_calc_sten.index(f_best_sten)
        u_best_sten = u_calc_sten[i_best_sten]
        #c_best_sten = c_calc_sten[i_best_sten]

        if ftol_abs is not None:
            if abs(f_best_sten - f_cent) <= abs(ftol_abs):

                # Objective improvement less than tolerance
                if sten_size < final_sten_size:
                    break
                else:
                    sten_size /= cf
                    ftol_abs /= cf

            # Objective improvement is relatively good
            if abs(f_best_sten - f_cent) > 2*abs(ftol_abs) and ftol_abs > 0:
                if sten_size <= 0.25:
                    sten_size *= cf
                    ftol_abs *= cf

        elif ftol is not None:
            if abs((f_best_sten - f_cent) / f_base) <= abs(ftol):

                # Objective improvement less than tolerance
                if sten_size < final_sten_size:
                    break
                else:
                    sten_size /= cf

        else:
            raise Exception('Either absolute or relative function tolerance must be specified!')
        
                
        # Move to new stencil center
        
        u_cent = u_best_sten
        f_cent = f_best_sten
        
        #c_cent = c_best_sten
    ##


    #f_best = min(f_calc)
    #i_best = f_calc.index(f_best)

    f_best = f_calc[0]
    for i in range(len(f_calc)):
        f = f_calc[i]
        #c = c_calc[i]

#        if use_constraints and c < c_base:
#            continue

        if f >= f_best:
            i_best = i
            f_best = f

    f_best = f_calc[i_best]
    u_best = u_calc[i_best]
    #c_best = c_calc[i_best]

    print('===============================================================')
    print('Optimization has finished!')
    print(' function evaluations: %s' %(len(f_calc)))
    print(' optimal function evaluation: %s' %(i_best))
   # print(' optimal controls: %s' %(denorm(u_best, xl, xu)))
    print(' optimal objective value: %s' %(f_best))
    print(' step size : %s' %(sten_size))
    print(' iteration number : %s' %(q1))
   # print(' optimal constraint value: %s' %(c_best))
    print('===============================================================')

    result = Result(fun=f_best, nit=None, nfev=len(f_calc),
                     status=None, success=None, message=None,
                     x=denorm(u_best, xl, xu))
    u_best=denorm(u_best, xl, xu)
    return u_best,f_best,u_calc,f_calc,f_best_list,sten_size,q,q1
###


def norm(x_real, xl, xu):

    if len(x_real) != len(xl) or len(x_real) != len(xu):
        raise Exception('Error! Inconsistent data. Lens of: x_real=%s, xl=%s, xu=%s' % (
            len(x_real), len(xl), len(xu)))

    x_math = []
    for i in range(len(x_real)):
        if float(xu[i]) - float(xl[i]) != 0:
            value = (float(x_real[i]) - float(xl[i])) / (float(xu[i]) - float(xl[i]))
        else:
            value = 0
        if value < 0 or value > 1:
            print('*** Warning! x_real[%s] is out of given ranges! ***\n x_real:%s\n xl:  %s\n xu:  %s' % (
            i, x_real, xl, xu))
            if value < 0:
                value = 0
            if value > 1:
                value = 1
        x_math.append(value)

    return x_math
##


def denorm(x_math, xl, xu):

    if len(x_math) != len(xl) or len(x_math) != len(xu):
        raise Exception('Error! Inconsistent data. Lens of: x_math=%s, xl=%s, xu=%s' % (
            len(x_math), len(xl), len(xu)))

    x_real = []
    for i in range(len(x_math)):
        value = xl[i] + x_math[i] * (xu[i] - xl[i])
        x_real.append(value)

    return x_real
##


def math_fun(fun, xs, xl, xu):

    xs = [denorm(x, xl, xu) for x in xs]

    #for x in xs:
        #print(x)

    results = [(x,fun(x)) for x in xs]
    return results
##


#if __name__ == '__main__':
   # math_fun(None, [[0.125, 0.25, 0.5, 0.625]], [0, 0, -4, -4], [4]*4)







def rastrigin(X):
    return -(10*len(X)+sum([(x**2 - 10 * np.cos(2 * math.pi * x)) for x in X]))


def rosenbrock(X):   
    return -sum([ ((1-X[i])**2 + 100*(X[i+1]-X[i]**2)**2) for i in range(0,len(X)-1)])

def himmelblau(X):
    return -((X[0]**2 + X[1] -11)**2 + (X[0]+X[1]**2-7)**2)

def levi_N13(x):
    term1 = (np.sin(3*math.pi*x[0]))**2
    term2 = (x[0]-1)**2 * (1+(np.sin(3*math.pi*x[1]))**2)
    term3 = (x[1]-1)**2 * (1+(np.sin(2*math.pi*x[1]))**2)
    return  term1 + term2 + term3
def shekel(x):
    return -sum([1/(i+sum([(x[j]-j*i)**2 for j in range(0,len(x)-1)]))  for i in range(1,100)])
    
#def ackley (x):
    #return 











 
