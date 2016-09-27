#!/usr/bin/env python
# encoding: utf-8
"""
chained_ops.py

Created by Matthew Loper on 2012-10-19.
Copyright (c) 2012 MPI. All rights reserved.
"""

import chained
import math
import scipy.sparse as sp
from body.matlab import row,col
import numpy as np

class Subtract (chained.Chained):
    def __init__(self, x1, x2):
        self.x1, self.x2 = self.differentiable_wrt(x1, x2)

    def compute_r(self):
        return self.x1.r - self.x2.r

    def compute_dr_wrt(self, wrt):
        x1 = self.x1.r
        if wrt == self.x1 and wrt == self.x2:
            return sp.csc_matrix((len(x1), len(x1)))
        elif wrt == self.x1:
            return sp.eye(len(x1), len(x1), format='csc')
        elif wrt == self.x2:
            return -sp.eye(len(x1), len(x1), format='csc')
        else:
            raise Exception('unknown object')

    def __repr__(self) :
        if hasattr(self, '_label') :
            return '<%s>' % self._label 
        return self.label_for_arg(self.x1) + '-' + self.label_for_arg(self.x2)

class MultiAnd(chained.Chained):
    def __init__(self, *x_args):
        args = []
        for x in x_args :
            args += x.concatenated_args if hasattr(x, 'concatenated_args') else [x]
        self.x_args = self.differentiable_wrt(*args)
    
    @property
    def concatenated_args(self):
        return self.x_args

    def compute_r(self):
        return np.concatenate([x.r.flatten() for x in self.x_args])
        
    def wrt_obj_overide(self, wrt_obj):
        wrt_id = id(wrt_obj)
        self.xhs = [len(x.r) for x in self.x_args]

        jacobians = []
        num_columns = None
        for i, x in enumerate(self.x_args) :
            jacobians.append(x.dr_wrt(wrt_obj) if (wrt_id in x.objids) else None)
            if jacobians[-1] != None : num_columns = jacobians[-1].shape[1]

        if not num_columns : raise Exception('unknown object')
        rows_to_stack = []
        for i, jac in enumerate(jacobians) :
            if jac != None :
                rows_to_stack.append(jac)
            elif self.xhs[i] :
                rows_to_stack.append(sp.csc_matrix((self.xhs[i], num_columns)))
        if np.sum([sp.issparse(r) for r in rows_to_stack]) :
            return sp.vstack(rows_to_stack).tocsc() if (len(rows_to_stack) > 1) else rows_to_stack[0]
        else :
            return np.vstack(rows_to_stack) if (len(rows_to_stack) > 1) else rows_to_stack[0]

    def __repr__(self) :
        return " & ".join([self.label_for_arg(x) for x in self.x_args])


class And_(chained.Chained):
    def __init__(self, x1, x2):
        self.x1, self.x2 = self.differentiable_wrt(x1, x2)
    
    def compute_r(self):
        return np.concatenate((self.x1.r.flatten(), self.x2.r.flatten()))
        
    def compute_dr_wrt(self, wrt):
        x1h = len(self.x1.r)
        x2h = len(self.x2.r)
        
        if wrt == self.x1 and x1h > 0:    
            return lambda x : sp.vstack((x, sp.csc_matrix((x2h, x1h))))
        elif wrt == self.x2 and x2h > 0:
            return lambda x : sp.vstack((sp.csc_matrix((x1h, x2h)), x))
    
    def __repr__(self) :
        return "%s & %s" % (self.label_for_arg(self.x1), self.label_for_arg(self.x2))
    
    
class And(chained.Chained):
    def __init__(self, x1, x2):
        self.x1, self.x2 = self.differentiable_wrt(x1, x2)
    
    @property
    def concatenated_args(self):
        return [self.x1, self.x2]

    def compute_r(self):
        return np.concatenate((self.x1.r.flatten(), self.x2.r.flatten()))
        
    def wrt_obj_overide(self, wrt_obj):
        wrt_id = id(wrt_obj)
        self.x_args = [self.x1, self.x2]
        self.xhs = [len(x.r) for x in self.x_args]

        jacobians = []
        num_columns = None
        for i, x in enumerate(self.x_args) :
            jacobians.append(x.dr_wrt(wrt_obj) if (wrt_id in x.objids) else None)
            if jacobians[-1] != None : num_columns = jacobians[-1].shape[1]

        if not num_columns : raise Exception('unknown object')
        rows_to_stack = []
        for i, jac in enumerate(jacobians) :
            if jac != None :
                rows_to_stack.append(jac)
            elif self.xhs[i] :
                rows_to_stack.append(sp.csc_matrix((self.xhs[i], num_columns)))

        return sp.vstack(rows_to_stack).tocsc() if (len(rows_to_stack) > 1) else rows_to_stack[0]

    def compute_dr_wrt(self, wrt):
        x1h = len(self.x1.r)
        x2h = len(self.x2.r)

        mtxs = []
        if wrt == self.x1 and x1h > 0:
            top = sp.eye(x1h, x1h)
            if x2h > 0:
                btm = sp.csc_matrix((x2h, x1h))
                mtxs.append(sp.vstack((top, btm)))
            else:
                mtxs.append(top)
        if wrt == self.x2 and x2h > 0:
            btm = sp.eye(x2h, x2h)
            if x1h > 0:
                top = sp.csc_matrix((x1h, x2h))
                mtxs.append(sp.vstack((top, btm)))
            else:
                mtxs.append(btm)

        result = np.sum(mtxs)
        return result.tocsc()
    
    def __repr__(self) :
        return "%s & %s" % (self.label_for_arg(self.x1), self.label_for_arg(self.x2))

        

class Divide (chained.Chained):
    def __init__(self, x1, x2):
        self.x1, self.x2 = self.differentiable_wrt(x1, x2)

    def compute_r(self):
        return self.x1.r / self.x2.r

    def compute_dr_wrt(self, wrt):
        result = []
        x1, x2 = (self.x1.r, self.x2.r)

        # inputs should be of equal size, or one/both of them could be scalars
        assert((x1.size == x2.size) or (x1.size==1) or (x2.size==1))

        if wrt == self.x1:
            if x2.size < x1.size:
                x2 = np.ones((x1.shape)) * x2
            result.append(sp.spdiags(1. / x2.flatten(), [0], x2.size, x2.size, format='csc'))
        if wrt == self.x2:
            if x1.size < x2.size:
                x1 = np.ones((x2.shape)) * x1
            result.append(sp.spdiags(-x1.flatten() / (x2*x2), [0], x1.size, x1.size, format='csc'))
        if not result:
            raise Exception('Cant compute dr wrt that variable')
            
            
   #   result = []
   #   x1, x2 = (self.x1.r, self.x2.r)
   #   if wrt == self.x1:
   #       result.append(sp.spdiags(1. / x2, [0], x2.size, x2.size))
   #   if wrt == self.x2:
   #       result.append(sp.spdiags(-x1 / (x2*x2), [0], x2.size, x2.size))
   #   if not result:
   #       raise Exception('Cant compute dr wrt that variable')
   
        if len(result)==1:
            return result[0]
        else:
            return np.sum(result).tocsc()
            
    # def compute_dr_wrt(self, wrt):
    #     if wrt == self.x1 and wrt == self.x2:
    #         return sp.csc_matrix((len(self.x1), len(self.x1)))
    #     elif wrt == self.x1: # should be 1 / x2
    #         return sp.spdiags(1. / self.x2, [0], self.x2.size, self.x2.size)
    #     elif wrt == self.x2: # should be -x1 / (x2 ** 2)
    #         return sp.spdiags(-self.x1 / (self.x2*self.x2), [0], self.x2.size, self.x2.size)
    #     else:
    #         raise Exception('unknown object')
    def __repr__(self) :
        if hasattr(self, '_label') :
            return '<%s>' % self._label 
        
        return "%s / %s" % (self.label_for_arg(self.x1), self.label_for_arg(self.x2))

            
class Multiply (chained.Chained):
    def __init__(self, x1, x2):
        self.x1, self.x2 = self.differentiable_wrt(x1, x2)

    def compute_r(self):
        return self.x1.r * self.x2.r

    def compute_dr_wrt(self, wrt):
        result = []
        x1, x2 = (self.x1.r, self.x2.r)

        # inputs should be of equal size, or one/both of them could be scalars
        assert((x1.size == x2.size) or (x1.size==1) or (x2.size==1))

        if wrt == self.x1:
            if x2.size < x1.size:
                x2 = np.ones((x1.shape)) * x2
            result.append(sp.spdiags(x2.flatten(), [0], x2.size, x2.size, format='csc'))
        if wrt == self.x2:
            if x1.size < x2.size:
                x1 = np.ones((x2.shape)) * x1
            result.append(sp.spdiags(x1.flatten(), [0], x1.size, x1.size, format='csc'))
        if not result:
            raise Exception('Cant compute dr wrt that variable')

        if len(result)==1:
            return result[0]
        else:
            return np.sum(result).tocsc()
            
    def __repr__(self) :
        if hasattr(self, '_label') :
            return '<%s>' % self._label 
        
        return "%s * %s" % (self.label_for_arg(self.x1), self.label_for_arg(self.x2))


class PowerOf (chained.Chained):
    """Given vector \f$x\f$, computes \f$x^2\f$ and \f$\frac{dx^2}{x}\f$"""
    def __init__(self, x, pow):
        self.x, self.pow = self.differentiable_wrt(x, pow)

    def compute_r(self):
        x, pow = (self.x.r, self.pow.r)
        epsv = 0 if math.fabs(pow) >= 1 else 4.7020e-38
        
        return col((x+epsv)**pow).flatten()

    def compute_dr_wrt(self, obj):
        return self.compute_dr_wrt(self.x)

    def compute_dr_wrt(self, wrt):
        x, pow = (self.x.r, self.pow.r)
        epsv = 0 if math.fabs(pow) >= 1 else 4.7020e-38

        result = []
        if wrt == self.x:
            result.append(sp.spdiags(pow*((x.flatten()+epsv)**(pow-1)), [0], x.size, x.size, format='csc'))
        if wrt == self.pow:
            result.append(sp.csc_matrix(col(np.array(np.log(x.flatten()) * (x.flatten()+epsv)**(pow)))))
        if not result:
            raise Exception('Cant compute dr wrt that variable')

        if len(result)==1:
            return result[0]
        else:
            return np.sum(result).tocsc()

    def __repr__(self) :
        if hasattr(self, '_label') :
            return '<%s>' % self._label 
        
        return "%s ** %s" % (self.label_for_arg(self.x), self.label_for_arg(self.pow))


class Sum(chained.Chained):
    def __init__(self, x1, x2):
        self.x1, self.x2 = self.differentiable_wrt(x1, x2)

    def compute_r(self):
        #print 'sum r entered:' + str(self.x1.r + self.x2.r)         
        return self.x1.r + self.x2.r

    def compute_dr_wrt(self, wrt):
        #print 'sum dr entered'
        x1 = self.x1.r
        if wrt == self.x1 and wrt == self.x2:
            return 2 * sp.eye(len(x1), len(x1), format='csc')
        elif wrt == self.x1 or wrt == self.x2:
            return sp.eye(len(x1), len(x1), format='csc')
        else:
            raise Exception('Cant compute dr wrt that variable')

    def __repr__(self) :
        if hasattr(self, '_label') :
            return '<%s>' % self._label 
        
        return self.label_for_arg(self.x1) + '+' + self.label_for_arg(self.x2)


def main():
    pass


if __name__ == '__main__':
    main()

