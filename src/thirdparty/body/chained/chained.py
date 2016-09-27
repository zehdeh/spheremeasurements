__all__ = ['Chained', 'Cw', 'ChainedNoOp', 'cached_property', 'MatVecMult', 'VecsDotVec', 'Pssq', 'Rewgt', 'Prmse', 'PrintR', 'SumOf', 'RepSumOf', 'Select', 'Sample', 'Regularize', 'Expand']

from ..matlab.matlab import row, col
from ..misc.sparse_id import sparse_type
import types
import itertools
import time
from functools import wraps
from copy import deepcopy

import math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg

def cached_property(func):
    @wraps(func)
    def with_caching(self, *args, **kwargs):
        if not hasattr(self, '_cached_properties'):
            self._cached_properties = {}
        if func.__name__ not in self._cached_properties:
            self._cached_properties[func.__name__] = func(self, *args, **kwargs)
        return self._cached_properties[func.__name__]
    return property(with_caching)




class Chained(object):
    """Allows chaining of value/derivative computation.
    
    Use this class by inheriting from it and by defining
    the methods compute_r() and compute_dr_wrt(), which should
    compute the output vector and the jacobian respectively.
    
    Note that x may be a vector or a matrix, and that the result
    from compute_r() may be a vector or a matrix, but the result
    of compute_dr_wrt() must be a matrix (a jacobian) so that it 
    may be multiplied by others.
    """

    def __new__(cls, *args, **kwargs):
        """Intercepts and processes first non-self argument to init (ie "x" vector).
        
        If "x" has attributes "r" and "dr", we store "dr" for later (to multiply
        according to the chain rule) and pass "r" to __init__ of child class. 
        If "x" is just a vector, then we pass that directly to __init__ of child class.        
        """
        
        # TODO: Make work later?
        #
        # if callable(x):
        #     def f(y):
        #         s = Chained.__new__(cls, x.__call__(y), *args, **kwargs)
        #         s.__init__(*args, **kwargs)
        #         return s
        #     return lambda y : f(y)
        
        label = False
        if 'label' in kwargs.keys() :
            label = kwargs['label']
            del kwargs['label']

        # Actually construct an instance of the child class
        result = super(Chained, cls).__new__(cls)# , *args, **kwargs)
        result._cached_properties = {}
        result._drs_multiplied = {}
        result._compute_dr_wrt_cached = {}       
        result._nonuser_attributes = set(dir(result))      

        # These contain all the arguments to init. When "differentiable_wrt"
        # is called, those particular arguments will be removed from these lists.
        # This is done so that we can manage equality/redundancy operations
        result._static_args = args
        result._static_kwargs = kwargs

        if label :
            result._label = label
            
        
        return result


    @cached_property
    def r(self):
        """Return a function of the input vector x."""
        while hasattr(self, 'redirect'):
            self = self.redirect        
        return self.compute_r()

    def call_on_step(self, success):
        if hasattr(self, 'step'):
            raise Exception('step() callback has been renamed to on_step(), it takes success as a parameter, and may return True if optimization should be aborted.')
        if hasattr(self, 'on_step') and callable(self.on_step) and self.on_step(success) == True:
            return True

        abort = False
        if hasattr(self, '_all_args'):
            for a in self._all_args:
                if a.call_on_step(success):
                    abort = True
        return abort

    def find_names(self, obj):
        import gc, sys
        frame = sys._getframe()
        for frame in iter(lambda: frame.f_back, None):
            frame.f_locals
        result = []
        for referrer in gc.get_referrers(obj):
            if isinstance(referrer, dict):
                for k, v in referrer.iteritems():
                    if v is obj:
                        result.append(k)
        return result
        
    def show_tree(self, depth=0):
        print(' ' * depth * 3 + '---> ' + str(self.__class__.__name__) + ' (id=%d)' % (id(self)))
        for a in self._all_args:
            a.show_tree(depth+1)


    def reassign_branches(self, branches):
        any_reassigned = False
        
        old_arg_ids = [id(obj) for obj in self._all_args]
        for i in range(len(self._all_args)):            
            arg = self._all_args[i]
            if hasattr(arg, '_all_args'):
                if arg.reassign_branches(branches):
                    any_reassigned = True     
            if self._all_args[i] in branches and id(self._all_args[i]) != id(branches[self._all_args[i]]):
                self._all_args[i].redirect = branches[self._all_args[i]]
                self._all_args[i] = branches[self._all_args[i]]
            
        new_arg_ids = [id(obj) for obj in self._all_args]
        if new_arg_ids != old_arg_ids:
            any_reassigned = True
            #print 'REASSIGNMENT OCCURRED'
        
        return any_reassigned
        
    def remove_redundancies(self):
        # while reassignment is successful, keep
        # retraversing the tree
        while True:
            branches = self.get_branches_in_tree()        
            if not self.reassign_branches(branches):
                break

    def get_branches_in_tree(self):
        d =  { self : self }
        if hasattr(self, '_all_args'):
            for a in self._all_args:
                if isinstance(a, Chained):
                    d = dict(list(d.items()) + list(a.get_branches_in_tree().items()))
            
        return d
                
                
            
    def dr_wrt(self, wrt_obj):                
        while hasattr(self, 'redirect'):
            self = self.redirect        

        # If we have the derivative with respect to "wrt_obj" already, 
        # return it and we're done.
        wrt_id = id(wrt_obj)        
        if wrt_id in self._drs_multiplied.keys():    
            result = self._drs_multiplied[wrt_id]
            assert(len(result.shape)==2)
            assert(result.shape[0] == self.r.size)
            assert(result.shape[1] == (wrt_obj.r.size if hasattr(wrt_obj, 'r') else wrt_obj.size))
            return result
        
        if hasattr(self, 'wrt_obj_overide'):
            return self.wrt_obj_overide(wrt_obj)

        #thread = threading.Thread(target=get_file, args=(url,))
        # For each "free" arg we received in __init__ which used
        # "wrt_obj" as input, multiply our jacobian by one
        # earlier in the chain. Note that our jacobian may be
        # either computed or retrieved from cache.
        self._drs_multiplied[wrt_id] = []
        for arg in set(self._all_args):

            # if this argument isn't differentiable wrt this input, 
            # skip this argument
            if wrt_id not in arg.objids:
                continue

            arg_id = id(arg)

            # compute dr_wrt_r: our jacobian
            if arg_id not in self._compute_dr_wrt_cached.keys():
#                print 'computing derivative of %s wrt %s...' % (self, arg)
                tm = time.time()
                dr_wrt_r = self.compute_dr_wrt(arg)
#                print 'Derivative computed in %.3fs' % (time.time() - tm)
                if sp.issparse(dr_wrt_r) and self.__class__.__name__ != 'PassAlong':
                    assert sp.isspmatrix_csc(dr_wrt_r), "Sparse jacobian is %s not csc at %s!" % (sparse_type(dr_wrt_r), self.__class__.__name__)
                self._compute_dr_wrt_cached[arg_id] = dr_wrt_r
            else:
                dr_wrt_r = self._compute_dr_wrt_cached[arg_id]
                
            # do matrix multiplication of our jacobian by that of arg
            if isinstance(arg, Cw):
#                print 'chaining derivative of %s wrt %s...' % (self, arg)
                tm = time.time()
                self._drs_multiplied[wrt_id].append(dr_wrt_r)
#                print 'chain computed in %.3fs' % (time.time() - tm)
            else:                        
                prev_jac = arg.dr_wrt(wrt_obj)
#                print 'chaining derivative of %s wrt %s...' % (self, arg)
                tm = time.time()
                assert(len(prev_jac.shape)==2)
                
                if isinstance(dr_wrt_r, sp.linalg.interface.LinearOperator):
                    new_jac = dr_wrt_r.matmat(prev_jac)
                elif isinstance(dr_wrt_r, type(lambda:None)):
                    new_jac = dr_wrt_r(prev_jac)
                else:
                    # if dr_wrt_r isn't sparse, and prev_jac is sparse, "dot" won't work,
                    # so we have to play some transposition games
                    if (not sp.issparse(dr_wrt_r)) and (sp.issparse(prev_jac)):
                        new_jac = prev_jac.T.dot(dr_wrt_r.T).T
                    else:
                        # if dr_wrt_r contains many rows with all zeros, it makes sense to keep the jacobian sparse
                        if sp.issparse(dr_wrt_r) and (0.05 > float((dr_wrt_r.sum(axis=1) != 0).sum())/float(dr_wrt_r.shape[0])) :
                            new_jac = dr_wrt_r.dot(sp.csc_matrix(prev_jac))
                        else :
                            new_jac = dr_wrt_r.dot(prev_jac)
                assert(new_jac.shape[0]==self.r.size)
                assert(new_jac.shape[1]==(wrt_obj.r.size if hasattr(wrt_obj, 'r') else wrt_obj.size))
                self._drs_multiplied[wrt_id].append(new_jac)
#                print 'chain computed in %.3fs' % (time.time() - tm)
                
                    
        if len(self._drs_multiplied[wrt_id]) == 0:
            raise Exception("Can't find derivative wrt variable with id=%d (name is maybe %s from set [%s])" % (wrt_id, (self.find_names(wrt_obj))[-1], str(self.find_names(wrt_obj))))

        # FIXME
        #pre = deepcopy(self._drs_multiplied[wrt_id])
#        print 'adding chained derivative of %s...' % self
        tm = time.time()
        if len(self._drs_multiplied[wrt_id]) == 1:            
            self._drs_multiplied[wrt_id] = self._drs_multiplied[wrt_id][0]
        else:
            assert(isinstance(self._drs_multiplied[wrt_id], list))
            self._drs_multiplied[wrt_id] = reduce(lambda x, y: x+y, self._drs_multiplied[wrt_id])
            
            
        result = self._drs_multiplied[wrt_id]
        if len(result.shape) != 2:
            import pdb; pdb.set_trace()
        assert(len(result.shape)==2)
        assert(result.shape[0] == self.r.size)
        assert(result.shape[1] == (wrt_obj.r.size if hasattr(wrt_obj, 'r') else wrt_obj.size))
#        print 'sum computed in %.3fs' % (time.time() - tm)

#        if hasattr(self, 'wrt_obj_overide'):
#            assert ((result - self.wrt_obj_overide(wrt_obj)).sum() == 0)
#            print 'wrt_obj_overide checks out!'

        return result                   


    def differentiable_wrt(self, *args):
        if not isinstance(args, tuple):
            args = (args,)

        # Remove some items from "static_args" and "static_kwargs":
        # specifically, we don't want "diff_wrt" args in 
        # static_args and static_kwargs                
        old_arg_ids = [id(obj) for obj in args]
        self._static_args = [a for a in self._static_args if id(a) not in old_arg_ids]
        for k in self._static_kwargs.keys():
            if id(self._static_kwargs[k]) in old_arg_ids:
                del self._kwargs[k]
#                del self._static_kwargs[k]                

        old_args = args
        new_args = []
        self.objids = []
        
        # wrap scalars and numpy arrays in Cw objects
        for arg in old_args:
            self.objids.append(id(arg))
            if not isinstance(arg, Chained) and not isinstance(arg, Cw):
                arg = Cw(arg)
            self.objids += arg.objids
            new_args.append(arg)
                
        
                    
        self._all_args = new_args

        # if the user just wants one output, don't return it as a tuple
        if len(args) == 1:
            return self._all_args[0]
        else:
            return tuple(self._all_args)
    
    def set_label(self, label):
        self._label = label
        return self
    
    def label_for_arg(self, arg):
        if isinstance( arg , ( int, long, float, str ) ) :
            return str(arg)
        if (arg.__class__.__name__ != 'Cw') :
            return arg.__repr__() 

        if (arg.__class__.__name__ == 'Cw') and isinstance( arg.raw_x , ( int, long, float, str ) ) :
            return str(arg.raw_x)
        return ('<%s>' % self.find_names(arg)[0]) if self.find_names else '<?>'
    
    def __repr__(self) :
        if hasattr(self, '_label') :
            return '<%s>' % self._label 
        else :
            return self.__class__.__name__ + '(' + ', '.join([self.label_for_arg(arg) for arg in self._all_args]) + ')'
    
    def __add__(self, other):
        return chained_ops.Sum(self, other)
    def __sub__(self, other):
        return chained_ops.Subtract(self, other)
    def __pow__(self, other):
        return chained_ops.PowerOf(self, other)
    def __div__(self, other):
        return chained_ops.Divide(self, other)
    def __and__(self, other):
#        return chained_ops.And(self, other)
        return chained_ops.MultiAnd(self, other)
    def __mul__(self, other):
        return chained_ops.Multiply(self, other)
    
    def __rshift__ (self, reciever):
        if not reciever :
            return self
        if hasattr(reciever, '__call__'):
            reciever = reciever(self)
        if hasattr(reciever, 'on_compute_r'):
            return PassAlong(self, reciever)
        else:
            raise Exception('The second argument of chained class operator >> must be: an object that implemnts on_compute_r, a function that returns such an object, or None')



    def _key(self):
        k = ['static_args']
        k += [id(obj) for obj in self._static_args]
        k.append('static_kwarg_keys')
        k += [s for s in self._static_kwargs.keys()]
        k.append('static_kwarg_values')
        k += [id(obj) for obj in self._static_kwargs.keys()]
        k.append('class')
        k.append(str(self.__class__))
        if hasattr(self, '_all_args'):
            k += [a._key() for a in self._all_args]
        result = tuple(k)
        foo = hash(result) # for testing
        return result

    def __eq__(self, other):
        
        # comparing id's is faster, but may not always work.
        # some objects are equivalent but have different id's.
        if id(self) == id(other):
            return True
        return self._key() == other._key()

    def __hash__(self):
        return hash(self._key())
        
        
        
            
class Cw(object):
    def __init__(self, x):
        original_id = id(x)
        self.raw_x = x
        if np.isscalar(x):
            self.x = np.array(x)
        else:
            self.x = x
        self.objids = list(set([id(self), id(self.x), original_id]))

    def show_tree(self, depth=0):
        print(' ' * depth * 3 + '---> ' + str(self.__class__.__name__) + ' (id=%d)' % (id(self)))
        print(' ' * (depth+1) * 3 + '---> ' + str(self.x.__class__.__name__) + ' (id=%d)' % (id(self.x)))

    def set_label(self, label):
        self._label = label
        return self
    
    def call_on_step(self, success): return False

    def __add__(self, other):
        return chained_ops.Sum(self, other)
    def __sub__(self, other):
        return chained_ops.Subtract(self, other)
    def __mul__(self, other):
        return chained_ops.Multiply(self, other)
    def __pow__(self, other):
        return chained_ops.PowerOf(self, other)
    def __div__(self, other):
        return chained_ops.Divide(self, other)
    def __and__(self, other):
#        return chained_ops.And(self, other)
        return chained_ops.MultiAnd(self, other)
    def __rshift__ (self, reciever):
        if not reciever :
            return self
        if hasattr(reciever, '__call__'):
            reciever = reciever(self)
        if hasattr(reciever, 'on_compute_r'):
            return PassAlong(self, reciever)
        else:
            raise Exception('The second argument of chained class operator >> must be: an object that implemnts on_compute_r, a function that returns such an object, or None')

    def _key(self):
        return (id(self.x),)

    @property
    def r(self):
        return self.x

    def dr_wrt(self, obj):
        raise Exception('Shouldnt get here')
        

    
class Sample(Chained):
    def __init__(self, x, indices_to_keep, entries_per_index=1):
        self.x = self.differentiable_wrt(x)
        self.indices_to_keep = np.array([ [entries_per_index*index + i for i in range(entries_per_index)] for index in indices_to_keep ]).flatten()

    def compute_r(self):
        return self.x.r[self.indices_to_keep].copy()

    def compute_dr_wrt(self, obj):
#        diag = np.zeros(len(self.x.r))
#        diag[self.indices_to_keep] = 1
#        return sp.spdiags(diag, [0], diag.size, diag.size)
        
        IS = np.arange(len(self.indices_to_keep)).flatten()
        JS = self.indices_to_keep.flatten()
        ij = np.vstack((row(IS), row(JS)))
        data = np.ones(len(self.indices_to_keep))
        return sp.csc_matrix((data, ij), shape=(len(self.indices_to_keep), len( self.x.r )))

class MatVecMult(Chained):
    def __init__(self, mtx, vec):
        self.vec = self.differentiable_wrt(vec)
        self.mtx = mtx

    def compute_r(self):
        return self.mtx.dot(col(self.vec.r)).flatten()

    def compute_dr_wrt(self, obj):
        return sp.csc_matrix(self.mtx)

class VecsDotVec(Chained):
    def __init__(self, x, vec):
        self.x, self.vec = self.differentiable_wrt(x, vec)

    def compute_r(self):
        return np.array([self.vec.r.dot(v) for v in self.x.r.reshape(-1, len(self.vec.r.flatten()))]).flatten()

    def compute_dr_wrt(self, wrt):
        if wrt == self.x :
            IS = np.array([range(self.x.r.size/3)]*3).T.flatten()
            JS = np.arange(self.x.r.size).flatten()
            ij = np.vstack([row(IS), row(JS)])
            return sp.csc_matrix((np.hstack([self.vec.r.flatten()]*(self.x.r.size/3)), ij), shape=(self.x.r.size/3, self.x.r.size))
        else :
            return sp.csc_matrix(np.array([v for v in self.x.r.reshape(-1, len(self.vec.r.flatten()))]).flatten())
    
class VertexTransform(Chained):
    def __init__(self, verts, three_by_three_mat, translation=[0.0, 0.0, 0.0]):
        self.verts = self.differentiable_wrt(verts)
        self.mtx = np.matrix(three_by_three_mat)
        self.trans = np.matrix(translation)

    def compute_r(self):
        vertices = self.verts.r.reshape(-1,3)
        return np.array(self.mtx*vertices.T + self.trans).T.flatten()

    def compute_dr_wrt(self, obj):
        vertices = self.verts.r.reshape(-1,3)
        return sp.block_diag(itertools.repeat(self.mtx, vertices.shape[0]), format='csc')

class VertexProjection(Chained):
    def __init__(self, verts, camera_matrix, rotation=np.eye(3), translation=np.zeros((3,1)), return_xyz=True):
        self.verts = self.differentiable_wrt(verts)
        self.mtx = np.matrix(camera_matrix)*np.matrix(rotation)
        self.trans = np.matrix(camera_matrix)*np.matrix(translation)
        self.return_xyz = return_xyz

    def compute_r(self):
        vertices = self.verts.r.reshape(-1,3).T
        transformed_vertices = (self.mtx*vertices + self.trans)
        xy = np.array(transformed_vertices[:2]/transformed_vertices[2,:])
        return np.vstack([xy, np.zeros(vertices.shape[1])]).T.flatten() if self.return_xyz else xy.T.flatten()

    def compute_dr_wrt(self, obj):
        vertices = self.verts.r.reshape(-1,3)
        m = np.array(self.mtx)
        num_inputs = len(self.verts.r)
        outputs_per_inputs = 3 if self.return_xyz else 2

        data = []
        rows = []
        cols = []
        ind = np.arange(vertices.shape[0])
        for n in range(3) :
            data.extend(m[0,n]/vertices[:,2] - vertices[:,0]*m[2,n]/(vertices[:,2]**2))
            rows.extend(outputs_per_inputs*ind)
            cols.extend(3*ind + n)
        for n in range(3) :
            data.extend(m[1,n]/vertices[:,2] - vertices[:,1]*m[2,n]/(vertices[:,2]**2))
            rows.extend(outputs_per_inputs*ind + 1)
            cols.extend(3*ind + n)
        derivatives = sp.coo_matrix((data, (rows, cols)), shape=(outputs_per_inputs*num_inputs/3, num_inputs))
#        derivatives = sp.lil_matrix((outputs_per_inputs*num_inputs/3, num_inputs))
#        for i, v in enumerate(vertices) :
#            derivatives[outputs_per_inputs*i, 3*i:3*(i+1)] = m[0,:]/v[2] - v[0]*m[2,:]/(v[2]**2)
#            derivatives[outputs_per_inputs*i + 1, 3*i:3*(i+1)] = m[1,:]/v[2] - v[1]*m[2,:]/(v[2]**2)
        return sp.csc_matrix(derivatives)



class Select(Chained):
    def __init__(self, parms, which):
        self.parms = self.differentiable_wrt(parms)
        self.which = np.array(which)

    def compute_r(self):
        return self.parms.r[self.which].copy()

    def compute_dr_wrt(self, obj):
        IS = np.arange(len(self.which)).flatten()
        JS = self.which.flatten()
        ij = np.vstack((row(IS), row(JS)))
        data = np.ones(len(self.which))
        return sp.csc_matrix((data, ij), shape=(len(self.which), len(self.parms.r)))

class ChainedPassthrough(Chained):
    def __init__(self, x, *args, **kwargs):
        self.x = self.differentiable_wrt(x)

    def compute_r(self):
        if hasattr(self, 'on_compute_r') :
            self.on_compute_r()
        return self.x.r

    def compute_dr_wrt(self, obj):
#        return sp.LinearOperator(vecmat = lambda x : x, matmat = lambda x : x)
        if hasattr(self, 'on_compute_dr_wrt') :
            self.on_compute_dr_wrt(obj)
        return sp.eye(len(self.x.r), len(self.x.r), format='csc')


class ChainedNoOp(Chained):
    def __init__(self):
        self.objids = []    
        self._all_args = []

    def compute_r(self):
        return np.array([])

    def compute_dr_wrt(self, obj):
        return np.array([])


# The first argument of PassAlong should be Chained object, the second can be any object that
# impliments on_compute_r. PassAlong should typically be constructed via the __rshift__ operator ('>>')
# applied to the Chained object.
class PassAlong(ChainedPassthrough):
    def __init__(self, x, reciever):
        self.x = self.differentiable_wrt(x)
        self.reciever = reciever
        
    def on_compute_r(self):
        self.reciever.on_compute_r(self.x)
    
    def on_compute_dr_obj(self, obj) :
        if hasattr(self.reciever, 'on_compute_dr') :
            self.reciever.on_compute_dr_wrt(self.x, obj)
    
    def on_step(self, success):
        if hasattr(self.reciever, 'on_step') :
            self.reciever.on_step(success);
            

# This class is designed specifically to print chained objects.
# If call_on_input returns something that is not a chained object, and does not have an 'r' attribute, the return value will be printed directly
class PrintR():
    def __init__(self, call_on_input=None, debug_mode=False, to_stderr=False):
        self.call_on_input = call_on_input if call_on_input else (lambda x : x)
        self.debug_mode = debug_mode
        self.stuff_to_print = [];
        self.to_stderr = to_stderr
    def on_compute_r(self, chained_obj):
        object_to_print = self.call_on_input(chained_obj)
        quantity_to_print = object_to_print.r if (isinstance(object_to_print, Chained) or isinstance(object_to_print, Cw) or hasattr(object_to_print, 'r')) else object_to_print
                
        if len(quantity_to_print) == 1:
            quantity_to_print = '%.3e' % quantity_to_print

        self.stuff_to_print.append('%s: %s' % (chained_obj.__repr__(), quantity_to_print))
        if self.debug_mode :
            print('r for %s: %s' % (chained_obj.__repr__(), quantity_to_print))
        
    def on_step(self, success):
        if success and self.stuff_to_print:
            if self.to_stderr:
                from sys import stderr
                stderr.write(' \t'.join(self.stuff_to_print)+'\n')
            else:
                print(' \t'.join(self.stuff_to_print))
        self.stuff_to_print = [];
            
class SumOf(Chained):
    def __init__(self, x):
        self.x = self.differentiable_wrt(x)

    def compute_r(self):
        return np.array([np.sum(self.x.r)])

    def compute_dr_wrt(self, obj):
        return sp.csc_matrix(np.ones(len(self.x.r)))
#        return sp.csc_matrix(np.array([np.sum(self.x.dr_wrt(obj))]))

class RepSumOf(Chained):
    def __init__(self, x, rep):
        self.x = self.differentiable_wrt(x)
        self.rep = rep
    def compute_r(self):
        return np.sum(self.x.r)*np.ones(self.rep)

    def compute_dr_wrt(self, obj):
        return sp.csc_matrix(np.ones((len(self.x.r), self.rep)))


class Regularize(Chained):
    def __init__(self, x, mean=None):
        self.x = self.differentiable_wrt(x)
        self.mean = Cw(mean) if mean else None

    def compute_r(self):
        return self.x.r - self.mean.r if self.mean else self.x.r

    def compute_dr_wrt(self, obj):
        return self.x.dr_wrt(obj)

class Pssq (Chained):
    def __init__(self, result, resultname):        
        self.result = self.differentiable_wrt(result)
        self.resultname = resultname

    def compute_r(self):
        print ('%s ssq: %e' % (self.resultname, np.sum(self.result.r.flatten()**2.)))
        return self.result.r

    def compute_dr_wrt(self, obj):
        return sp.eye(len(self.result.r),len(self.result.r), format='csc')

class Prmse (Chained):
    def __init__(self, result, resultname):        
        self.result = self.differentiable_wrt(result)
        self.resultname = resultname

    def compute_r(self):
        print ('%s rmse: %e' % (self.resultname, np.sqrt(np.mean(self.result.r.flatten()**2.))))
        return self.result.r

    def compute_dr_wrt(self, obj):
        return sp.eye(len(self.result.r),len(self.result.r), format='csc')

class Rewgt(Chained):

    def __init__(self, x1, baseweight=1.):
        self.x1 = self.differentiable_wrt(x1)
        self.baseweight = baseweight

    def compute_r(self):
        return self.x1.r * self.baseweight / len(self.x1.r)

    def compute_dr_wrt(self, obj):
        return sp.eye(len(self.x1.r), len(self.x1.r), format='csc') * self.baseweight / len(self.x1.r)

class Expand(Chained):
    def __init__(self, subvec, i_subvec, fullvec):
        self.subvec = self.differentiable_wrt(subvec)
        self.i_subvec = list(i_subvec)
        self.fullvec = fullvec

    def compute_r(self):
        tmp = self.fullvec.copy()
        tmp[self.i_subvec] = self.subvec.r
        return tmp

    def compute_dr_wrt(self, obj):
        n_sub = len(self.i_subvec)
        return sp.csc_matrix((np.ones(n_sub), np.vstack((self.i_subvec, range(n_sub)))),
                              shape=(self.fullvec.size, self.subvec.r.size))
# This must be here (not at top of file) b/c
# of circular dependency issues.
import chained_ops.py

def main():
    import numpy as np

    # TEST SUM
    print("should print 2 x identity, and [2 4 6]")
    c1 = Cw(np.array([1,2,3]))
    c2 = c1 + c1
    print(c2.dr_wrt(c1).todense())
    print(c2.r)
    # 
    # print "should print identity twice, and [5 7 9]"
    # c1 = Cw(np.array([1,2,3]))
    # c2 = Cw(np.array([4,5,6]))
    # c3 = c1 + c2
    # print c3.dr_wrt(c1).todense()
    # print c3.dr_wrt(c2).todense()
    # print c3.r
    # 
    # # TEST SUBTRACTION
    # print 'Should print a zero 3x3 matrix, and [0 0 0]'
    # c1 = Cw(np.array([1,2,3]))
    # c2 = c1 - c1
    # print c2.dr_wrt(c1).todense()
    # print c2.r
    # 
    # print 'Should print a positive 3x3 identity, a negative 3x3 identity, and [-3 -3 -3]'
    # c1 = Cw(np.array([1,2,3]))
    # c2 = Cw(np.array([4,5,6]))
    # c3 = c1 - c2
    # print c3.dr_wrt(c1).todense()
    # print c3.dr_wrt(c2).todense()
    # print c3.r
    
    # # TEST MULTIPLICATION
    # print 'Should print a [2 4 6] 3x3 and [1 4 9]'
    # c1 = Cw(np.array([1,2,3]))
    # c2 = c1 * c1
    # print c2.dr_wrt(c1).todense()
    # print c2.r
    # 
    # print 'Should print a [4 5 6] and [1 2 3] 3x3s and a [ 4 10 18]'
    # c1 = Cw(np.array([1,2,3]))
    # c2 = Cw(np.array([4,5,6]))
    # c3 = c1 * c2
    # print c3.dr_wrt(c1).todense()
    # print c3.dr_wrt(c2).todense()
    # print c3.r
    # 
    # TEST POWEROF
    # c1 = Cw(np.array([1,2,3]))
    # p1 = Cw(np.array([.5]))
    # c2 = c1 ** p1
    # print c2.dr_wrt(c1).todense()
    # print c2.dr_wrt(p1).todense()
    # print c2.r    
    # 
    # 
    # TEST DIVISION
    #print 'Should print a ???'
    # c1 = Cw(np.array([1,2,3]))
    # c2 = c1 / c1
    # print c2.dr_wrt(c1).todense()
    # print c2.r
    
    #print 'Should print a ???'
    # c1 = Cw(np.array([1.,2.,3.]))
    # c2 = Cw(np.array([4.,5.,6.]))
    # c3 = c1 / c2
    # print c3.dr_wrt(c1).todense()
    # print c3.dr_wrt(c2).todense()
    # print c3.r
    # # 

    # General nested case
    
    # for k in [0, 1]:
    #     a1 = Cw(np.array([1.,2.,3.]))
    #     a2 = Cw(np.array([4.,5.,6.]))
    #     a3 = Cw(np.array([7.,10.,12.]))
    # 
    #     #c3.show_tree()
    # 
    #     b1 = a1 + a1
    #     
    #     c1 = a1 + a1
    #     
    #     d1 = b1 + c1
    # 
    #     if k == 1:
    #         print 'beginning replacement'
    #         d1.remove_redundancies()
    #         print 'done with replacement'
    #         
    #     print d1.r
    #     print d1.dr_wrt(a1).todense()
    #     print d1.dr_wrt(a1).todense()
    #     print d1.dr_wrt(a1).todense()
    #     print d1.dr_wrt(a1).todense()
    #     print d1.dr_wrt(a1).todense()
    #     print d1.dr_wrt(a1).todense()
    #     print d1.dr_wrt(a1).todense()
    #     print d1.dr_wrt(a1).todense()
    #     print d1.dr_wrt(a1).todense()
    #     
    #     d1.show_tree()
    #     
    #     print '-----------------'


    a1 = Cw(np.array([1.,2.,3., 4.]))
    b1 = Cw(np.array([1.,2.,3., 5.]))
    c1 = a1 + a1
    d1 = b1 + b1
    
    d = c1 &  d1
    print(d.r)
    print(d.dr_wrt(b1).todense())


    

if __name__ == '__main__':
    main()
