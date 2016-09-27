__all__ = ["sparse_type"]

import scipy.sparse as sp

def sparse_type(sp_mat):
    if sp.issparse(sp_mat):
        if sp.isspmatrix_csc(sp_mat):
            return 'csc'
        if sp.isspmatrix_csr(sp_mat):
            return 'csr'
        if sp.isspmatrix_bsr(sp_mat):
            return 'bsr'
        if sp.isspmatrix_lil(sp_mat):
            return 'lil'
        if sp.isspmatrix_dok(sp_mat):
            return 'dok'
        if sp.isspmatrix_coo(sp_mat):
            return 'coo'
        if sp.isspmatrix_dia(sp_mat):
            return 'dia'
    else:
        return 'non-sparse'
