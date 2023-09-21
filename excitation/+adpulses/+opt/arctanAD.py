# arctan.py

import torch
import os
from adpulses_sgd import io, optimizers, metrics, penalties
os.environ["CUDA_VISIBLE_DEVICES"] = "1"
if __name__ == "__main__":
    import sys

    if len(sys.argv) <= 1:  # mode DEBUG
        import os
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

    m2pName = ('m2p.mat' if len(sys.argv) <= 1 else sys.argv[1])
    p2mName = ('p2m.mat' if len(sys.argv) <= 2 else sys.argv[2])
    gpuID = ('0' if len(sys.argv) <= 3 else sys.argv[3])
    
    print('gpuid!!!\t',gpuID,'gpuid!!!\t')
    
    # %% load
    if gpuID == '-1':
        dkw = {'device': torch.device('cpu'), 'dtype': torch.float32}
    else:
        dkw = {'device': torch.device('cuda:'+gpuID), 'dtype': torch.float32}
    
    
    target, cube, cube1, cube2, pulse, arg = io.m2p(m2pName, **dkw)

    def dflt_arg(k, v, fn):
        return (fn(k) if ((k in arg.keys()) and (arg[k].size > 0)) else v)

    arg['doRelax'] = dflt_arg('doRelax', True, lambda k: bool(arg[k].item()))

    arg['b1Map_'] = dflt_arg('b1Map_', None,
                             lambda k: f_tensor(f_c2r_np(arg[k], -2)))

    arg['niter'] = dflt_arg('niter', 1, lambda k: arg[k].item())
    arg['niter_gr'] = dflt_arg('niter_gr', 0, lambda k: arg[k].item())
    arg['niter_rf'] = dflt_arg('niter_rf', 0, lambda k: arg[k].item())

    eta = dflt_arg('eta', 4, lambda k: float(arg[k].item()))
    print('eta: ', eta)

    err_meth = dflt_arg('err_meth', 'l2xy', lambda k: arg[k].item())
    pen_meth = dflt_arg('pen_meth', 'l2', lambda k: arg[k].item())

    err_hash = {'null': metrics.err_null,
                'l2xy': metrics.err_l2xy, 'ml2xy': metrics.err_ml2xy,
                'l2z': metrics.err_l2z}
    pen_hash = {'null': penalties.pen_null, 'l2': penalties.pen_l2}

    fn_err, fn_pen = err_hash[err_meth], pen_hash[pen_meth]

    # %% pulse design
    kw = {k: arg[k] for k in ('b1Map_', 'niter', 'niter_gr', 'niter_rf',
                              'doRelax')}

    pulse, optInfos = optimizers.arctanLBFGS(target, cube, cube1, cube2, pulse,
                                             fn_err, fn_pen, eta=eta, **kw)

    # %% saving
    #print('python saved pulse',pulse.gr.shape)
    io.p2m(p2mName, pulse, optInfos)
