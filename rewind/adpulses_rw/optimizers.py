from typing import Tuple, Callable, Optional
from time import time
from numbers import Number

import numpy as np
from torch import optim, Tensor
import mrphy_rw
from mrphy_rw.mobjs import SpinCube, Pulse
import torch

import numpy as np
import random

from adpulses_rw.penalties import pen_l2_mat, pen_inf_mat, rf_smooth, crush_stop

from scipy.io import savemat
import os

def arctanLBFGS(
        target: dict, cube: SpinCube, cube1: SpinCube, cube2: SpinCube, pulse: Pulse,
        fn_err: Callable[[Tensor, Tensor, Optional[Tensor]], Tensor],
        fn_pen: Callable[[Tensor], Tensor],
        niter: int = 1, niter_gr: int = 1, niter_rf: int = 1,
        eta: Number = 4., b1Map_: Optional[Tensor] = None, doRelax: bool = True
        ) -> Tuple[Pulse, dict]:
    r"""Joint RF/GR optimization via direct arctan trick

    Usage:
        ``arctanLBFGS(target, cube, pulse, fn_err, fn_pen; eta=eta)``

    Inputs:
        - ``target``: dict, with fields:
            ``d_``: `(1, nM, xy)`, desired excitation;
            ``weight_``: `(1, nM)`.
        - ``cube``: mrphy_sgd.mobjs.SpinCube.
        - ``pulse``: mrphy_sgd.mobjs.Pulse.
        - ``fn_err``: error metric function. See :mod:`~adpulses.metrics`.
        - ``fn_pen``: penalty function. See :mod:`~adpulses.penalties`.
    Optionals:
        - ``niter``: int, number of iterations.
        - ``niter_gr``: int, number of LBFGS iters for updating *gradients*.
        - ``niter_rf``: int, number of LBFGS iters for updating *RF*.
        - ``eta``: `(1,)`, penalization term weighting coefficient.
        - ``b1Map_``: `(1, nM, xy,(nCoils))`, a.u., transmit sensitivity.
        - ``doRelax``: [T/f], whether accounting relaxation effects in simu.
    Outputs:
        - ``pulse``: mrphy_sgd.mojbs.Pulse, optimized pulse.
        - ``optInfos``: dict, optimization informations.
    """
    # Set up: Interior mapping
    print(pulse)
    tsl_rw =  mrphy_rw.utils.s2ts(mrphy_rw.utils.g2s(pulse.rw, pulse.dt), pulse.smax)    

    opt_sl = optim.LBFGS([tsl_rw], lr=0.005, max_iter=20, history_size=100,
                         tolerance_change=1e-6,
                         line_search_fn='strong_wolfe') # lr = 0.1, 20? ## 4mm lr = 0.5 / 0.1 # lr = 0.005
                         
    scheduler_gr = optim.lr_scheduler.MultiStepLR(opt_sl, milestones=[240, 280], gamma=0.5)
      
                         
    tsl_rw.requires_grad = True


    # Set up: optimizer
    loss_hist = np.full((niter*(niter_gr+niter_rf),), np.nan)
    Md_, w_ = target['d_'], target['weight_'].sqrt()  # (1, nM, 3), (1, nM)
    Md1_, w1_ = target['d1_'], target['weight1_'].sqrt()  # (1, nM, 3), (1, nM)
    Md2_, w2_ = target['d2_'], target['weight2_'].sqrt()  # (1, nM, 3), (1, nM)
    
    
    print(Md_.shape, Md1_.shape, Md2_.shape)
    print(w_.shape, w1_.shape, w2_.shape)
    
    rfmax, smax = pulse.rfmax, pulse.smax

    def fn_loss(cube, pulse, idx, offset_loc_ = None):
    
        Mr_ = cube.applypulse(pulse, b1Map_=b1Map_, doRelax=doRelax, offset_loc_ = offset_loc_)
        if idx == 1:
            loss_err, loss_pen, loss_pengr, loss_edge, loss_rfsm = fn_err(Mr_, Md_, w_=w_), fn_pen(pulse.rf), pen_inf_mat(pulse.gr), pen_l2_mat(pulse.gr), pen_l2_mat(pulse.gr)#pen_l2(pulse.rf)
        if idx == 2:
            loss_err, loss_pen, loss_pengr, loss_edge, loss_rfsm = fn_err(Mr_, Md1_, w_=w_), fn_pen(pulse.rf), pen_inf_mat(pulse.gr), pen_l2_mat(pulse.gr), pen_l2_mat(pulse.gr)#pen_l2(pulse.rf)
        if idx == 3:
            loss_err, loss_pen, loss_pengr, loss_edge, loss_rfsm = fn_err(Mr_, Md2_, w_=w2_), fn_pen(pulse.rf), pen_inf_mat(pulse.rw), pen_l2_mat(pulse.rw), pen_l2_mat(pulse.rw)
            #loss_err, loss_pen, loss_pengr, loss_edge, loss_rfsm = fn_err(Mr_, Md2_, w_=w2_), fn_pen(pulse.rf), pen_inf_mat(pulse.gr), pen_l2_mat(pulse.gr), pen_l2_mat(pulse.gr)#pen_l2(pulse.rf)
        return loss_err, loss_pen, loss_pengr, loss_edge, loss_rfsm, Mr_
    
    def crush_mod_loss(Mx, My, w):
        return crush_stop(Mx, My, w)
    

    log_col = '\n#iter\t ‖ elapsed time\t ‖ error\t ‖ penalty\t ‖ total loss'

    def logger(i, t0, loss_err1, loss_err2, loss_err3, loss_pen, loss_pengr, loss_edge, loss_rfsm):
        loss = 1 * loss_err1 + 1*loss_err2 +  0.2* loss_err3 + eta*loss_pen + 0.0*loss_pengr + 10.0*loss_edge + 0.0 * loss_rfsm
        print("%i\t | %.1f  \t | %.3f\t | %.3f\t | %.3f\t | %.3f\t | %.3f\t | %.3f\t| %.3f| %.3f" %(i, time()-t0, loss_err1.item(), loss_err2.item(), loss_err3.item(),loss_pen.item(), loss_pengr.item(), loss_edge.item(), loss_rfsm.item(), loss.item()))
        return loss

    # Optimization
    t0 = time()
    for i in range(niter):

        if not (i % 5):
            print(log_col)

        log_ind = 0
        
        off_coef = 1
        offset_loc_zero = torch.ones((1, 3,3,480, 3)) * off_coef
        #pos_choice = random.choice(range(0,500))
        #print(pos_choice)
        offset_loc_f1 = torch.rand((1, 3, 3, 480, 3))#torch.rand((1, 88, 88, 60, 3))#offset_loc_f1 = torch.rand((1, 88, 88, 60, 3))   #torch.zeros((1, 30, 30, 20, 3)) + pos_choice * 0.002#torch.rand((1, 30, 30, 20, 3))        
        #offset_loc_f1 = torch.linspace(0,1,steps = 861)
        def closure():
            opt_sl.zero_grad()

            pulse.rw = mrphy_rw.utils.s2g(mrphy_rw.utils.ts2s(tsl_rw, smax), pulse.dt)
            
            pulse.gr = torch.cat((pulse.lr,torch.cat((pulse.nr, pulse.rw), -1)), 1)
            

            ''' # refocus
            loss_err1, loss_pen, loss_pengr, loss_edge, loss_rfsm, M1 = fn_loss(cube, pulse, 1, offset_loc_f1)
            loss_err2, loss_pen, loss_pengr, loss_edge, loss_rfsm, M2 = fn_loss(cube1, pulse, 2, offset_loc_f1)
            loss_crushx, loss_crushy = crush_mod_loss(M1, M2, w1_)
            
            #loss = 1 * loss_err1 + 1 *loss_err2 + 0.1* loss_err3 + eta*loss_pen + 0.0*loss_pengr + 5.0*loss_edge + 0.0 * loss_rfsm
            loss = 1 * loss_err1 + 1 *loss_err2 + 1* (loss_crushx + loss_crushy) + eta*loss_pen + 0.0*loss_pengr + 100.0*loss_edge + 0.0 * loss_rfsm
            '''
            
            loss_err1, loss_pen, loss_pengr, loss_edge, loss_rfsm, M3 = fn_loss(cube2, pulse, 3, offset_loc_f1)
            loss = 1 * loss_err1 + 1 *loss_pengr + 100.0*loss_edge

            loss.backward()
            return loss
        

        
        print('gr-loop: ', niter_gr)
        for _ in range(niter_gr):
            #opt_sl.step(closure)
            '''
            off_coef = (1/niter_gr) * _
            offset_loc_zero = torch.ones((1, 3,3,480, 3)) * off_coef
            


            loss_err1, loss_pen, loss_pengr, loss_edge, loss_rfsm, M1 = fn_loss(cube2, pulse, 3,offset_loc_zero)
            loss = logger(i, t0, loss_err1, loss_err1, loss_err1, loss_pen, loss_pengr, loss_edge, loss_rfsm)
              
            MM1 = M1.cpu().detach().numpy()
            ww_ = w2_.cpu().detach().numpy()
            savemat(os.path.join('/home/molin/shimm_nick/op_rewind/results/minor_revision/pulserw13car/', str(off_coef)+'ex.mat'), {'Mx': MM1, 'mask_in':ww_})
              

            '''
            for j in range(8):
                for k in range(8):
                    #opt_sl.step(closure)
            
                    off_coef = (1/niter_gr) * _
                    #offset_loc_zero = torch.ones((1, 88,88,60, 3)) * off_coef
                    offset_loc_zero = torch.ones((1, 3,3,480, 3))
                    offset_loc_zero[:,:,:,:,0] = offset_loc_zero[:,:,:,:,0] * off_coef
                    offset_loc_zero[:,:,:,:,1] = offset_loc_zero[:,:,:,:,1] * j/8
                    offset_loc_zero[:,:,:,:,2] = offset_loc_zero[:,:,:,:,2] * k/8
            
      
                    loss_err1, loss_pen, loss_pengr, loss_edge, loss_rfsm, M1 = fn_loss(cube2, pulse, 3,offset_loc_zero)
                    loss = logger(i, t0, loss_err1, loss_err1, loss_err1, loss_pen, loss_pengr, loss_edge, loss_rfsm)
            
                    MM1 = M1.cpu().detach().numpy()
                    ww_ = w2_.cpu().detach().numpy()
                    #savemat(os.path.join('/home/molin/shimm_nick/decompose_refocus/demo/ISMRM/compare/refocus_test/selective/excitation/trash/pulse11/', str(off_coef)+'ex.mat'), {'Mx': MM1, 'mask_in':ww_})
                    savemat(os.path.join('/home/molin/shimm_nick/op_rewind/results/minor_revision/pulserw13car/', "{:.4f}".format(off_coef) + "{:.4f}".format(j) + "{:.4f}".format(k)+'ex.mat'), {'Mx': MM1, 'mask_in':ww_})
            
                    #loss_hist[i*(niter_gr+niter_rf)+log_ind] = loss.item()
                    #log_ind += 1
            
        scheduler_gr.step()
        
        
    print('\n== Results: ==')
    print(log_col)
    #loss = logger(i, t0, loss_err1,loss_err2,loss_err3, loss_pen, loss_pengr, loss_edge, loss_rfsm)

    optInfos = {'loss_hist': loss_hist}
    print(torch.max(pulse.gr[:,4:]))
    return pulse, optInfos
