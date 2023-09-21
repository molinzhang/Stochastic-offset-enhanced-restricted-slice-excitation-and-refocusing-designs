import torch
from torch import Tensor

import numpy as np

def pen_null(rf: Tensor) -> Tensor:
    """
    *INPUTS*
    - `rf`  (1, xy, nT, (nCoils))
    *OUTPUTS*
    - `pen` (1,)
    """
    return rf.new_zeros([])


def pen_l2(rf: Tensor) -> Tensor:
    pen = torch.norm(rf)**2
    return pen

def pen_l2_mat(gr: Tensor) -> Tensor:
    w = torch.zeros_like(gr,dtype = gr.dtype)
    w[:,:,0] = 1
    w[:,:,-1] = 1
    
    pen = torch.mean(torch.norm(gr * w,dim = 1) ** 2)
    return pen
    
def pen_inf_mat(gr: Tensor) -> Tensor:
    pen = torch.norm(gr, float('inf'))
    return pen

def crush_stop(Mx, My, w): # 1, nM, 3 # w is the weight of stop band.
    pen1 = Mx[:,:,0] - My[:,:,1]
    pen2 = Mx[:,:,1] + My[:,:,0]
    pen1 = pen1 * w
    pen2 = pen2 * w
    return pen1.norm()**2, pen2.norm()**2
    
    
def rf_smooth(rf: Tensor) -> Tensor: # rf shape [1,2,500]
    ### regularizer on fft of rf or rf gradient
    '''
    rf = rf.view(1, 500, 2)
    #element = torch.view_as_complex(rf)#torch.complex(rf[:,0,:], rf[:,1,:])
    fft_rf = torch.fft(rf, 1)
    weight = np.zeros(fft_rf.shape)
    weight[:,:250,:] = np.arange(0,1,250)
    weight[:,250:,:] = 1 - np.arange(0,1,250)
    weight = torch.tensor(weight, dtype = fft_rf.dtype).to(fft_rf.device)
    rf_weight = weight*fft_rf
    pen = torch.sqrt(rf_weight[:,:,0] ** 2 + rf[:,:,1] ** 2)
    '''
    
    #pen = torch.mean(torch.norm(fft_rf, p = 1,dim = 1) ** 2)
    ### regularizer on rf or rf gradient
    grad = rf[:,:,0:499] - rf[:,:,1:]
    #pen = torch.sqrt(rf[:,0,:] ** 2 + rf[:,1,:] ** 2)
    pen = torch.sqrt(grad[:,0,:] ** 2 + grad[:,1,:] ** 2)
    #pen = torch.norm(pen, float('inf'))
    pen = torch.mean(torch.norm(pen, dim = 1) ** 2)
    #pen = torch.norm(rf[:,0,:], p = 2)
    
    return pen