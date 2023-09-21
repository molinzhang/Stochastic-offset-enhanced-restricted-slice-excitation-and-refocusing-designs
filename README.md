# Stochastic-offset-enhanced-restricted-slice-excitation-and-refocusing-designs

This repository contains the official code for the MRM paper titled ["Stochastic‐offset‐enhanced restricted slice excitation and 180° refocusing designs with spatially non‐linear ΔB0 shim array fields"](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.29827). 

The code in this repository is based on and extends the work presented in another MRM paper titled ["Selective RF excitation designs enabled by time‐varying spatially non‐linear ΔB0 fields with applications in fetal MRI."](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.29114) Please note that this repository is currently under construction.

## About the Paper

The main contributions of this paper are summarized into three key parts:

### `1. Stochastic Offset Strategy`

We have adopted a stochastic offset strategy to address the 'fixed-point' artifacts in the resultant magnetization profile that can arise as a result of performing optimizations on fixed points as the spatial locations.

### `2. Refocusing Designs under current framework`

Our work enables refocusing designs through the utilization of the decomposition property of the Bloch equation.

### `3. Equivalent Crusher Gradient Formula`

To avoid extensive sub-voxel simulations for the effect of crusher gradients in refocusing designs, we have incorporated an equivalent crusher gradient formula into the loss function.

## Usage

The code is modified and derived from [AutoDiffPulses](https://github.com/tianrluo/AutoDiffPulses) which provides the optimization framework for RF pulse and B0 fields (linear gradient fields and shim array fields) with auto-differentiation.

We followed the same manner of [AutoDiffPulses](https://github.com/tianrluo/AutoDiffPulses). The actual optimization part is performed with Pytorch and wrapped with MATLAB.
  
Our work enables both excitation and refocusing designs by optimizing Rf pulse and time-varying $\Delta B_0$ shim array fields. Note that we used additional linear gradient fields but fixed them during the optimization. Optimizing linear gradient fields yields worse results. 

For excitation, we propose a two-stage optimization strategy. please first use codes under folder `excitation` to optimize RF pulse and shim current get desired magnitude of excited magnetizations without considering phase distributions. Next, please use the codes under folder `rewinding` to rewind the phase of excited magnetizations. Load optimized RF pulse and shim current from the first stage in the right place. Note that only shim current are optimized and RF pulse is set to zero for the rewinding stage. For more details, please check our paper.

For refocusing, please first use codes under folder `refocusing` to optimize RF pulse and shim current get desired refocusing magentization.

### Key features for usage.

`excitation/adpulses_sgd` and `excitation/mrphy_sgd` are python packages. A simple way to install it is to move them under the anaconda enviroment, `~/anaconda3/envs/<your env name>/lib/<your python version>/site-packages/`. Same for `rewinding` and  `refocusing`.

You could change `.requires_grad` for RF pulse and shim current to optimize the variables you want in `optimizer.py` under `adpulses`.

If you want to used shim array fields, you have to load the array field maps (Hz/A) in `mrphy/beffective.py` within function `rfgr2beff`. 



## Citation

If you find this work helpful, please consider citing the following two papers:

1. This Paper: "Stochastic‐offset‐enhanced restricted slice excitation and 180° refocusing designs with spatially non‐linear ΔB0 shim array fields."
```bibtex
@article{zhang2023stochastic,
  title={Stochastic-offset-enhanced restricted slice excitation and 180° refocusing designs with spatially non-linear $\Delta$B0 shim array fields},
  author={Zhang, Molin and Arango, Nicolas and Arefeen, Yamin and Guryev, Georgy and Stockmann, Jason P and White, Jacob and Adalsteinsson, Elfar},
  journal={Magnetic Resonance in Medicine},
  year={2023},
  publisher={Wiley Online Library}
}
```

2. Original Paper: "Selective RF excitation designs enabled by time‐varying spatially non‐linear ΔB0 fields with applications in fetal MRI."
```bibtex
@article{zhang2022selective,
  title={Selective RF excitation designs enabled by time-varying spatially non-linear $\Delta$ B 0 fields with applications in fetal MRI},
  author={Zhang, Molin and Arango, Nicolas and Stockmann, Jason P and White, Jacob and Adalsteinsson, Elfar},
  journal={Magnetic Resonance in Medicine},
  volume={87},
  number={5},
  pages={2161--2177},
  year={2022},
  publisher={Wiley Online Library}
}
```



## License

This project is licensed under the [MIT License](LICENSE).
