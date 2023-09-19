# Stochastic-offset-enhanced-restricted-slice-excitation-and-refocusing-designs

This repository contains the official code for the MRM paper titled "Stochastic‐offset‐enhanced restricted slice excitation and 180° refocusing designs with spatially non‐linear ΔB0 shim array fields". Please note that this repository is currently under construction.

## About the Paper

The code in this repository is based on and extends the work presented in another MRM paper titled "Selective RF excitation designs enabled by time‐varying spatially non‐linear ΔB0 fields with applications in fetal MRI." The paper is divided into three key parts:

`### 1. Stochastic Offset Strategy`

We have adopted a stochastic offset strategy to address the 'fixed-point' artifacts in the resultant magnetization profile that can arise as a result of performing optimizations on fixed points as the spatial locations.

`### 2. Refocusing Designs`

Our work enables refocusing designs through the utilization of the decomposition property of the Bloch equation.

`### 3. Equivalent Crusher Gradient Formula`

To avoid extensive sub-voxel simulations for the effect of crusher gradients in refocusing designs, we have incorporated an equivalent crusher gradient formula into the loss function.

## Usage

The code is modified and derived from [AutoDiffPulses](https://github.com/tianrluo/AutoDiffPulses) which provides the optimization framework for RF pulse and B0 fields (linear gradient fields and shim array fields) with auto-differentiation.

We followed the same manner of [AutoDiffPulses](https://github.com/tianrluo/AutoDiffPulses). The actual optimization part is performed with Pytorch and wrapped with MATLAB.

### Key features for usage.

In `+mrphy/+beffective/rfgr2beff.m`, if you want to 



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
