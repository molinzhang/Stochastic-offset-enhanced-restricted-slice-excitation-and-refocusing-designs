# Stochastic-offset-enhanced-restricted-slice-excitation-and-refocusing-designs

This repository contains the official code for the MRM paper titled "Stochastic‐offset‐enhanced restricted slice excitation and 180° refocusing designs with spatially non‐linear ΔB0 shim array fields". Please note that this repository is currently under construction.

## About the Paper

The code in this repository is based on and extends the work presented in another MRM paper titled "Selective RF excitation designs enabled by time‐varying spatially non‐linear ΔB0 fields with applications in fetal MRI." The paper is divided into three key parts:

### 1. Stochastic Offset Strategy

We have adopted a stochastic offset strategy to address the 'fixed-point' artifacts that can arise during the optimization of fixed points in the magnetization profile.

### 2. Refocusing Designs

Our work enables refocusing designs through the utilization of the decomposition property of the Bloch equation.

### 3. Equivalent Crusher Gradient Formula

To avoid extensive sub-voxel simulations for the effect of crusher gradients in refocusing designs, we have incorporated an equivalent crusher gradient formula into the loss function.

## Usage

Instructions on how to use the code and replicate the experiments from the paper will be provided here once the repository is complete.

## Citation

If you find this paper helpful, please consider citing the following two papers:

1. Original Paper: "Selective RF excitation designs enabled by time‐varying spatially non‐linear ΔB0 fields with applications in fetal MRI."
```
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

2. Extended Paper: "Stochastic‐offset‐enhanced restricted slice excitation and 180° refocusing designs with spatially non‐linear ΔB0 shim array fields."
```
@article{zhang2023stochastic,
  title={Stochastic-offset-enhanced restricted slice excitation and 180° refocusing designs with spatially non-linear $\Delta$B0 shim array fields},
  author={Zhang, Molin and Arango, Nicolas and Arefeen, Yamin and Guryev, Georgy and Stockmann, Jason P and White, Jacob and Adalsteinsson, Elfar},
  journal={Magnetic Resonance in Medicine},
  year={2023},
  publisher={Wiley Online Library}
}
```

## License

This project is licensed under the [MIT License](LICENSE).
