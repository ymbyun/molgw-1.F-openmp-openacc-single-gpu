[------------------------------------------------]: #
# MOLGW 1.F with OpenMP and OpenACC (Single GPU)
[------------------------------------------------]: #

[Multi-platform shared-memory parallelized MOLGW]: #

[## OpenMP (CPU) and OpenACC (GPU) implementations]: #

Hi there,

I parallelized and accelerated MOLGW 1 [https://molgw.org] using OpenMP (CPU) and OpenACC (GPU), respectively.
I already merged my OpenMP implemenation into MOLGW 2 [https://github.com/molgw/molgw/pull/5], and I'll merge my OpenACC implementation into MOLGW 3 soon.
For now, the OpenACC version of MOLGW runs on a single NVIDIA GPU.

I used my OpenMP-parallelized and OpenACC-accelerated MOLGW in the following papers:
- **Young-Moo Byun** and Jejoong Yoo, GPU acceleration of many-body perturbation theory methods in MOLGW with OpenACC, *Int. J. Quantum Chem.* **124**, e27345 (2024) [https://doi.org/10.1002/qua.27345] **NOTE:** This paper has been recognized as a top viewed article in IJQC in 2025.
- **Young-Moo Byun** and Serdar Ogut, Practical *GW* scheme for electronic structure of 3*d*-transition-metal monoxide anions: ScO<sup>-</sup>, TiO<sup>-</sup>, CuO<sup>-</sup>, and ZnO<sup>-</sup>, *J. Chem. Phys.* **151**, 134305 (2019) [https://doi.org/10.1063/1.5118671]

Currently, I'm working on a few things:
- A manual for how to compile and run OpenMP and OpenACC versions
- Multi-GPU support
- AMD GPU support
- OpenMP offloading support
- Replacing BLAS/LAPACK with NVBLAS or cuBLAS/cuSOLVER using OpenACC-CUDA interoperability
- Replacing Libint with LibintX using OpenACC-CUDA interoperability
- GPU acceleration of DOT_PRODUCT and MATMUL (There are multiple ways to do it)
  
[A paper for OpenMP and OpenACC implementations]: #

Please let me know if you have any questions.

Thanks,

Young-Moo Byun

Last updated on 2025/05/10
