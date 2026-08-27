
# Polarization-aware Reconfigurable Antenna Aided Wireless Communications

This repository contains the MATLAB code for the paper:  
**"Polarization-aware Reconfigurable Antenna Aided Wireless Communications"**, which can be found on the link: https://ieeexplore.ieee.org/document/11667052

![image](sys.png)

The code implements **Alternating Optimization (AO)** algorithms combined with **semidefinite programming (SDP)** and **Riemannian Manifold Optimization (RCG)** for Reconfigurable Antenna systems. It compares the following schemes:

* **Proposed RA Scheme (Proposed Rot+Pol RA + DBF):** This is the proposed architecture where the BS antennas support 3D rotation, and both the BS and user antennas support flexible polarization reconfiguration. Consequently, the digital beamforming matrix $\mathbf{W}$, the rotation matrices $\{\mathbf{R}_m\}$, and the polarization states $(\mathbf{V}, \mathbf{U})$ are jointly optimized to minimize the total transmit power.

* **Rotation-Only RA Scheme (Rotation-Only RA + DBF):** In this scheme, the BS antennas are rotatable but lack polarization reconfigurability. The digital beamforming $\mathbf{W}$ and rotation matrices $\{\mathbf{R}_m\}$ are optimized, while the transmit and receive polarization states are fixed. This baseline isolates the power saving gain provided by rotation.

* **Boresight-Only RA Scheme (Boresight-Only RA + DBF):** The BS antennas optimize only their boresight directions $\mathbf{r}_{m,3}$. The orientation of the polarization directions, $\mathbf{r}_{m,1}$ and $\mathbf{r}_{m,2}$, is deterministically coupled to $\mathbf{r}_{m,3}$ assuming zero rotation around the boresight axis. This benchmark represents conventional rotatable antennas that prioritize only directional gain.

* **Conventional Fixed-Antenna Scheme (Fixed UPA + DBF):** This represents the conventional fixed UPA system. Both the rotation and polarization states are fixed. Only the digital beamforming matrix $\mathbf{W}$ is optimized. This scheme serves as the performance lower bound.


## 1. Reproduce Main Results

To generate the simulation figure for **Transmit Power versus Required Rate**, simply run the following script:
```matlab
Pvs_require.m
```
Note: For other simulations mentioned in the paper, please modify the system parameters accordingly.
## Bibtex
```python
@ARTICLE{zhang2026polarization,
  author={Zhang, Chao and Zhou, Hu and Long, Ruizhe and Liang, Ying-Chang and Soong, Boon-Hee},
  journal={IEEE Transactions on Wireless Communications}, 
  title={Polarization-aware Reconfigurable Antenna Aided Wireless Communications}, 
  year={2026},
  volume={},
  number={},
  pages={1-1},
  keywords={Polarization;Antennas;Gain;Optimization;Modeling;Matrices;Arrays;Design methodology;Array signal processing;Vectors;Reconfigurable antenna (RA);Polarization and 3D rotation optimization;directional gain pattern},
  doi={10.1109/TWC.2026.3724761}}
```

If you have any questions, please feel free to reach me at: zhang_chao@std.uestc.edu.cn
