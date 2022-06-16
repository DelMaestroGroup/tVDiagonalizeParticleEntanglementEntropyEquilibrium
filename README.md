# tVDiagonalizeParticleEntanglementEntropyEquilibrium
Equilibrium version of the ED code for the t-V model

### Example
 
    julia .\tVVp_Vdependence_q0R1PH1_IntF.jl --out "./out/" --tmp "./tmp/" --spatial --V-start -3.0 --V-end 3.0 --V-step 0.01 --t 1.0 --ee 1 24 12

### Required julia packages:
- ArgParse
-  Arpack  
-  MKL
-  ProgressBars 

### Usage
    julia tVVp_Vdependence_q0R1PH1_IntF.jl 
                        [--out FOLDER]
                        [--tmp FOLDER]
                        [--out-obdm FOLDER] 
                        [--g2] 
                        [--obdm]
                        [--spatial] 
                        [--skip-hoffdiag-saving]
                        [--skip-hoffdiag-loading] 
                        [--no-flush]
                        [--no-recompute-structure-matrix] 
                        [--pbc]
                        [--obc] 
                        [--V-start V_start] 
                        [--V-end V_end]
                        [--V-step V_step] 
                        [--V-list [V_LIST...]]
                        [--Vp Vp] 
                        [--t t] 
                        --ee ℓ M N

| positional arguments | Description |
| --- | ----------- |
| M | number of sites (type: Int64) |
| N | number of particles (type: Int64) |



| optional arguments | Description |
| ----- | ---------- |
| --out FOLDER | path to output folder |
| --tmp FOLDER | folder for hamiltonian storage location |
| --out-obdm FOLDER | folder for obdm storage location (if --obdm provided) |
| --out-obdm FOLDER | folder for obdm storage location (if --obdm provided) |
| --g2              |    output the pair correlation function ⟨ρ_iρ_0⟩ |
| --obdm            |    output the spatial dependence of the OBDM |
| --spatial         |    output the spatial entanglement entropy for ℓ= M/2 |
| --skip-hoffdiag-saving | do not save offdiagonal terms of Hamiltonian for V=0 (if already saved or should never be saved in combination with --skip-hoffdiag-loading) |
| --skip-hoffdiag-loading | do not load offdiagonal terms of Hamiltonian from V=0 |
| --no-recompute-structure-matrix | only compute structure matrix once and keep it in memory (can cause memory problems for large systems but can speed up calculation) |
| --no-flush       |     do not flush write buffer to output files in after computation for each V |
| -h, --help     |       show help message and exit| 


| boundary conditions | Description |
| ----- | ---------- |
| --pbc | periodic boundary conditions (default) |
| --obc | open boundary conditions |
 

| tV parameters | Description |
| ----- | ---------- |
| --V-start | start V (type: Float64, default: -2.0) |
| --V-end | end V (type: Float64, default: 2.0) |
| --V-step | step in V (type: Float64, default: 0.1) |
| --V-list | space separated list of interaction values, if set ignore other V params |
|--Vp| nnn interaction Vp (type: Float64, default: 0.0) |
|--t| t hopping value (type: Float64, default: 1.0)|
  
| entanglement entropy | Description |                  
| ----- | ---------- |
|--ee ℓ|compute all EEs with partition size ℓ (type: Int64)| 