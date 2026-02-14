<p align="center">
  <img src="docs/assets/lirion_logo.png" alt="Lirion.jl" 
  width="600"/>
</p>

<h3 align="center">
  Large-scale Integrated Research on Innovative Optimization of Networks

</h3>

A Julia package for solving Optimal Power Flow (OPF) problems using both centralized and distributed algorithms.  
Current focus: advanced ADMM (Alternating Direction Method of Multipliers) implementations, with plans to add more distributed methods.


## Overview

Lirion.jl provides a framework for power system optimization, supporting both traditional centralized OPF solutions and distributed optimization algorithms. The package is built on top of [JuMP.jl](https://github.com/jump-dev/JuMP.jl) and [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) for robust mathematical modeling and optimization.

### Key Features

- **Centralized OPF**: DC, AC, and SOCP (Second-Order Cone Programming) formulations
- **Distributed ADMM**: Multiple ADMM variants with extensive customization options
- **Advanced Partitioning**: Both spectral and METIS-based network partitioning
- **Multi-period Support**: Time-series optimization capabilities
- **Storage Integration**: Battery and prosumer modeling
- **Comprehensive Testing**: Extensive test suite with PGLib cases