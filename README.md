<p align="center">
  <img src="docs/assets/logo_lirion.png" alt="Lirion.jl" width="420"/>
</p>

<h3 align="center">
  Large-scale Integrated Research in Intelligent Optimization of Networks
</h3>


# Lirion.jl

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

## Acknowledgments

Lirion.jl initially emerged as part of the undergraduate thesis (TCC) of its principal developer in the field of power systems. The project later evolved into a research project funded by the **CNPq (Brazilian National Council for Scientific and Technological Development)**.  

The package was conceived, designed, and developed by **Mateus Kaiss**, as its principal author and maintainer, with support from:

- **Kristie Kaminsky Kuster** — academic advisor during the thesis, providing theoretical and practical guidance, as well as emotional support throughout the development process. I will be forever grateful for her unwavering support and mentorship. 
- **Daniel Gebbran** — contributions to distributed ADMM implementations and testing.

We also acknowledge the inspiration and foundational work of the following open-source projects:  
- [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) — solid foundation for power system modeling  
- [PowerModelsADA.jl](https://github.com/mkhraijah/PowerModelsADA.jl) — reference for distributed algorithms  
- [JuMP.jl](https://github.com/jump-dev/JuMP.jl) — mathematical optimization modeling framework  
- [PGLib](https://github.com/power-grid-lib/pglib-opf) — benchmark test cases


## Comparison with Related Packages

### vs. [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl)

**Lirion.jl** builds upon PowerModels.jl's solid foundation but offers several key advantages:

- **Distributed Algorithms**: While PowerModels focuses on centralized optimization, Lirion provides extensive distributed implementations, currently focusing on ADMM

**Similarities**:
- Both use JuMP.jl for mathematical modeling
- Both include variations for relaxed variations (DC, SOCP)
- Compatible with the same optimizers (Ipopt, Gurobi, OSQP, SCS)
- Support for standard MATPOWER case formats
- Similar API design patterns

### vs. [PowerModelsADA.jl](https://github.com/mkhraijah/PowerModelsADA.jl)

**Lirion.jl** offers more depth in ADMM implementations but less breadth in distributed methods:

**Lirion.jl Advantages**:
- **Deeper ADMM Implementation**: More sophisticated ADMM variants (hybrid, branch consensus)
- **Advanced Partitioning**: Both spectral and METIS-based partitioning with power system-specific optimizations
- **Centralized OPF**: Also implements centralized solutions (DC, AC, SOCP) and allows for easy comparisons across methods
- **Active Development**: Currently under active development with recent improvements
- **Extensive Customization**: Fine-grained control over ADMM parameters and convergence criteria

**Current PowerModelsADA.jl Advantages**:
- **Broader Algorithm Support**: Implements multiple distributed algorithms (ADMM, ATC, APP, ALADIN)

**Similarities**:
- Both focus on distributed optimization for power systems
- Built on PowerModels.jl foundation
- Support for standard test cases

## Installation

### Prerequisites

- Julia 1.6 or later
- Required system dependencies for optimizers (e.g., Ipopt, Gurobi)

### Basic Installation

```julia
using Pkg
Pkg.add("Lirion")
```

### Development Installation

```julia
using Pkg
Pkg.develop("Lirion")
```

### Dependencies

The package requires several Julia packages including:
- JuMP.jl (mathematical modeling)
- PowerModels.jl (power system data structures)
- Ipopt.jl, Gurobi.jl, OSQP.jl, SCS.jl (optimizers)
- Metis.jl (graph partitioning)
- PlotlyLight.jl (visualization)

## Quick Start Guide

### Setting up the Environment

1. **Clone the repository**:
```bash
git clone https://github.com/your-username/Lirion.jl.git
cd Lirion.jl
```

2. **Activate the Julia environment**:
```julia
julia --project=.
```

3. **Instantiate dependencies**:
```julia
using Pkg
Pkg.instantiate()
```

### Basic Usage

#### Centralized OPF

```julia
using Lirion

# Solve a DC OPF problem
args = [
    "path/to/case.m",
    "--model", "dc",
    "--optimizer", "Ipopt.Optimizer"
]
result = Lirion.main(args)
```

#### Distributed ADMM

```julia
# Solve using hybrid ADMM with DC base model
args = [
    "path/to/case.m",
    "--model", "distributed-admm-hybrid-dc",
    "--optimizer", "Ipopt.Optimizer",
    "--admm-rho", "1000.0",
    "--admm-max-iter", "350"
]
result = Lirion.main(args)
```

### Available Models

#### Centralized Models
- `dc` - DC Optimal Power Flow
- `ac` - AC Optimal Power Flow  
- `socp` - Second-Order Cone Programming

#### Distributed ADMM Models
- `distributed-admm-dc` - Standard ADMM with DC
- `distributed-admm-ac` - Standard ADMM with AC
- `distributed-admm-socp` - Standard ADMM with SOCP
- `distributed-admm-branch-dc` - Branch-based ADMM with DC
- `distributed-admm-branch-ac` - Branch-based ADMM with AC
- `distributed-admm-branch-socp` - Branch-based ADMM with SOCP
- `distributed-admm-hybrid-dc` - Hybrid ADMM with DC
- `distributed-admm-hybrid-ac` - Hybrid ADMM with AC
- `distributed-admm-hybrid-socp` - Hybrid ADMM with SOCP

### Command Line Interface

The package provides a comprehensive command-line interface:

```bash
julia main.jl case.m --model distributed-admm-hybrid-dc --optimizer Ipopt.Optimizer
```

#### Available Arguments

- `--model`: Model type (see above)
- `--optimizer`: Optimizer (Ipopt.Optimizer, Gurobi.Optimizer, OSQP.Optimizer, SCS.Optimizer)
- `--plot-outputs`: Generate plots
- `--outputs`: Generate CSV output files
- `--admm-rho`: Default penalty parameter (default: 1000.0)
- `--admm-max-iter`: Maximum ADMM iterations (default: 350)
- `--admm-eps-primal`: Primal tolerance (default: 1e-7)
- `--admm-eps-dual`: Dual tolerance (default: 1e-7)
- `--admm-alpha`: Over-relaxation parameter (1 < α < 2)

## Network Partitioning

### Spectral Partitioning

The package includes a sophisticated spectral partitioning algorithm optimized for power systems:

```julia
using Lirion
include("scripts/partition_cases.jl")

# Partition cases with spectral method
partition_cases(
    "path/to/cases/";
    ks = [2, 3, 4, 5],
    seed = 42,
    balance_tol = 0.10,
    min_zone_ratio = 0.15
)
```

### METIS Partitioning

For large-scale systems, METIS-based partitioning is available:

```julia
using Lirion
include("scripts/partition_cases.jl")

# Partition cases with METIS
partition_cases_metis(
    "path/to/cases/";
    ks = :auto,  # or specific numbers like [2, 3, 4, 5]
    max_zone_size = 1500,
    target_edgecut = 0.02,
    balance_tol = 0.10
)
```

#### METIS Setup Requirements

1. **Install METIS system library**:
   - Ubuntu/Debian: `sudo apt-get install libmetis-dev`
   - macOS: `brew install metis`
   - Windows: Download from [METIS website](http://glaros.dtc.umn.edu/gkhome/metis/metis/download)

2. **Install Julia METIS package**:
```julia
using Pkg
Pkg.add("Metis")
```

3. **Usage**:
```julia
# Automatic partitioning based on system size
partition_cases_metis(
    "test/matpower_cases/";
    ks = :auto,
    max_zone_size = 5000,
    target_edgecut = 0.02
)

# Manual partitioning for specific zone counts
partition_cases_metis(
    "test/matpower_cases/";
    ks = [2, 3, 4, 5],
    max_zone_size = 2000,
    target_edgecut = 0.03
)
```

## Examples

### Example 1: Basic DC OPF

```julia
using Lirion

# Solve IEEE 30-bus case with DC OPF
args = [
    "test/matpower_cases/pglib_opf_case30_ieee.m",
    "--model", "dc",
    "--optimizer", "Ipopt.Optimizer"
]
result = Lirion.main(args)
```

### Example 2: Distributed ADMM with Custom Parameters

```julia
using Lirion

# Solve with hybrid ADMM and custom parameters
args = [
    "test/matpower_cases/pglib_opf_case118_2zones_ieee.m",
    "--model", "distributed-admm-hybrid-dc",
    "--optimizer", "Ipopt.Optimizer",
    "--admm-rho", "2000.0",
    "--admm-max-iter", "500",
    "--admm-eps-primal", "1e-6",
    "--admm-eps-dual", "1e-6",
    "--admm-alpha", "1.5"
]
result = Lirion.main(args)
```

### Example 3: AC OPF with Plotting

```julia
using Lirion

# Solve AC OPF and generate plots
args = [
    "test/matpower_cases/pglib_opf_case57_ieee.m",
    "--model", "ac",
    "--optimizer", "Ipopt.Optimizer",
    "--plot-outputs"
]
result = Lirion.main(args)
```

## Testing

Run the test suite to verify installation:

```julia
using Pkg
Pkg.test("Lirion")
```

Or run specific tests:

```julia
julia test/runtests.jl
```

## Performance Considerations

### Threading Configuration

The package automatically configures threading for optimal performance:

```julia
# Check current configuration
println("Julia threads: ", Threads.nthreads())
println("BLAS threads: ", BLAS.get_num_threads())
```

### Memory Management

For large-scale problems:
- Use appropriate partitioning to control subproblem sizes
- Consider using `--admm-max-iter` to limit iterations
- Monitor memory usage with large test cases

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for:

- Bug fixes
- New features
- Performance improvements
- Documentation enhancements
- Additional test cases

### Development Setup

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Submit a pull request

## Citation

If you use Lirion.jl in your research, please cite:

```bibtex
@software{Lirion2024,
  title={Lirion.jl: A Julia Package for Distributed Optimal Power Flow},
  author={Kaiss, Mateus and Gebbran, Daniel},
  year={2024},
  url={https://github.com/your-username/Lirion.jl}
}
```

## Related Packages

- [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) - Core power system modeling
- [PowerModelsADA.jl](https://github.com/mkhraijah/PowerModelsADA.jl) - Alternative distributed algorithms
- [JuMP.jl](https://github.com/jump-dev/JuMP.jl) - Mathematical optimization modeling
- [InfrastructureModels.jl](https://github.com/lanl-ansi/InfrastructureModels.jl) - Infrastructure modeling utilities

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


## Roadmap

- [ ] Additional distributed algorithms (ATC, APP, ALADIN)
- [ ] GPU acceleration support
- [ ] Real-time optimization capabilities
- [ ] Enhanced visualization tools
- [ ] More comprehensive documentation
- [ ] Performance benchmarking suite

---

For more information, please visit our [documentation](docs/) or open an [issue](https://github.com/your-username/Lirion.jl/issues) for questions and support.
