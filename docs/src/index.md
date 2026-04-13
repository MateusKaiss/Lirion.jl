```@raw html
<img src="assets/logo.png" alt="Lirion.jl logo"
     style="max-width: 520px; width: 100%; margin: 0 auto 1.25rem; display: block;" />

```


```@meta
# This prevents certain rendering issues in Documenter
```

# Introduction

Lirion is a Julia framework for modeling and solving Optimal Power Flow (OPF) problems under both **centralized** and **distributed** optimization paradigms, providing a practical environment for implementing established formulations and investigating novel decomposition and consensus strategies.

At present, Lirion unifies DC, AC, and SOCWR OPF formulations within a modular architecture that spans centralized models and distributed implementations based on ADMM and ALADIN, enabling systematic experimentation with partitioning strategies and coupling-variable design.

!!! important "Mathematical Foundations and Derivations"
    Lirion is accompanied by a dedicated technical website maintained by Mateus Kaiss, the principal developer of the framework.

    The site consolidates mathematical formulations, derivations, and structural studies developed throughout his academic training and research activities, covering both fundamental optimization concepts and advanced theoretical developments that underpin Lirion’s design.

    Researchers, students, and practitioners interested in optimization and Optimal Power Flow — as well as anyone passionate about mathematical modeling — are encouraged to explore this resource.

    [Visit the Lirion Mathematical Foundations Website](https://your-link-here.com)

## How this documentation is structured

This documentation is organized into three complementary layers:

- **Overview and Guides**: high-level explanations of Lirion’s modeling workflow and recommended usage patterns.
- **Modeling and Algorithms**: centralized formulations, distributed variants (ADMM/ALADIN), and how decomposition choices are expressed in Lirion.
- **API Reference**: the complete list of functions and types for advanced usage and extension.

## Centralized Models

Lirion implements classical centralized OPF formulations:

- DC OPF
- AC OPF
- SOCWR OPF

---

## Distributed Formulations

Beyond centralized models, Lirion provides distributed variants based on the **Alternating Direction Method of Multipliers (ADMM)**.

Distributed formulations are structured along three main dimensions:

- The base model (DC, AC, or SOCWR);
- The choice of **coupling variables**;
- The **system partitioning strategy**.

---
