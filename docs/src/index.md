# Introduction

Lirion is a Julia framework for modeling and solving Optimal Power Flow (OPF) problems under both **centralized** and **distributed** optimization paradigms, providing a practical environment for implementing established formulations and investigating novel decomposition and consensus strategies.

At present, Lirion unifies DC, AC, and SOCP OPF formulations under a modular architecture that supports both centralized models and distributed implementations based on ADMM and ALADIN, enabling systematic experimentation with network partitioning and coupling-variable choices.

!!! important "Mathematical Foundations and Derivations"
    Lirion is accompanied by a dedicated technical website created and maintained by Mateus Kaiss, the principal developer of the framework.

    The site consolidates mathematical formulations, derivations, and structural studies developed throughout his academic training and research activities, covering both fundamental optimization concepts and advanced theoretical developments that underpin Lirion’s design.

    Researchers, students, and practitioners interested in optimization and Optimal Power Flow, as well as anyone passionate about mathematical modeling, are encouraged to explore this resource.

    [Visit the Lirion Mathematical Foundations Website](https://your-link-here.com)


## Centralized Models

Lirion implements classical centralized OPF formulations:

- DC OPF
- AC OPF
- SOCP OPF

---

## Distributed Formulations

Beyond centralized models, Lirion provides distributed variants based on the **Alternating Direction Method of Multipliers (ADMM)**.

Distributed formulations are structured along three main dimensions:

- The base model (DC, AC, or SOCP);
- The choice of **coupling variables**;
- The **system partitioning strategy**.

---
