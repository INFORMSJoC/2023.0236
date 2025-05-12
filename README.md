[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Tightening Quadratic Convex Relaxations for the AC Optimal Transmission Switching Problem (2023.0236)

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper [Tightening Quadratic Convex Relaxations for the AC Optimal Transmission Switching Problem](https://doi.org/10.1287/ijoc.2023.0236) by C. Guo, H. Nagarajan and M. Bodur.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0236

https://doi.org/10.1287/ijoc.2023.0236.cd

Below is the BibTex for citing this snapshot of the repository.

```bibtex
@article{ACOTS_GuoNagarajanBodur2025,
  author =        {Cheng Guo, Harsha Nagarajan, Merve Bodur},
  publisher =     {INFORMS Journal on Computing},
  title =         {Tightening Quadratic Convex Relaxations for the AC Optimal Transmission Switching Problem},
  year =          {2025},
  doi =           {https://doi.org/10.1287/ijoc.2023.0236.cd},
  url =           {https://github.com/INFORMSJoC/2023.0236},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0236},
}  
```

## Description
The goal of this repository is to demonstrate the efficiency and effectiveness of our strengthened *Quadratic Convex (QC) relaxation* formulation and algorithms for the **AC Optimal Transmission Switching (ACOTS)** problem.

To showcase the performance of our proposed method, we conducted three sets of comparative experiments:  
- One against the **state-of-the-art QC relaxation** implemented in PowerModels.jl ("P" formulation),  
- Another against the **basic extreme-point enhanced formulation** ("E" formulation),  
- And a third using the **lifted cycle constraints and optimization-based bound tightening (OBBT)** ("ECB" and "ECB*" formulations).

The experiments were performed on a server running a 64-bit Linux operating system, equipped with an Intel Core i9-9900K CPU running at 3.6 GHz and 128 GB of RAM (for medium-scale instances), and a larger server with Intel CPUs and 250 GB of RAM (for large-scale instances). The test cases are based on the PGLib-OPF v20.07 benchmark library, covering instances from small networks (3 buses) to very large systems (up to 2,312 buses).

To ensure robustness in the presented results, each setup involved running the algorithms on a variety of test instances under different operating conditions (TYP, SAD, and API), with key metrics such as **optimality gap** and **runtime** recorded and compared.

In particular:
- We compared the basic PowerModels.jl relaxation ("P") against our *extreme-point enhanced* ("E") formulation, observing that "E" consistently provides tighter relaxations and closes larger optimality gaps, especially for congested (API) and small-angle (SAD) cases.
- We evaluated the impact of further strengthening through **lifted cycle constraints** and **bound tightening** ("ECB" formulation), which showed dramatic improvements, closing as much as **17.6%** more gap compared with "E" alone.
- We also implemented a **branch-and-cut variant** ("ECB\*") where lifted cycle constraints are dynamically generated during the branch-and-bound process. This approach retains nearly the same solution quality while significantly improving computational efficiency.
- When paired with the proposed maximum spanning tree heuristic, our solution methods find good-quality solutions for large-scale instances up to **2,312 buses**.

The experimental results confirm that our strengthened formulations not only tighten the relaxation bounds compared to previous benchmarks but also allow for faster and more scalable solutions to realistic, large-scale ACOTS problems.

## Replicating
Run the code and algorithm with the softwares and library of these versions:
- Julia - v1.6
- Juniper - v0.7.0
- Knitro - v13.0.1
- Ipopt - v3.13.4
- Gurobi - v9.0.0
- PGLib instances repo version v20.07 (provided in the [data](data) folder)

The Project.toml and Manifest.toml files specify package dependencies of the Julia code.

### Figure 2

### UB in Table 2 and Table 4
To obtain the values in the "UB" column of Table 2 and Table 4, we obtain the minimum of the following 4 values:
- Non-convex ACOTS model (1) with the Knitro solver: In "params.jl", set `params["ots_exact_pm"] = true`, then run the "test.jl" file.
- Non-convex ACOTS model (1) with the Juniper solver: In "params.jl", set `params["ots_exact_pm"] = true` and `params["minlp_solver"] = "juniper"`, then run the "test.jl" file.
- Non-convex ACOPF model with all lines switched on: In "params.jl", set `params["model"] = "opf"` and `params["opf_exact_pm"] = true`, then run the "test.jl" file.
- Non-convex ACOPF model with the set of lines switched off, as indicated by the ACOTS-QC solutions: In "params.jl", set `params["model"] = "opf"`, `params["opf_exact_pm"] = true`, and `params["turn_off_lines"] = true`, then run the "test.jl" file.

The result is available at [results/ub.xlsx](results/ub.xlsx)

### Other columns in Table 2 and Table 4
- Column "P": In line 54 of "tests.jl", set `args[4] = "rmc"` and `args[5] = true`, then run the "test.jl" file.
- Column "E": In line 54 of "tests.jl", set `args[4] = "tri"` and `args[5] = true`, then run the "test.jl" file. 
- Column "EC":
- Column "ECB":
- Column "ECB*":

Note that instances are run with warm-start files, which are located in [results/opf_sol](results/opf_sol). To generate warm-start files for a specific instances, ???.
