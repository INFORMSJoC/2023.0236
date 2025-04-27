# Tightening Quadratic Convex Relaxations for the AC Optimal Transmission Switching Problem (2023.0236)

This archive is distributed in association with the INFORMS Journal on Computing under the MIT License.

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper "Tightening Quadratic Convex Relaxations for the AC Optimal Transmission Switching Problem" by C. Guo, H. Nagarajan and M. Bodur.

## Cite 

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

*** Enter DOI link for paper ***

*** Enter DOI link for code ***

Below is the BibTex for citing this snapshot of the repository.

```bibtex
@article{ACOTS_GuoNagarajanBodur2025,
  author =        {Cheng Guo, Harsha Nagarajan, Merve Bodur},
  publisher =     {INFORMS Journal on Computing},
  title =         {Tightening Quadratic Convex Relaxations for the AC Optimal Transmission Switching Problem},
  year =          {2025},
  doi =           {},
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
- We evaluated the impact of further strengthening through **lifted cycle constraints** and **bound tightening** ("ECB" formulation), which showed dramatic improvements, often closing over **17%** more gap compared to "E" alone.
- We also implemented a **branch-and-cut variant** ("ECB*") where lifted cycle constraints are dynamically generated during the branch-and-bound process. This approach retains nearly the same solution quality while significantly improving computational efficiency, enabling the solution of large-scale instances up to **2,312 buses**.

The experimental results confirm that our strengthened formulations not only tighten the relaxation bounds compared to previous benchmarks but also allow for faster and more scalable solutions to realistic, large-scale ACOTS problems.

## Running
Run the code and algorithm with the packages of these versions: 
- Julia - vXXX 
- JuMP - vXXX 
- Gurobi - vXXX 
- PGLib instances repo version
