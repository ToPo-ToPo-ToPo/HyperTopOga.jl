# HyperTopOga.jl
## Description
**HyperTopOga.jl** (**H**yperelastic **T**opology **O**ptimization using **G**radient-based Methods and **A**utomatic Differentiation) is a topology optimization library that accounts for hyperelastic material models. By applying automatic differentiation to part of the sensitivity analysis, it enables topology optimization regardless of the choice of objective function or hyperelastic model.
This library is based on the following paper:

**Sensitivity analysis of any hyperelastic evaluation functions coupled with adjoint method and automatic differentiation, Finite Elements in Analysis and Design
Vol. 252 (2025), 104440**. [DOI:10.1016/j.finel.2025.104440](https://doi.org/10.1016/j.finel.2025.104440)

When using this library, we would appreciate it if you could cite this paper and give it a star.
If you would like to make improvements, please feel free to fork it if necessary.

## Installation
1. Set up a virtual environment.
    ```
    (@v1.10) pkg> activate .
    ```
2. Install the dependent packages.
    ```
    (HyperTopOga) pkg> instantiate
    ```

## Usage
1. Set up a virtual environment.
    ```
    (@v1.10) pkg> activate .
    ```
2. Run the program from the command. hoge.jl is the main program you want to run. Change it as needed. A is the number of multi-threaded parallelism.
    ```
    julia -t A --project=. hoge.jl
    ```

## Citation
### Paper Reference
[DOI:10.1016/j.finel.2025.104440](https://doi.org/10.1016/j.finel.2025.104440)


## License
This project is licensed under the MIT License.  
See the [LICENSE](./LICENSE.txt) file for details.
