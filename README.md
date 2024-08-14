# reweight_residue_secstruct
## Reweighting of GaMD (Amber)-sampled residue secondary structure propensities of peptides

>[!CAUTION]
>**Beta version/release candidate!**

Please see the following publication for more information on the applied **discrete residue-based secondary structure propensity reweighting** approach: 

> "Simulation of the Positive Inotropic Peptide S100A1ct in Aqueous Environment by Gaussian Accelerated Molecular Dynamics"; M. Glaser, N. J. Bruce, S. B. Han, R. C. Wade, J Phys Chem B. 2021 May 13;125(18):4654-4666 (https://pubs.acs.org/doi/10.1021/acs.jpcb.1c00902)

The reweighting procedure is implemented in the Python module **`src/reweight_residue_secstruct/reweight_residue_secstruct.py`** (please also check the introductory comment in `reweight_residue_secstruct.py` and its options for more information).

### Python libraries needed for `reweight_residue_secstruct.py`

To run `reweight_residue_secstruct.py`, please first create the necessary Conda environment `reweight_residue_secstruct` via the supplied YAML file (`environment_reweight_residue_secstruct.yml`)
```
$ conda env create -f environment_reweight_residue_secstruct.yml
```

### Example case on how to use `reweight_residue_secstruct.py`

The example case provided in the directory `example` allows: 
**i.** to create the input files for `reweight_residue_secstruct.py` from Amber GaMD output files and 
**ii.** to run `reweight_residue_secstruct.py`.

**The `example` directory contains:** 
* Amber GaMD production run trajectory files (`gaMD_params_stripped.nc`, `gaMD_prod_1_stripped_cut.nc`)
* the corresponding Amber GaMD log files (`gamd.log`, `gaMD_prod_1_cut.log`), containing the GaMD boost potential information for the respective MD steps the frames were outputted at
* two bash scripts: `prepare_input_reweight_residue_secstruct.bsh` and `run_reweight_residue_secstruct.bsh`

> [!CAUTION]
> The example Amber GaMD output files were generated with **Amber18**. 
> Please note that the example bash scripts as well as `reweight_residue_secstruct.py` were **only used in combination with GaMD output from Amber18**.

> [!IMPORTANT]
> To rerun the example, just remove the directories `input_reweight_residue_secstruct` and `output_reweight_residue_secstruct`.
> Also **make sure you have a working Amber installation**, as `cpptraj` is needed.
> The example bash scripts were setup and tested in the context of the above mentioned `reweight_residue_secstruct` Conda environment as well as `cpptraj` from Amber18.

**To run the example, do the following:**

*1) Setup your Amber installation, if not already done:*

```
$ source /your-path/amber18/amber.sh 
```

*2) Activate the Conda environment `reweight_residue_secstruct`:* 

```
$ conda activate reweight_residue_secstruct
```

*3) Next, execute the following two bash scripts within the `example` directory:*

```
$ ./prepare_input_reweight_residue_secstruct.bsh
$ ./run_reweight_residue_secstruct.bsh
```

**About `prepare_input_reweight_residue_secstruct.bsh`:**

This bash script will generate the necessary input files for `reweight_residue_secstruct.py`, `dssp.out` and `weights.dat`, from the GaMD trajectory and GaMD log files, respectively.
The bash script will write its output to the directory `input_reweight_residue_secstruct`.
The file `dssp.out` contains the secondary structure information for each peptide residue of the all the GaMD frames, obtained via **`cpptraj`**.
The file `weights.dat` contains the corresponding boost potential information parsed from the GaMD log files **in the format that is also used by the standard GaMD reweighting scripts `PyReweighting-*D.py`**, which are used for reweighting of continuous order parameters.

> [!NOTE]
> Please consult the comments in `prepare_input_reweight_residue_secstruct.bsh` for more information on its input and the generated output etc.

In the example case, the script `prepare_input_reweight_residue_secstruct.bsh` aggregates the information from a consecutive GaMD production run (i.e., a single replicon), split over two trajectory files.
In case of the first trajectory file (i.e., `gaMD_params_stripped.nc`), only the last third of frames is taken (as the remaining part corresponds to preparatory steps).

> [!IMPORTANT]
> * If your trajectory is split over more than two trajectory files, you have to extend the script accordingly, in any case, **make sure that the order of your GaMD frames and their respective GaMD log file information are in sync.**
> * A **dual boost potential** was applied during GaMD simulations (relevant for parsing of GaMD boost potential information to generate the file `weights.dat`).

> [!TIP]
> By making the corresponding adaptions to `prepare_input_reweight_residue_secstruct.bsh`, in principle, one can also aggregate information from different GaMD replica for reweighting.

**About `run_reweight_residue_secstruct.bsh`:**

The bash script is a wrapper-script to run `reweight_residue_secstruct.py`.
It will collect all the output from `reweight_residue_secstruct.py` in the directory `output_reweight_residue_secstruct`.

**About `reweight_residue_secstruct.py`:**

> [!NOTE]
> * Please check comments in `run_reweight_residue_secstruct.bsh` for information on the `reweight_residue_secstruct.py` output files.
> * The script `reweight_residue_secstruct.py` checks for common standard Amber amino acid three-letter codes, if your residue is not listed there, the `allowed_residues`-list needs to be extended accordingly. 

> [!TIP]
> The script `reweight_residue_secstruct.py` can be run with automatic parallelization (`@njit`) by using its `-m parallel` option.
> With this option, depending on the resources, as many residue secondary structure reweightings as possible are then performed simultaneously.
