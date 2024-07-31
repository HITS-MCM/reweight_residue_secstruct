# reweight_residue_secstruct
Reweighting of GaMD (Amber)-sampled residue secondary structure propensities of peptides.

Please see the following publication for more information on the applied discrete residue-based secondary structure propensity reweighting approach: 

    "Simulation of the Positive Inotropic Peptide S100A1ct in Aqueous Environment by Gaussian Accelerated Molecular Dynamics"<br/> 
    M. Glaser, N. J. Bruce, S. B. Han, R. C. Wade, J Phys Chem B. 2021 May 13;125(18):4654-4666 (https://pubs.acs.org/doi/10.1021/acs.jpcb.1c00902)

The reweighting procedure is implemented in the Python module 'src/reweight_residue_secstruct/reweight_residue_secstruct.py'.<br/>
Please also check the introductory comment in 'reweight_residue_secstruct.py' for more information.<br/>

To run 'reweight_residue_secstruct.py', please first created the necessary "environment", either by creating

1) the Conda environment 'reweight_residue_secstruct' via the supplied YAML file (environment_reweight_residue_secstruct.yml)

    $ conda env create -f environment_reweight_residue_secstruct.yml<br/>
    $ conda activate reweight_residue_secstruct<br/>

2) or alternatively, by creating the corresponding Conda environment via the following set of commands

    $ conda create -n reweight_residue_secstruct python=3.6<br/>
    $ conda activate reweight_residue_secstruct<br/>
    $ conda install -c conda-forge pandas<br/>
    $ conda install -c conda-forge matplotlib<br/>
    $ conda install -c numba numba<br/>

3) or by use of any other means to setup Python and the corresponding libraries.

An example case on how to (i) create the input files for 'reweight_residue_secstruct.py' from examplary Amber GaMD output files and to (ii) run 'reweight_residue_secstruct.py' is provided in the directory 'example'. The directory contains examplary Amber GaMD production run trajectory files (gaMD_params_stripped.nc, gaMD_prod_1_stripped_cut.nc) as well as the corresponding Amber GaMD log files (gamd.log, gaMD_prod_1_cut.log), containing the GaMD boost potential information for the respective MD steps the frames were written. The examplary files were generated with Amber18. Furthermore, 'example' contains two bash scripts: prepare_input_reweight_residue_secstruct.bsh and run_reweight_residue_secstruct.bsh. Please note that the example bash scripts as well as 'reweight_residue_secstruct.py' were only used in combination with GaMD output from Amber18.

To rerun the example, just remove the directories 'input_reweight_residue_secstruct' and 'output_reweight_residue_secstruct' and make sure you have a working Amber installation.

Then execute the following commands within the 'example' directory:

1) $ ./prepare_input_reweight_residue_secstruct.bsh 

    This bash script will generate the necessary input files for 'reweight_residue_secstruct.py', 'dssp.out' and 'weights.dat', from the GaMD trajectory and GaMD log files, respectively.
    The bash script will write its output to the directory 'input_reweight_residue_secstruct'.
    The file 'dssp.out' contains the secondary structure information for each peptide residue of the all the GaMD frames, obtained via cpptraj.
    The file 'weights.dat' contains the corresponding boost potential information parsed from the GaMD log files in the format that is also used by the standard GaMD reweighting scripts "PyReweighting-*D.py", used for continuous order parameters.
    Please consult the comments in 'prepare_input_reweight_residue_secstruct.bsh' for more information on its input and the generated output etc.
    The script 'prepare_input_reweight_residue_secstruct.bsh' aggregates the information from a consecutive GaMD production run (i.e., a single replicon), split over two trajectory files.
    In case of the first trajectory file (i.e., gaMD_params_stripped.nc), only the last third of frames is taken (as the remaining part corresponds to preparatory steps).
    If your trajectory is split over more than two trajectory files, you have to extend the script accordingly.
    By making the corresponding adaptions to 'prepare_input_reweight_residue_secstruct.bsh', in principle, one can also aggregate information from different GaMD replica for reweighting.
    In any case, make sure that the order of your GaMD frames and their respective GaMD log file information are in sync.

3) $ ./run_reweight_residue_secstruct.bsh 

    This is a wrapper-script for 'reweight_residue_secstruct.py'.
    It will collect all the output from 'reweight_residue_secstruct.py' in the directory 'output_reweight_residue_secstruct'.
    Please check comments in 'run_reweight_residue_secstruct.bsh' for information on the 'reweight_residue_secstruct.py' output files.
    The script 'reweight_residue_secstruct.py' can be run with automatic parallelization (@njit) by using its '-m parallel' option.

Note that the example bash scripts were setup in the context of the above mentioned 'reweight_residue_secstruct' Conda environment.
If you want to run the example yourself, you need to properly setup your Amber and Conda installations by changing the respective lines in the bash scripts (indicated by the comment '# PLEASE ADAPT THIS!').
