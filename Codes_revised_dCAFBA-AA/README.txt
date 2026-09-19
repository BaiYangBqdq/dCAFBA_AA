%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This package is provided as Supporting Material software to the article:

A global resource constrained model for predicting metabolic flux dynamics during amino acid shifts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Contents:
1) DESCRIPTION
2) PREREQUISITES
3) dCAFBA-AA WORKFLOW
4）FVA analysis


=== 1: DESCRIPTION ===

This package provides COBRA-compatible MATLAB code for dCAFBA-AA simulations
of amino-acid (AA) shifts in Escherichia coli.

The package contains the following simulation cases:

1. AA downshift, glycerol, iJR904
       dCAFBA_AADown_Glyc_iJR904.m
       dCAFBA_AADown_Glyc_iJR904_pre.m
       dCAFBA_AADown_Glyc_iJR904_post.m

2. AA downshift, glucose, iJR904
       dCAFBA_AADown_Glc_iJR904.m
       dCAFBA_AADown_Glc_iJR904_pre.m
       dCAFBA_AADown_Glc_iJR904_post.m

3. AA downshift, glycerol, iML1515
       dCAFBA_AADown_Glyc_iML1515.m
       dCAFBA_AADown_Glyc_iML1515_pre.m
       dCAFBA_AADown_Glyc_iML1515_post.m

4. AA upshift, glycerol, iJR904
       dCAFBA_AAUp_Glyc_iJR904.m
       dCAFBA_AAUp_Glyc_iJR904_pre.m
       dCAFBA_AAUp_Glyc_iJR904_post.m

For each case, the main script performs the dynamic simulation, while the
_pre and _post functions solve the corresponding time-specific dCAFBA
problems before and after the nutrient shift.

The package also contains the following main utility functions:

CAFBA_OptimizeCbModel_glpk.m
   Solves the CAFBA linear programming problem using GLPK.


setWeights.m
   Adds or modifies the proteome allocation weights.

CAFBA_buildLP.m
   Constructs the CAFBA linear programming problem.

CAFBA_FVA_glpk.m
   Performs  flux variability analysis (FVA) with the proteome allcoation constraints.


=== 2: PREREQUISITES ===

dCAFBA-AA uses metabolic models in the COBRA MATLAB structure and is based
on the CAFBA framework of Mori et al. (2016).

The main requirements are:

    MATLAB
    COBRA Toolbox
    GLPK


The corresponding metabolic-model and experimental-data MAT files required
for each simulation are loaded in the individual main scripts.


=== 3: dCAFBA-AA WORKFLOW ===

The general workflow is as follows:

1. Load the genome-scale metabolic model and define the pre-shift and
   post-shift growth conditions.

2. Set up the model parameters, amino-acid biosynthesis groups, and initial
   proteome fractions.

3. Initialize the cellular proteome state and simulate dynamic proteome
   reallocation during the AA shift.

4. At each time point, update the proteome-dependent CAFBA constraints and
   solve the corresponding linear programming problem.

5. Store the resulting metabolic fluxes and cellular-state variables for
   subsequent analysis.


Detailed model formulation, parameterization, dynamic equations, and
simulation procedures are described in the accompanying manuscript.

=== 4: Perform FVA ===
The time-specific CAFBA models saved during the dCAFBA-AA simulation can be
used for  flux variability analysis (FVA).

FVA is performed using:

    cafvaResult = CAFBA_FVA_glpk(model_t, ...);

Here, model_t is the time-specific CAFBA model saved during the simulation.


See the accompanying manuscript for details of the CAFVA procedure and
parameter settings.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
