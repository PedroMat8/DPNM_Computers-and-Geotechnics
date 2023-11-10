# Please cite the following paper: 
"Drying-induced volumetric behaviour of clays interpreted via binary pore-scale modellings"

Matteo Pedrotti, Long Xu, Ian Murray, Alessandro Tarantino

Computers and Geotechnics 2023

https://doi.org/10.1016/j.compgeo.2023.105911
 
# DPNM
Deformable pore network model

## pn_classes.py
File containing all the pore network classes:
 - Soil: soil parameters
 - Fluid: wetting fluid parameters
 - Test: drying test parameters
 - Macro: macro variable
 - PN: pore network
 - MainFigure

## pn_functions.py
  - PN loop:
    - cycle: iteration over the pn
    - coordination_count: counts adjacient pores saturation
  - Import/Export:
    - save_psd
    - save_pickle
    - import_pickle

## pn_generate_pn.py --> all parameters are here
  - get_pn: initialize pore network classes

## pn_run_pn.py
  - run_pn: import pore network classes and run

## pn_main.py

For additional information please consult the paper or contact matteo.pedrotti@strath.ac.uk
