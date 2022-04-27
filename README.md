Three-box model of atmospheric mercury (Hg) cycle in Python. The model includes two boxes for each hemisphere in the troposphere, and one stratospheric box. Elemental mercury (Hg0) and oxidized mercury (Hg2+) emissions, transformations, transport, and deposition are simulated. Rate coefficients of the model are tuned from GEOS-Chem results. The model can be used to tune photo-reduction of Hg2+ in GEOS-Chem, by identifying the necessary reduction rate needed to balance the Hg0 burden. 

Files:

atm_3box_model_Hg.py - Code to initialize, run, and visualize atmospheric Hg box model

Originally described in Feinberg et al. (2022): Evaluating atmospheric mercury (Hg) uptake by vegetation in a chemistry-transport model, Environmental Science: Processes & Impacts, https://doi.org/10.1039/D2EM00032F. The code builds off of an earlier Hg box model Python implementation from Noelle Selin, https://github.com/noelleselin/sixboxmercury.

[![DOI](https://zenodo.org/badge/451642669.svg)](https://zenodo.org/badge/latestdoi/451642669)
