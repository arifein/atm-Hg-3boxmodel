Three-box model of atmospheric mercury (Hg) cycle in Python. The model includes two boxes for each hemisphere in the troposphere, and one stratospheric box. Elemental mercury (Hg0) and oxidized mercury (Hg2+) emissions, transformations, transport, and deposition are simulated. Rate coefficients of the model are tuned from GEOS-Chem results. The model can be used to tune photo-reduction of Hg2+ in GEOS-Chem, by identifying the necessary reduction rate needed to balance the Hg0 burden. 

Files:

atm_3box_model_Hg.py - Code to initialize, run, and visualize atmospheric Hg box model

The code builds off of an earlier Hg box model Python implementation from Noelle Selin, https://github.com/noelleselin/sixboxmercury.

