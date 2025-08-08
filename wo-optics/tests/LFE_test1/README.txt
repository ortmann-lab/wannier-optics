In this test we setup a matrix where we specify the single entries by local field effects.
The LFE are only specified in the unit cell (0,0,0) but the supercell is larger.
A onsite energy of 2eV is added to the diagonals to make sure we do not have zero eigenvalues
(which are for some reasons a problem for the Lanczos algorithm?!)
Single particle TIs are not used

The matrix can be solved with numpy (see peak_analysis script) and with the exciton code.
