gfx read node Micropipette.part0.exnode
gfx read elem Micropipette.part0.exelem
gfx define faces

gfx define field deformed composite Deformation.1 Deformation.2

gfx define field F gradient coordinate Coordinate field deformed
gfx define field F_transpose transpose source_number_of_rows 2 field F
gfx define field Identity2 composite 1 0 0 1
gfx define field C matrix_multiply number_of_rows 2 fields F_transpose F
gfx define field E2 add fields C Identity2 scale_factors 1 -1
gfx define field E scale field E2 scale_factors 0.5 0.5 0.5 0.5
gfx define field principal_strains eigenvalues field E
gfx define field principal_strain_vectors coordinate_system rectangular_cartesian eigenvectors eigenvalues principal_strains
gfx define field deformed_principal_strain_vectors coordinate_system rectangular_cartesian matrix_multiply number_of_rows 3 fields principal_strain_vectors F_transpose

gfx define field E11 composite E.1
gfx define field E12 composite E.2
gfx define field E21 composite E.3
gfx define field E22 composite E.4

