
# this is a comparison of the results of seismic_CPML_2D_velocity_and_stress_fourth_order_viscoelastic.f90
# and of the analytical solution of analytical_solution_viscoelastic_2D_plane_strain_Carcione_correct_with_1_over_L.f90  seismic_CPML_2D_velocity_and_stress_fourth_order_viscoelastic.f90

set term x11

set xrange [0:1.2]

plot "Vx_file_001.dat" w l lc 1, "Vx_time_analytical_solution_viscoelastic.dat" w l lc 3
pause -1 "Hit any key..."

plot "Vy_file_half_a_grid_cell_away_from_Vx_001.dat" w l lc 1, "Vz_time_analytical_solution_viscoelastic.dat" w l lc 3
pause -1 "Hit any key..."

