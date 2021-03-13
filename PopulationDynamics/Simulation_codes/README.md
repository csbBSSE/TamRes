# Algorithm for population dynamics
- The algorithm for simulating population dynamics is as follows- 

  1. Fix the values of Heterogenity, Plasticity (P_SR and P_RS), proliferation rate, death rate
  2. Initialize the simulation with a list of cells (for eg. initialize the simulation with 100 sensitive cells, or 50 sensitive and 50 resistant cells). This step involves randomly choosing scores of cells from a given gaussian with mean -2 or +2 and given heterogeneity.
  3. Choose a cell from this list.
  4. Randomly decide the order of proliferation, death and switching.
  5. Generate random numbers and check if the criteria for proliferation, death or switching is satified.
  6. If yes, carry out that process. If no, then continue.
  7. For proliferation, create a new cell and sample its score from the distribution.
  8. For switching, flip the `current_score_of_cell` parameter. This is 1 for sensitive and 2 for resistant. 
  9. For death, remove the cell from the list of cells and break the loop. 
  10. Repeat process 3-9 for each cells in the list of cells.
  11. Repeat process 3-10 for each time step. 
  12. Calculate the diversity of the population, population size, number of switches, number of deaths for each timestep.
  
- For the case of Drug induced plasticity - 
  - If a sensitive cell does not die due to drug, then with some fixed probability switch it to resistant cell.
  
- In the presence of MET - 
  - If a resistant cell does not die, then with some fixed probability switch it to sensitive cell.