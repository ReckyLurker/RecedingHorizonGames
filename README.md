# RecedingHorizonGames
Contains the code for the Battery Charging Receding Horizon Game. 

Run `receding-horizon/src/simulation.m` from the project root to run the simulation. The Agent class contains VGNE computational logic, constraint tightening logic, terminal set, and steady state computation logic. The optimization is performed using [Sedumi-v1.3](https://github.com/sqlp/sedumi) and [MOSEK](https://www.mosek.com/). The vGNE computation involves solving the SOCP projection problem, which is solved using MOSEK. To switch to Sedumi, change the solver defined [here](https://github.com/ReckyLurker/RecedingHorizonGames/blob/main/receding-horizon/src/Agent.m#L325) to 'sedumi' in place of 'mosek'. 

If required, delete `libraries/mpt3` and run the `libraries/install_mpt3.m` script to reinstall MPT3. 
