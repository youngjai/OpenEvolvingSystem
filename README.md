# Open evolving system model
&nbsp;An ecosystem consists of various species.
The system has two features: the one is species living in the system interact with each other.
The other is new species invade the system constantly.
As a result, the population composition of the system changes over time.
To examine how the population composition of the system changes depending on these features, we suggest a network model with a generalized Lotka-Volterra type equation in an Open Evolving System.
In the model, a new species appears in the system at an interval time.
The dynamics of our model follows as the Lotka-Volterra type equation.
The equation is given by

![equation](https://latex.codecogs.com/gif.image?\dpi{110}\frac{df_i}{dt}&space;=&space;G_i(\textbf{f})f_i\left(1-\frac{\sum&space;f_j}{K}\right)&plus;D_i(\textbf{f})f_i .)

where x_i is an abundance of ith species, G_i and D_i are the growth and death rates depending on the interaction structure of ith species, respectively.
Especially, the species living in the model share resources which are limited such as habitat and water resources, so-called Resource Limitation, K.

# How to run the code
To run our code, we need six input parameters. 
```bash
g++ -g OES_model.cpp OES_packages/*.c -o OES_model.out
./OES_model.out <date(yymmdd)> <sigma> <mutant event time> <saturation time> <dt or K> <alpha>
```
\<date(yymmdd)\> is a folder name that will be saved output files.
\<sigma\> and \<alpha\> mean the interaction strength among species and the invasion rate of new species.
\<mutant event time\> and \<saturation time\> are related to our simulation time. 
The former indicates how many species invade the system. 
The latter means how long the code runs more without the invasion after the invasion finishes.
\<dt or K\> is a unit time step of the code.

# The output files
The code generates four types of output files when to finish.
The output files are written as below
- 'log.dat': parameter set and run-time
- 'simul_log.txt': normalized abundance of all species over time
- 'network_edges.txt': information of the network structure(source id, target id, link id, weight)
- 'num_of_S.txt': the number of species over time
