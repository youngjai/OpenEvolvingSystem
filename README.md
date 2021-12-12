# Open evolving system model
&nbsp;An ecosystem consists of various species.
The system has two features: the one is species living in the system interact with each other.
The other is new species invade the system constantly.
As a result, the population composition of the system changes over time.
To examine how the population composition of the system changes depending on these features, we suggest a network model with a generalized Lotka-Volterra type equation in an Open Evolving System.
In the model, a new species appears in the system at an interval time.
The dynamics of our model follows as the Lotka-Volterra type equation.
The equation is given by

![equation](https://latex.codecogs.com/gif.latex?%5Cfrac%7Bd%7Bx_i%7D%7D%7Bdt%7D%3DG_i%28%5Cmathbf%7Bx%7D%29%5C%3Ax_i%5Cleft%281-%5CSigma_jx_j/K%5Cright%29&plus;D_i%28%5Cmathbf%7Bx%7D%29%5C%3Ax_i%20%5C%3B%5C%3B.)

where x_i is an abundance of ith species, G_i and D_i are the growth and death rates depending on the interaction structure of ith species, respectively.
Especially, the species living in the model share resources which are limited such as habitat and water resources, so-called Resource Limitation, K.

# How to run the code
To run our code, we need six input parameters. 
```bash
g++ -g OES_model.cpp OES_packages/*.c -o OES_model.out
./OES_model.out <date(yymmdd)> <sigma> <mutant event time> <saturation time> <dt or K> <alpha>
```
\<date(yymmdd)\> is a folder name that will be saved output files.
\<sigma\> and \<alpha\> mean the interaction stregnth among species and the invasion rate of new species.
<mutant event time> and <saturation time> are related to our simulation time. 
The former indicates how many species invade the system. 
The latter means how long the code runs more without the invasion after the invasion finishes.
<dt or K> is an unit time step of the code.
d

# The output files
dfdfd
