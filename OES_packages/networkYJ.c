/*
	This code is a collection of fucntions
	for calculating the standard score of motifs.

	It was written by Youngjai in June 20, 2022.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
	struct __node *head;
	struct __node *tail;
} Node_list;

typedef struct __node {
	int node_name;
	double abundance;
	int degree;
	int *neighbors;
	int *direction; // incoming : -1, outgoing : 1
	double *sign;
	double self_loop;
	struct __node *next;
} Node;

typedef struct {
	int nodes;
	int links;
	int node_label;
	Node_list *node_list;
} Network;

typedef struct {
	int N0; // the number of initial nodes in the system
	double f_th; // the minimum abundance density
	double f0; // the initial abundnace density of each node (scaled by carrying capacity K)
	int m; // the number of initial interaction with the residents
	double sigma; // the standard deviation of interaction weight distribution
	double alpha; // the rate of a new node
	int mutant_event_time;
	int saturation_time;
	int simul_step;
	double simul_dt;
	unsigned long seed; // the random seed
} Parameters;

// from mt19937ar.c
void init_genrand(unsigned long s);
double genrand_real2(void);

// from youngjai_packages.c
double normal_dist(double avg, double stdev, unsigned long seed);
int *shuffle_array_int(int length, unsigned long seed);
void swap_two_integers(int *num_i, int *num_j);
void swap_two_floats(double *num_i, double *num_j);

void free_node_list(Node_list *&node_list)
{
	Node *cursor = node_list->head;
	Node *temp;

	while(cursor != NULL){
		temp = cursor;
		cursor = cursor->next;
		free(temp->neighbors);
		free(temp->direction);
		free(temp->sign);
		free(temp);
	}
	cursor = NULL; free(cursor);
	temp = NULL; free(temp);
	node_list = NULL; free(node_list);
}

void create_node(Network& network, Parameters& paras)
{
	Node *new_node = (Node *) calloc(1, sizeof(Node));
	new_node->node_name = network.node_label++;
	new_node->abundance = paras.f0;
	new_node->degree = 0;
	new_node->neighbors = (int *) calloc((6*paras.m), sizeof(int));
	new_node->direction = (int *) calloc((6*paras.m), sizeof(int));
	new_node->sign = (double *) calloc((6*paras.m), sizeof(double));
	// new_node->self_loop = normal_dist(0., paras.sigma, paras.seed++);
	new_node->self_loop = 1.;
	// new_node->self_loop = paras.sigma;
	new_node->next = NULL;

	if(network.node_list->head==NULL && network.node_list->tail==NULL)
		network.node_list->head = network.node_list->tail = new_node;
	else{
		network.node_list->tail->next = new_node;
		network.node_list->tail = new_node;
	}
	network.nodes++;
}

void create_connection(Network& network, Parameters& paras, int name_i, int name_j)
{
	int i;
	Node *cursor = network.node_list->head;
	Node *node_i, *node_j;
	int count = 0;
	double weight_temp = normal_dist(0., paras.sigma, paras.seed++);
	int direct_temp;
	if(genrand_real2() < 0.5) direct_temp = -1;		// Do not allow multiedges
	else direct_temp = 1;

	while(cursor->next != NULL){
		cursor = cursor->next;
		if(cursor->node_name == name_i){
			node_i = cursor;
			count++;
		}
		if(cursor->node_name == name_j){
			node_j = cursor;
			count++;
		}
		if(count == 2) break;
	}

	count = 0;
	for(i=0; i<node_i->degree; i++) if(node_i->neighbors[i] == name_j) count++;

	if(count < 1){
		if((node_i->degree+1)%(6*paras.m) == 0){
			node_i->neighbors = (int *) realloc(node_i->neighbors, (node_i->degree+(6*paras.m))*sizeof(int));
			node_i->direction = (int *) realloc(node_i->direction, (node_i->degree+(6*paras.m))*sizeof(int));
			node_i->sign = (double *) realloc(node_i->sign, (node_i->degree+(6*paras.m))*sizeof(double));
		}
		if((node_j->degree+1)%(6*paras.m) == 0){
			node_j->neighbors = (int *) realloc(node_j->neighbors, (node_j->degree+(6*paras.m))*sizeof(int));
			node_j->direction = (int *) realloc(node_j->direction, (node_j->degree+(6*paras.m))*sizeof(int));
			node_j->sign = (double *) realloc(node_j->sign, (node_j->degree+(6*paras.m))*sizeof(double));			
		}
		node_i->neighbors[node_i->degree] = name_j;
		node_j->neighbors[node_j->degree] = name_i;

		node_i->sign[node_i->degree] = weight_temp;
		node_j->sign[node_j->degree] = weight_temp;

		node_i->direction[node_i->degree++] = -direct_temp;
		node_j->direction[node_j->degree++] = direct_temp;
		network.links++;
	}
	else if(count < 2){
		if((node_i->degree+1)%(6*paras.m) == 0){
			node_i->neighbors = (int *) realloc(node_i->neighbors, (node_i->degree+(6*paras.m))*sizeof(int));
			node_i->direction = (int *) realloc(node_i->direction, (node_i->degree+(6*paras.m))*sizeof(int));
			node_i->sign = (double *) realloc(node_i->sign, (node_i->degree+(6*paras.m))*sizeof(double));
		}
		if((node_j->degree+1)%(6*paras.m) == 0){
			node_j->neighbors = (int *) realloc(node_j->neighbors, (node_j->degree+(6*paras.m))*sizeof(int));
			node_j->direction = (int *) realloc(node_j->direction, (node_j->degree+(6*paras.m))*sizeof(int));
			node_j->sign = (double *) realloc(node_j->sign, (node_j->degree+(6*paras.m))*sizeof(double));			
		}
		node_i->neighbors[node_i->degree] = name_j;
		node_j->neighbors[node_j->degree] = name_i;

		node_i->sign[node_i->degree] = weight_temp;
		node_j->sign[node_j->degree] = weight_temp;

		for(i=0; i<node_i->degree; i++) if(node_i->neighbors[i] == name_j) break;		
		node_i->direction[node_i->degree++] = -1*node_i->direction[i];
		node_j->direction[node_j->degree++] = 1*node_i->direction[i];
		network.links++;
	}
	cursor = NULL; free(cursor);
	node_i = NULL; free(node_i);
	node_j = NULL; free(node_j);
}

void initial_SPP_network(Network& network, Parameters& paras)
{
	int i;
	int name_i, name_j;
	int m_temp;
	int *shuffled = shuffle_array_int(paras.N0*paras.N0, paras.seed++);

	//Node_list pointer define start
	network.nodes = -1;
	network.links = 0;
	network.node_label = 0;
	network.node_list = (Node_list *) calloc(1, sizeof(Node_list));
	network.node_list->head = NULL; network.node_list->tail = NULL;
	create_node(network, paras); // define head pointer
	//Node_list pointer define end

	/* incubation rule. */
	if(paras.N0 <= paras.m) m_temp = paras.N0-1;
	else m_temp = paras.m;

	for(i=0; i<paras.N0; i++) create_node(network, paras);

	i = 0;
	while(network.links < paras.N0*m_temp/2){
		name_i = shuffled[i]/paras.N0%paras.N0;
		name_j = shuffled[i++]%paras.N0;

		if(name_i != name_j) create_connection(network, paras, ++name_i, ++name_j);
	}
	free(shuffled);
}

double *RK4(double *f, double t, double h, Network network, double *(*dfdt)(double *, double, Network))
{
	int i;
	double h2 = 0.5*h;

	double **k = (double **) calloc(4, sizeof(double *));

	// RK4
	// k1 = f(x,t)
	// k2 = f(x+h2*k1, t+h2)
	// k3 = f(x+h2*k2, t+h2)
	// k4 = f(x+h*k3, t+h)
	// x += h/6 * (k1 + 2k2 + 2k3 + k4)
	double *f_rk4 = (double *) calloc(network.nodes, sizeof(double));

	// k1
	k[0] = dfdt(f, t, network);

	// k2
	for(i=0; i<network.nodes; i++) f_rk4[i] = f[i]+h2*k[0][i];
	k[1] = dfdt(f_rk4, t+h2, network);

	// k3
	for(i=0; i<network.nodes; i++) f_rk4[i] = f[i]+h2*k[1][i];
	k[2] = dfdt(f_rk4, t+h2, network);

	// k4
	for(i=0; i<network.nodes; i++) f_rk4[i] = f[i]+h*k[2][i];
	k[3] = dfdt(f_rk4, t+h, network);

	// RK4 array
	for(i=0; i<network.nodes; i++)
		f_rk4[i] = f[i] + h/6.*(k[0][i]+2.*k[1][i]+2.*k[2][i]+k[3][i]);

	for(i=0; i<4; i++) free(k[i]);
	free(k);

	return f_rk4;
}

double *Adaptive_RK4(double *f, double t, double &h, Network network, double *(*dfdt)(double *, double, Network))
{
	int i;
	double *f1 = NULL;
	double *f_tmp, *f2;
	double f_del;
	double rho = 0.;
	// If the difference between two results is less than 'precision', these are same.
	double acc = 1e-8; // target accuracy
	double pre = 1e-14; // 소수 14번째 자리까지 정도 확인

	while(rho<1.-pre){
		f_del = 0.; 

		if(f1 != NULL){ delete f1; delete f_tmp; delete f2; }
		f_tmp = RK4(f, t, h, network, dfdt);
		f1 = RK4(f_tmp, t+h, h, network, dfdt);
		f2 = RK4(f, t+2.*h, 2.*h, network, dfdt);

		// To measure the Euclidean distance in f dynamics space
		for(i=0; i<network.nodes; i++) f_del += (f1[i]-f2[i])*(f1[i]-f2[i]);
		f_del = sqrt(f_del);

		// If f_del is too small, we replace f_del as precision.
		if(f_del < pre){ 
			rho = 10.; h *= 2.;
		}
		else{
			// If f_del is diverse, we replace f_del as one over precision.
			if((isnan(f_del)) || (isinf(f_del))) f_del = 1./pre;
			rho = 30.*h*acc/f_del;

			if(rho<1.-pre) h *= pow(rho, 0.25);
			else h *= 2.;
		}
	}

	// adaptive RK4 array
	double *f_adap = (double *) calloc(network.nodes, sizeof(double));
	for(i=0; i<network.nodes; i++) f_adap[i] = f1[i] + (f1[i]-f2[i])/15.;

	free(f1); free(f_tmp); free(f2);

	return f_adap;
}
