/*
	This code is a collection of fucntions
	for calculating the standard score of motifs.

	It was written by Youngjai in Oct. 12, 2020.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

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
	double y_th; // the minimum abundance density
	double y0; // the initial abundnace density of each node (scaled by carrying capacity K)
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

clock_t time_begin(char *path)
{	
	clock_t begin = clock();
	FILE *log = fopen(path, "wt");
	fprintf(log, "Log file is opened by youngjai.\r\n");
	fclose(log);

	return begin;
}

void time_end(clock_t begin, char *path)
{
	FILE *log = fopen(path, "at");
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC; /* unit : second */
	int dd, hh, mm;
	double ss;
	dd = time_spent/86400;
	hh = (time_spent-dd*86400)/3600;
	mm = (time_spent-dd*86400-hh*3600)/60;
	ss = time_spent-dd*86400-hh*3600-mm*60;

	fprintf(log, "\r\nThe runing time is %d days %02d hours %02d minutes %02.4lf seconds.\r\n", dd, hh, mm, ss);
	fclose(log);
}

void swap_two_integers(int *num_i, int *num_j)
{
	int temp = *num_i;	/* temporary storage */
	*num_i = *num_j; 
	*num_j = temp;
}

void swap_two_floats(double *num_i, double *num_j)
{
	double temp = *num_i;	/* temporary storage */
	*num_i = *num_j; 
	*num_j = temp;
}

int call = 0;
double normal_dist(double avg, double stdev, unsigned long seed)
{
	init_genrand(seed);
	/* generate a number from the normal distribution. */
	double v1, v2, s, temp;
	static double x1, x2;
	// static int call = 0;

	if(call == 1){
		call = !call;
		return (avg + stdev * (double) x2);
	}

	do{
		v1 = 2. * genrand_real2() - 1.;
		v2 = 2. * genrand_real2() - 1.;
		s = v1 * v1 + v2 * v2;
	} while(s >= 1. || s == 0.);

	temp = sqrt((-2. * log(s)) / s);
	x1 = v1 * temp;
	x2 = v2 * temp;

	call = !call;

	return (avg + stdev * (double) x1);
}

int *shuffle_array_int(int length, unsigned long seed)
{
	init_genrand(seed);
	int i, j;
	int *shuffle = (int *) calloc(length+1, sizeof(int));
	for(i=0; i<length; i++) shuffle[i] = i;

	for(i=length-1; i>0; i--){
		j = (int) ((i+1)*genrand_real2());
		swap_two_integers(&shuffle[i], &shuffle[j]);
	}
	return shuffle;
}

double RK4(double x, double t, double h, double *&values, double (*f)(double, double, double*))
{
	double h2 = 0.5*h;
	double k1 = f(x, t, values);
	double k2 = f(x+h2*k1, t+h2, values);
	double k3 = f(x+h2*k2, t+h2, values);
	double k4 = f(x+h*k3, t+h, values);

	return x + h/6.*(k1+2.*k2+2.*k3+k4);
}

void Adaptive_RK4_h(double *&x, int x_size, double t, double& h, double **&values, double (*f)(double, double, double*))
{
	int i;
	double *x_next1 = (double *) calloc(x_size, sizeof(double));
	double *x_next2 = (double *) calloc(x_size, sizeof(double));
	double x_del;
	double rho = 0.;
	double accuracy = 1e-12; // target accuracy
	// If the difference between two results is less than 'precision', these are same.
	double precision = 1e-6; // 소수 6번째 자리까지 정도 확인

	while(rho<1.-precision){
		x_del = 0.; 
		for(i=0; i<x_size; i++){
			x_next1[i] = RK4(RK4(x[i], 0, h, values[i], f), 0, h, values[i], f);
			x_next2[i] = RK4(x[i], 0, 2.*h, values[i], f);
			x_del += fabs(x_next1[i]-x_next2[i]);
		}
		if(x_del < precision){
			h *= 2.;
			rho = 10.;
		}
		else{
			rho = 30.*h*accuracy/x_del;

			if(rho<1.-precision) h *= pow(rho, 0.25);
			else h *= 2.;
		}
	}
	for(i=0; i<x_size; i++) x[i] = x_next1[i]+(x_next1[i]-x_next2[i])/15.;
	free(x_next1); free(x_next2);
}

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
	new_node->abundance = paras.y0;
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
