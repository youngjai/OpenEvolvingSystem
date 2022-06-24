/*
	This code is a collection of fucntions
	for calculating the standard score of motifs.

	It was written by Youngjai in Feb. 18, 2021.

	input variables
	date(yymmdd), sigma, mutant event time, saturation time, dt(K), alpha (m=5, alpha=0.1)

	* x_i/K -> y_i , y_th = 1/K, y_0 = 10*y_th (To remove the carrying capacity)
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
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
	int K; // the carring capacity
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

// from youngjai_packages.c
clock_t time_begin(char *path);
void time_end(clock_t begin, char *path);
int *shuffle_array_int(int length, unsigned long seed);
extern int call; // double normal_dist(double avg, double stdev, unsigned long seed)

// from networkYJ.c
void initial_SPP_network(Network& network, Parameters& paras);
void free_node_list(Node_list *&node_list);
void create_node(Network& network, Parameters& paras);
void create_connection(Network& network, Parameters& paras, int name_i, int name_j);

double *RK4(double *f, double t, double h, Network network, double *(*dfdt)(double *, double, Network));
double *Adaptive_RK4(double *f, double t, double &h, Network network, double *(*dfdt)(double *, double, Network));

// In this c file
void simul_SPP_model(Parameters& paras, double **&results, double **&network_edges, int *&num_of_S);
void competitive_dynamics(Network& network, Parameters& paras, double& h, double t_i, double t_f);
void appear_new_node(Network& network, Parameters& paras);
void remove_a_node(double t, Network& network, Parameters& paras);
double *dfdt(double *f, double t, Network network);
void save_parameters(Parameters& paras, char *path);
void print_results(Parameters& paras, double **&results, double **&network_edges, \
	int len_rows, int len_columns, int len_links, char *path);
void print_num(int *&num_of_S, int len_columns, char *path);

// temporary functions
void print_network_test(Network& network);

int main(int argc, char **argv)
{
	if(argc != 7) {
		printf("Error! wrong number of arguments.\r\n");
		return 1;
	}
	// make folder which to save results
	char *str = (char *) calloc(100, sizeof(char)); strcpy(str, argv[1]);
	mkdir(str, 0755);
	// To write simulation parameters and running time
	char *log_path = (char *) calloc(100, sizeof(char)); 
	sprintf(log_path, "%s/log.dat", str);
	clock_t runtime = time_begin(log_path);

	// Input the initial parameters
	Parameters paras;
	paras.N0 = 100; paras.y_th = 1./atof(argv[5]); paras.y0 = 10.*paras.y_th;
	paras.sigma = pow(10, atof(argv[2])*0.5); paras.m = 5; 
	paras.alpha = pow(10, atof(argv[6]));
	paras.mutant_event_time = atoi(argv[3]); paras.saturation_time = atoi(argv[4]);
	paras.simul_step = paras.mutant_event_time+paras.saturation_time; 
	// d\tau (simul_dt) = K dt (K: carrying capacity, tau: reduced time)
	paras.simul_dt = atof(argv[5]);
	// Random seed is loaded depending current time when the code is started.
	paras.seed = time(NULL);

	int i, j;
	// the total number of species: N0 + mutant evet time
	int len_rows = paras.N0+paras.mutant_event_time;
	// [0]: name, [1]: initial abundances
	// [2-end]: abundance time series (mutant_event_time + saturation time)
	int len_columns = 1+1+paras.mutant_event_time+paras.saturation_time;
	// the total links: initial links(0.5*N0*m) + invasion links(T*m) + selfloop(N0+T)
	int len_links = (0.5*paras.N0+paras.mutant_event_time)*paras.m+(paras.N0+paras.mutant_event_time);

	int *num_of_S = (int *) calloc(len_columns-2, sizeof(int));
	double **results = (double **) calloc(len_rows, sizeof(double*));
	double **network_edges = (double **) calloc(len_links, sizeof(double*));
	for(i=0; i<len_rows; i++){
		// [0]: name, 
		// [1-end]: abundance time series (mutant_event_time + saturation time)
		results[i] = (double *) calloc(len_columns, sizeof(double));
		results[i][0] = i+1; 
		// initialization of the 'results' array
		for(j=1; j<len_columns; j++) results[i][j] = 0.;
	}
	// To save the network structure
	for(i=0; i<len_links; i++){
		network_edges[i] = (double *) calloc(3, sizeof(double));
		network_edges[i][0] = -1.;
	}

	save_parameters(paras, log_path);
	simul_SPP_model(paras, results, network_edges, num_of_S);
	print_results(paras, results, network_edges, \
		len_rows, len_columns, len_links, str);
	print_num(num_of_S, len_columns, str);

	for(i=0; i<len_rows; i++) free(results[i]);
	for(i=0; i<len_links; i++) free(network_edges[i]);
	free(num_of_S);
	free(results);
	free(network_edges);
	time_end(runtime, log_path);
	free(str); free(log_path);
	return 0;
}

void simul_SPP_model(Parameters& paras, double **&results, double **&network_edges, int *&num_of_S)
{
	int i, j, t;
	double h = paras.simul_dt/paras.alpha; // integration step size
	double simul_t = 0.;
	Network network;
	Node *cursor;

	call = 0;
	j = 0;
	// abundance time series region is from 1 to end
	initial_SPP_network(network, paras);
	for(i=0; i<paras.N0; i++) results[i][1] = paras.y0;
	cursor = network.node_list->head;
	while(cursor->next != NULL){
		cursor = cursor->next;
		network_edges[j][0] = cursor->node_name;
		network_edges[j][1] = cursor->node_name;
		network_edges[j++][2] = cursor->self_loop;
		for(i=0; i<cursor->degree; i++) if(cursor->direction[i] > 0){
			network_edges[j][0] = cursor->node_name;
			network_edges[j][1] = cursor->neighbors[i];
			network_edges[j++][2] = cursor->sign[i];
		}
	} cursor = NULL;

	// a new species appear every 1/alpha steps
	for(t=0; t<paras.mutant_event_time; t++){
		num_of_S[t] = network.nodes;
		paras.seed++;
		while(simul_t<(t+1)*paras.simul_dt/paras.alpha){
			competitive_dynamics(network, paras, h, simul_t, (t+1)*paras.simul_dt/paras.alpha);
			remove_a_node(simul_t, network, paras);
			simul_t += h;
		}
		// record the abundance dist.
		cursor = network.node_list->head;
		while(cursor->next != NULL){
			cursor = cursor->next;
			results[cursor->node_name-1][t+2] = cursor->abundance;
		} cursor = NULL;

		appear_new_node(network, paras);

		// record the edge list of a new node
		cursor = network.node_list->tail;
		network_edges[j][0] = cursor->node_name;
		network_edges[j][1] = cursor->node_name;
		network_edges[j++][2] = cursor->self_loop;
		for(i=0; i<cursor->degree; i++){
			if(cursor->direction[i] > 0){
				network_edges[j][0] = cursor->node_name;
				network_edges[j][1] = cursor->neighbors[i];
			}
			else{
				network_edges[j][0] = cursor->neighbors[i];
				network_edges[j][1] = cursor->node_name;						
			}
			network_edges[j++][2] = cursor->sign[i];
		} cursor = NULL;
	}

	// saturation process
	for(t=paras.mutant_event_time; t<paras.simul_step; t++){
		num_of_S[t] = network.nodes;
		paras.seed++;
		while(simul_t<(t+1)*paras.simul_dt/paras.alpha){
			competitive_dynamics(network, paras, h, simul_t, (t+1)*paras.simul_dt/paras.alpha);
			remove_a_node(simul_t, network, paras);
			simul_t += h;
		}
		// record the abundance dist.
		cursor = network.node_list->head;
		while(cursor->next != NULL){
			cursor = cursor->next;
			results[cursor->node_name-1][t+2] = cursor->abundance;
		} cursor = NULL;
	}
	// print_network_test(network);
	free_node_list(network.node_list);
	free(cursor);
}

void competitive_dynamics(Network& network, Parameters& paras, double& h, double t_i, double t_f)
{
	int i;

	double *f = (double *) calloc(network.nodes, sizeof(double));
	double *f_tmp;
	Node *cursor = network.node_list->head; i = 0; 
	while(cursor->next != NULL){
		cursor = cursor->next;
		f[i] = cursor->abundance; i++;
	} cursor = NULL;

	// Adaptive RK4
	f_tmp = Adaptive_RK4(f, t_i, h, network, dfdt);
	// If h is larger than t_f-t_i, it just recalculate abundances with RK4
	if(h > t_f-t_i){
		h = t_f-t_i;
		f_tmp = RK4(f, t_i, h, network, dfdt);
	}
	cursor = network.node_list->head; i = 0;
	while(cursor->next != NULL){
		cursor = cursor->next;
		cursor->abundance = f_tmp[i]; i++;
	} cursor = NULL;

	// // RK4
	// if(h > t_f-t_i) h = t_f-t_i;
	// else h = 1./paras.alpha;
	// f_tmp = RK4(f, t_i, h, network, dfdt);
	// cursor = network.node_list->head; i = 0;
	// while(cursor->next != NULL){
	// 	cursor = cursor->next;
	// 	cursor->abundance = f_tmp[i]; i++;
	// } cursor = NULL; 

	free(f); free(f_tmp); free(cursor);
}

void appear_new_node(Network& network, Parameters& paras)
{
	int i;
	int count = 0;
	int m_temp, name_i;
	int *shuffled = shuffle_array_int(network.nodes*2, paras.seed);
	Node *cursor;

	/* incubation rule. */
	if(network.nodes <= paras.m) m_temp = network.nodes;
	else m_temp = paras.m;

	create_node(network, paras);

	while(count < m_temp){
		cursor = network.node_list->head->next;
		for(i=0; i<shuffled[count]/2; i++) cursor = cursor->next;
		name_i = cursor->node_name;

		create_connection(network, paras, network.node_label-1, name_i);
		count++;
	}
	free(shuffled);
	cursor = NULL; free(cursor);
}

void remove_a_node(double t, Network& network, Parameters& paras)
{
	int i, j, k;
	Node *cursor = network.node_list->head;
	Node *cursor1;

	while(cursor->next != NULL){
		if(cursor->next->abundance < paras.y_th){
			for(i=0; i<cursor->next->degree; i++){
				cursor1 = network.node_list->head->next;
				while(cursor1->node_name != cursor->next->neighbors[i])
					cursor1 = cursor1->next;
				j = 0;
				while(cursor1->neighbors[j++] != cursor->next->node_name);
				cursor1->degree--; j--;
				cursor1->neighbors[j] = cursor1->neighbors[cursor1->degree];
				cursor1->direction[j] = cursor1->direction[cursor1->degree];
				cursor1->sign[j] = cursor1->sign[cursor1->degree];
			}
			network.nodes--;
			network.links -= cursor->next->degree;
			cursor1 = cursor->next;
			cursor->next = cursor1->next;
			if(cursor1->node_name == network.node_list->tail->node_name){
				network.node_list->tail = cursor;
			}

			free(cursor1->neighbors);
			free(cursor1->direction);
			free(cursor1->sign);
			free(cursor1);
		}
		if(network.node_list->tail != cursor) cursor = cursor->next;
	}
	cursor = NULL; free(cursor);
	cursor1 = NULL; free(cursor1);
}

double *dfdt(double *f, double t, Network network)
{
	int i, j, k;
		
	// the summation of all abundances
	double f_sum = 0.;
	for(i=0; i<network.nodes; i++) f_sum += f[i];

	// the calculation of growth rate G_i and death rate D_i
	// The calloc() function in C is used to allocate 
	// a specified amount of memory and then initialize it to zero.
	double *G = (double *) calloc(network.nodes, sizeof(double));
	double *D = (double *) calloc(network.nodes, sizeof(double));
	int *index = (int *) calloc(network.nodes, sizeof(int)); 
	Node *cursor = network.node_list->head; i = 0;
	while(cursor->next != NULL){
		cursor = cursor->next;
		index[i] = cursor->node_name;
		G[i] = cursor->self_loop * f[i]; i++;
	}
	cursor = network.node_list->head; i = 0;
	while(cursor->next != NULL){
		cursor = cursor->next;
		for(j=0; j<cursor->degree; j++) if(cursor->direction[j] < 0){
			k = -1; while(index[++k] != cursor->neighbors[j]);
			if(cursor->sign[j] > 0) G[i] += cursor->sign[j] * f[k];
			else D[i] += cursor->sign[j] * f[k];
		} i++;
	} cursor = NULL;

	// df/dt array
	double *df = (double *) calloc(network.nodes, sizeof(double));
	for(i=0; i<network.nodes; i++) df[i] = f[i]*(G[i]*(1-f_sum)+D[i]);

	free(G); free(D); free(index); free(cursor);
	return df;
}

void save_parameters(Parameters& paras, char *path)
{
	FILE *file = fopen(path, "at");
	fprintf(file, "N0 : %d\n", paras.N0);
	fprintf(file, "y_th : %.2le\n", paras.y_th);
	fprintf(file, "y0 : %.2le\n", paras.y0);
	fprintf(file, "m : %d\n", paras.m);
	fprintf(file, "sigma : %.2le\n", paras.sigma);
	fprintf(file, "alpha : %.2le\n", paras.alpha);
	fprintf(file, "mutant event time : %d\n", paras.mutant_event_time);
	fprintf(file, "saturation time : %d\n", paras.saturation_time);
	fprintf(file, "simulation steps : %d\n", paras.simul_step);
	fprintf(file, "simulation unit time : %lf\n", paras.simul_dt);
	fprintf(file, "seed : %lu\n", paras.seed);
	fclose(file);
}

void print_results(Parameters& paras, double **&results, double **&network_edges, \
	int len_rows, int len_columns, int len_links, char *path)
{
	int i, j;
	char *str = (char *) calloc(100, sizeof(char));
	FILE *file1, *file2;

	sprintf(str, "%s/simul_log.txt", path); file1 = fopen(str, "wt");
	fprintf(file1, "name,abundance_series\n");
	for(i=0; i<len_rows; i++){
		fprintf(file1, "%d,", (int) results[i][0]);
		for(j=1; j<len_columns-1; j++) fprintf(file1, "%le,", results[i][j]);
		fprintf(file1, "%le\n", results[i][len_columns-1]);
	} fclose(file1);

	sprintf(str, "%s/network_edges.txt", path); file2 = fopen(str, "wt");
	fprintf(file2, "Source,Target,Label,Weight\n");
	for(i=0; i<len_links; i++){
		if(network_edges[i][0]>0){
			fprintf(file2, "%d,%d,", (int) network_edges[i][0], (int) network_edges[i][1]);
			if(network_edges[i][2]>0) fprintf(file2, "1,%.4le\n", network_edges[i][2]);
			else fprintf(file2, "-1,%.4le\n", network_edges[i][2]);
		}
	}fclose(file2); 
	free(str);
}

void print_num(int *&num_of_S, int len_columns, char *path)
{
	int i;
	char *str = (char *) calloc(100, sizeof(char));
	sprintf(str, "%s/num_of_S.txt", path); FILE *file = fopen(str, "wt");

	for(i=0; i<len_columns-3; i++) fprintf(file, "%d,", num_of_S[i]);
	fprintf(file, "%d\n", num_of_S[len_columns-3]);
	fclose(file); free(str);
}


void print_network_test(Network& network)
{
	int i;
	Node *cursor = network.node_list->head;
	printf("\n");
	while(cursor->next != NULL){
		cursor = cursor->next;
		printf("%d %d %.4lf :\t", cursor->node_name, cursor->degree, cursor->abundance);
		for(i=0; i<cursor->degree; i++)
			printf("%d[%d, %.2le]\t", cursor->neighbors[i], cursor->direction[i], cursor->sign[i]);
		printf("\n");
	}
	printf("nodes:%d links:%d\n", network.nodes, network.links);
	cursor = NULL; free(cursor);
}
