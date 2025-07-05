#include <iostream>
#include <algorithm>
#include <cmath>
#include <random>
#include <numeric>
#include <limits>
#include <string>
#include <vector>
#include <bitset>
#include <unordered_map>
#include "geometria.h"

using namespace std;
using Polygon = vector<Point>;

#define MAX_CAGES 64

void show_data(const string &name, const Point &departure, const Point &arrival, const Polygon &area, const vector<Polygon> &cages){
    cerr << "==========================================================================" << endl;
    cerr << name << endl;
    cerr << "departure: "; departure.show(); cerr << endl;
    cerr << "arrival: "; arrival.show(); cerr << endl;
    cerr << "Area: ";
    for(int i=0; i<area.size(); i++){
        area[i].show();
        if(i!=area.size()-1)
            cerr << ", ";
    }
    cerr << endl;
    cerr << "Cages:" << endl;;
    for(int i=0; i<cages.size(); i++){
        cerr << "Cage " << i << ": ";
        for(int j=0; j<cages[i].size(); j++){
            cages[i][j].show();
            if(j!=cages[i].size()-1)
                cerr << ", ";
        }
        cerr << endl;
    }
    cerr << "==========================================================================" << endl;
}

bitset<MAX_CAGES> create_mask(const vector<int> &numbers){
    bitset<MAX_CAGES> mask; //creates a bitset of MAX_CAGES bits, all initialized to 0
    
    for(int num : numbers){
        mask.set(num);
    }
    return mask;
}

/*
Important variables:
verificar codigo de interseccao segmento-poligono para poligonos nao convexos
num_vertices_area : Number of vertices of the big area
num_cages: number of cages
cages: vector that stores in each position a cage (vector that stores in each position a pair (x, y) that represents a cage vertex)
departure: pair (x, y) that represents a departure point =)
arrival: pair (x, y) that represents a arrival point point =)
nodes: vector of pairs with all the points, the position 0 is the departure and the last position is the arrival, between then are all the points
node_to_cage_map: maps a node to his cage (departure e arrival belongs to cage -1)
distances: matrix that stores in the position i, j the size (double) of the segment ij
cages_intersections: matrix that stores in the position i, j a bitset that is 1 in the position k if the k cage is intersected by the segment i, j and 0 otherwise
ant_cages_intersections: matrix that stores in the position i, j a bitset that is 1 in the position k if the ant already visited the k cage
*/

int main(int argc, char* argv[]){
    ///-----------------------------Ant colony parameters-----------------------------//
    //Default values ​​if no parameters are passed
    int NUM_ITERATIONS = 1000;
    int NUM_ANTS       = 53;
    double ALFA        = 5.1877; //pheromone influence
    double BETA        = 3.8149; //heuristic influence
    double RHO         = 0.2053; //evaporation rate
    double Q           = 179.5837; //pheromone deposit factor

    //Loop to read command line parameters
    for(int i=1; i<argc; i+=2){
        string param = argv[i];

        if(param == "--iter"){
            NUM_ITERATIONS = stoi(argv[i + 1]);

        }else if (param == "--ants"){
            NUM_ANTS = stoi(argv[i + 1]);

        }else if (param == "--alfa"){
            ALFA = stod(argv[i + 1]);

        }else if (param == "--beta"){
            BETA = stod(argv[i + 1]);

        }else if (param == "--rho"){
            RHO = stod(argv[i + 1]);

        }else if (param == "--q"){
            Q = stod(argv[i + 1]);
        }
    }
    //-------------------------------------------------------------------------------//

    //----------------Reading the instances and storing the polygons----------------//
    string name;
    double x, y;
    int num_vertices_area, num_cages;
    Polygon area;
    vector<Polygon> cages; //vector<vector<Point>> cages; cages[i] is a vector<Point>

    getline(cin, name);
    cin >> x >> y;
    Point departure(x, y);
    cin >> x >> y;
    Point arrival(x, y);
    cin >> num_vertices_area;
    area.resize(num_vertices_area);
    for(int i = 0; i < num_vertices_area; i++){
        cin >> x >> y;
        area[i] = Point(x, y);
    }
   
    cin >> num_cages;
    cages.resize(num_cages);
    for(int i = 0; i < num_cages; i++){
        int num_vertices_cage;
        cin >> num_vertices_cage;
        cages[i].resize(num_vertices_cage);
        for(int j = 0; j < num_vertices_cage; j++){
            cin >> x >> y;
            cages[i][j] = Point(x, y);
        }
    }

    show_data(name, departure, arrival, area, cages);
    //-------------------------------------------------------------------------------//

    //--------------------------Initialize data structures--------------------------//
    vector<Point> nodes; //node list: 0=departure, then each cage's vertices, then arrival
    nodes.push_back(departure);
    vector<int> node_to_cage_map;
    node_to_cage_map.push_back(-1); //mapping each vertex to the index of his cage

    for(int i=0; i<num_cages; i++){ //for each cage
        for(Point &vertex : cages[i]){ //for each vertex's cage
            nodes.push_back(vertex); //add in nodes
            node_to_cage_map.push_back(i); //maps the vertex in the box i
        }
    }
    node_to_cage_map.push_back(-1); //the last node doest have a cage

    nodes.push_back(arrival);
    int arrival_index = nodes.size() - 1;
    int num_nodes     = nodes.size();

    vector<vector<double>> distances(num_nodes, vector<double>(num_nodes)); //matrix of distances between any pair of points
    for(int i=0; i<num_nodes; i++){ 
        for(int j=0; j<num_nodes; j++){ 
            distances[i][j] = distance(nodes[i], nodes[j]);
        }
    }

    vector<vector<bitset<MAX_CAGES>>> cages_intersections(num_nodes, vector<bitset<MAX_CAGES>>(num_nodes)); //matrix that stores intersected cages(non visited) by segment
    for(int i=0; i<num_nodes; i++){ 
        for(int j=0; j<num_nodes; j++){
            vector<int> intersections;
            for(int k=0; k<num_cages; k++){
                if(is_intersection_segment_polygon(nodes[i], nodes[j], cages[k])){
                    intersections.push_back(k);
                }
            }
            cages_intersections[i][j] = create_mask(intersections);
        }
    }

    vector<vector<double>> pheromones(num_nodes, vector<double>(num_nodes, 1.0)); //pheromone matrix

    mt19937 gen(random_device{}());
    uniform_real_distribution<> dis(0.0, 1.0);

    double best_length = numeric_limits<double>::infinity();
    vector<int> best_tour; //stores the tour with minimum cost
    //-------------------------------------------------------------------------------//

    //-----------------------------ACO in its pure state-----------------------------//
    for(int iter=0; iter<NUM_ITERATIONS; iter++){
        vector<vector<int>> all_tours(NUM_ANTS);
        vector<double> tour_lengths(NUM_ANTS);

        for(int ant=0; ant<NUM_ANTS; ant++){
            vector<int> tour; //each ant has a tour
            tour.push_back(0); //each tour starts with the departure
            bitset<MAX_CAGES> ant_visited_cages; //controls visited cages
            vector<bool> visited_nodes(num_nodes, false); //controls visited cages
            int current_node = 0;

            while(ant_visited_cages.count() < num_cages){
                vector<int> candidates;
                vector<double> probabilities;
                double total_prob = 0.0;

                //identify all candidate nodes (unvisited vertices)
                for(int next_node=1; next_node<arrival_index; next_node++){
                    if(visited_nodes[next_node]) continue;
            
                    int num_new_cages   = ((ant_visited_cages ^ cages_intersections[current_node][next_node]) & cages_intersections[current_node][next_node]).count();
                    double distance_val = distances[current_node][next_node];
                    double pheromone    = pow(pheromones[current_node][next_node], ALFA);
                    double heuristic    = pow(0.1*(double)num_new_cages, BETA); //pow((double)num_new_cages / distance_val, BETA);
                    double prob         = pheromone * heuristic;

                    candidates.push_back(next_node);
                    probabilities.push_back(prob);
                    total_prob += prob;
                }

                //roulette wheel selection
                double r = dis(gen) * total_prob;
                double prob_sum = 0.0;
                int chosen_node = candidates.back();

                for(int i=0; i<candidates.size(); i++){
                    prob_sum += probabilities[i];
                    if (r <= prob_sum){
                        chosen_node = candidates[i]; //next node
                        break;
                    }
                }

                //update the ant's state
                tour.push_back(chosen_node);
                int chosen_cage = node_to_cage_map[chosen_node];
                int last_node   = current_node;
                current_node    = chosen_node;
                visited_nodes[chosen_node] = true;
                ant_visited_cages.set(chosen_cage);

                bitset<MAX_CAGES> new_cages_visited = (ant_visited_cages ^ cages_intersections[last_node][current_node]) & cages_intersections[last_node][current_node];
                ant_visited_cages |= new_cages_visited; //or
            }
            tour.push_back(arrival_index); //ends the tour

            //compute length
            double L = 0;
            for(int i=0; i<tour.size()-1; i++)
                L += distances[tour[i]][tour[i+1]];

            all_tours[ant]    = tour;
            tour_lengths[ant] = L;

            if(L<best_length){
                best_length = L; 
                best_tour   = tour;
            }
            
        }

        //pheromones evaporation
        for(int i=0; i<num_nodes; i++) 
            for(int j=0; j<num_nodes; j++)
                pheromones[i][j] *= (1.0 - RHO);

        //update pheromones
        for(int k=0; k<NUM_ANTS; k++){
            double delta      = Q / tour_lengths[k];
            for(int i = 0; i<all_tours[k].size()-1; i++){
                int u = all_tours[k][i];
                int v = all_tours[k][i+1];
                pheromones[u][v] += delta;
                pheromones[v][u] += delta;
            }
        }
        if(!((iter+1) % 100)) cerr << "Iter " << iter+1 << ", best = " << best_length << "\n";
    }
    //-------------------------------------------------------------------------------//

    //-----------------------Print the final result and route-----------------------//
    cerr << endl << "====== Final results ======" << endl;
    cerr << "Approximate optimal distance: "<< best_length <<"\n";
    cerr << "Route: [";
    for(int i=0; i<best_tour.size(); i++){
        nodes[best_tour[i]].show();
        if(i!=best_tour.size()-1) cerr << ", ";
    }
    cerr << "]" << endl;

    cout << best_length << endl;

    return 0;
    //-------------------------------------------------------------------------------//
}