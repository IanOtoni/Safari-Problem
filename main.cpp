#include <iostream>
#include <algorithm>
#include <cmath>
#include <random>
#include <numeric>
#include <limits>
#include <string>
#include <vector>
#include <unordered_map>
#include "geometria.h"

using namespace std;
using Polygon = vector<Point>;

void show_data(const string &name, const Point &departure, const Point &arrival, const Polygon &area, const vector<Polygon> &cages){
    cout << "==========================================================================" << endl;
    cout << name << endl;
    cout << "departure: "; departure.show(); cout << endl;
    cout << "arrival: "; arrival.show(); cout << endl;
    cout << "Area: ";
    for(int i=0; i<area.size(); i++){
        area[i].show();
        if(i!=area.size()-1)
            cout << ", ";
    }
    cout << endl;
    cout << "Cages:" << endl;;
    for(int i=0; i<cages.size(); i++){
        cout << "Cage " << i << ": ";
        for(int j=0; j<cages[i].size(); j++){
            cages[i][j].show();
            if(j!=cages[i].size()-1)
                cout << ", ";
        }
        cout << endl;
    }
    cout << "==========================================================================" << endl;
}

void show_tour(vector<int> tour, vector<Point> nodes){
    for(int i : tour){
        nodes[i].show();
        cout << " ";
    }
    cout << endl;
}

int main(){
    //-----------------------------Ant colony parameters-----------------------------//
    const int NUM_ITERATIONS = 1000;
    const int NUM_ANTS       = 40;
    const double ALFA        = 1.0; //pheromone influence
    const double BETA        = 1.0; //heuristic influence
    const double RHO         = 0.2; //evaporation rate
    const double Q           = 100.0; //pheromone deposit factor
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
    vector<vector<int>> cage_vertices(num_cages); //mapping each cage to his vertices

    for(int i=0; i<num_cages; i++){ //for each cage
        for(Point &vertex : cages[i]){ //for each vertex's cage
            nodes.push_back(vertex); //add in nodes
            node_to_cage_map.push_back(i); //maps the vertex in the box i
            cage_vertices[i].push_back(nodes.size()-1);
        }
    }
    node_to_cage_map.push_back(-1); //the last node doest have a cage

    nodes.push_back(arrival);
    int arrival_index = nodes.size() - 1;
    int num_nodes     = nodes.size();

    vector<vector<double>> dist(num_nodes, vector<double>(num_nodes)); //matrix of distances between any pair of points
    for(int i=0; i<num_nodes; i++){ 
        for(int j=0; j<num_nodes; j++){ 
            dist[i][j] = distance(nodes[i], nodes[j]);
        }
    }

    vector<vector<double>> pheromones(num_nodes, vector<double>(num_nodes, 1.0)); //pheromone matrix

    mt19937 gen(random_device{}());
    uniform_real_distribution<> dis(0.0, 1.0);

    double best_length = numeric_limits<double>::infinity();
    vector<int> best_tour; //stores the tour with minimum cost
    //-------------------------------------------------------------------------------//

    //-----------------------------ACO in his pure state-----------------------------//
    for(int iter=0; iter<NUM_ITERATIONS; iter++){
        vector<vector<int>> all_tours(NUM_ANTS);
        vector<double> all_lengths(NUM_ANTS);

        for(int ant=0; ant<NUM_ANTS; ant++){
            vector<int> tour; //each ant has a tour
            tour.push_back(0); //each tour starts with the departure
            vector<bool> visited_cages(num_cages, false); //controls visited cages
            vector<bool> visited_nodes(num_nodes, false); //controls visited cages
            int current_node        = 0;
            int visited_cages_count = 0;

            while(visited_cages_count < num_cages){
                vector<int> candidates;
                vector<double> probabilities;
                double total_prob = 0.0;

                //identify all candidate nodes (vertices of unvisited cages)
                for(int next_node=1; next_node<arrival_index; ++next_node){
                    int cage_of_node = node_to_cage_map[next_node];
                    //if(visited_cages[cage_of_node]) continue;
                    if(visited_nodes[next_node]) continue;
            
                    //-----------------------------------------------------------------------------------//
                    double num_cages_visited = 0;
                    for(int i=0; i<cages.size(); i++){
                        if(visited_cages[i]) continue;

                        if(is_intersection_segment_polygon(nodes[current_node], nodes[next_node], cages[i])){
                            num_cages_visited++;
                        }
                    }
                    //-----------------------------------------------------------------------------------//

                    double distance_val = dist[current_node][next_node];
                    double pheromone    = pow(pheromones[current_node][next_node], ALFA);
                    double heuristic    = pow((1e-9  + num_cages_visited * 1.0)/ (distance_val + 1e-9), BETA); // 1e-9 to avoid division by 0
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
                if(!visited_cages[chosen_cage]){
                    visited_cages[chosen_cage] = true;
                    visited_cages_count++;
                }

                //check if the new segment intercepts some cage
                
                for(int i=0; i<cages.size(); i++){
                    if(!visited_cages[i] && is_intersection_segment_polygon(nodes[last_node], nodes[current_node], cages[i])){
                        visited_cages[i] = true;
                        visited_cages_count++;
                        //cout << "entrou aqui";
                    }
                }
                
            }
            tour.push_back(arrival_index); //ends the tour

            //compute length
            double L = 0;
            for(int i=0; i<tour.size()-1; i++)
                L += dist[tour[i]][tour[i+1]];

            all_tours[ant]   = tour;
            all_lengths[ant] = L;

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
            double delta      = Q / all_lengths[k];
            vector<int> &tour = all_tours[k];
            for(int i = 0; i<tour.size()-1; i++){
                int u = tour[i], v = tour[i+1];
                pheromones[u][v] += delta;
                pheromones[v][u] += delta;
            }
        }
        if(!((iter+1) % 100)) cout << "Iter " << iter+1 << ", best = " << best_length << "\n";
    }
    //-------------------------------------------------------------------------------//

    //-----------------------Print the final result and route-----------------------//
    cout << endl << "=== Final results ===" << endl;
    cout << "Approximate optimal distance: "<< best_length <<"\n";
    cout << "Route: [";
    for(int i=0; i<best_tour.size(); i++){
        nodes[best_tour[i]].show();
        if(i!=best_tour.size()-1) cout << ", ";
    }
    cout << "]" << endl;

    return 0;
    //-------------------------------------------------------------------------------//
}