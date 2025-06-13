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
}

int main() {
    //-----------------------------Ant colony parameters-----------------------------//
    const int NUM_ITERATIONS = 1000;
    const int NUM_ANTS       = 40;
    const double ALFA        = 1.0; //pheromone influence
    const double BETA        = 5.0; //heuristic influence
    const double RHO         = 0.2; //evaporation rate
    const double Q           = 100.0; //pheromone deposit factor
    //-------------------------------------------------------------------------------//

    //----------------Reading the instances and storing the polygons----------------//
    string name;
    double x, y;
    int num_vertices_area, num_cages;
    Polygon area;
    vector<Polygon> cages; //vector<vector<Point>> cages; cages[i] eh vector<Point>

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
    for(int i = 0; i < num_cages; i++) {
        int num_vertices_cage;
        cin >> num_vertices_cage;
        cages[i].resize(num_vertices_cage);
        for(int j = 0; j < num_vertices_cage; j++){
            cin >> x >> y;
            cages[i][j] = Point(x, y);
        }
    }
    //-------------------------------------------------------------------------------//

    //--------------------------Initialize data structures--------------------------//
    vector<Point> nodes; //node list: 0=departure, then each cage's vertices, then arrival
    nodes.push_back(departure);
    vector<int> start_index_of_cages(num_cages); //stores the starting index(in nodes vec) of each box i in the i position

    for(int i=0; i<num_cages; i++) { //for each cage
        start_index_of_cages[i] = nodes.size(); //the starting index will be the sum of all nodes before(nodes.size())
        for(Point &vertex : cages[i]) nodes.push_back(vertex); //add all the vertices of de i-cage to nodes
    }
    
    nodes.push_back(arrival);
    int arrival_index = nodes.size() - 1;
    int num_nodes     = nodes.size();


    vector<vector<double>> dist(num_nodes, vector<double>(num_nodes)); //matrix of distances between any pair of points
    for(int i=0; i<num_nodes; i++) 
        for(int j=0; j<num_nodes; j++) 
            dist[i][j] = distance(nodes[i], nodes[j]);

    vector<vector<double>> pheromones(num_nodes, vector<double>(num_nodes, 1.0)); //pheromone matrix

    mt19937 gen(random_device{}());
    uniform_real_distribution<> dis(0.0, 1.0);

    double best_length = numeric_limits<double>::infinity();
    vector<int> best_tour; //stores the tour with minimum cost
    //-------------------------------------------------------------------------------//

    //-----------------------------ACO in his pure state-----------------------------//
    for(int iter=0; iter<NUM_ITERATIONS; iter++){
        vector<vector<int>> all_tours;
        vector<double> all_lengths;

        for(int ant=0; ant<NUM_ANTS; ant++){
            vector<int> tour; //each ant has a tour
            tour.push_back(0); //each tour starts with the departure
            int current = 0;   //ant actual position
            unordered_map<int, bool> visited_cages; 
       
            for(int c=0; c<num_cages; c++){ //visit each cage in order
                // [start, end)
                int start = start_index_of_cages[c]; //first vertex of the cage
                int end   = (c+1 < num_cages) ? start_index_of_cages[c+1] : arrival_index; //last vertex + 1 of the cage

                vector<int> cand; //vec to store candidates for the next step
                vector<double> probs; //vec to stores the weights of the candidates (pheromones·heuristic)
                double sum_prob = 0;

                for(int v=start; v<end; v++){ //for each vertex of the cage
                    double pher = pow(pheromones[current][v], ALFA);
                    double heur = pow(1.0/(dist[current][v]+1e-6), BETA);
                    double p = pher * heur;
                    cand.push_back(v);
                    probs.push_back(p);
                    sum_prob += p;
                }

                //draw the next step
                double r = dis(gen) * sum_prob;
                double acc = 0;
                int chosen = cand.back();
                for(int i=0; i<cand.size(); i++){
                    acc += probs[i];
                    if(r <= acc) {
                        chosen = cand[i];
                        break; 
                    }
                }
                tour.push_back(chosen);
                current = chosen;
            }

            //ends the tour
            tour.push_back(arrival_index);

            //compute length
            double L = 0;
            for(int i=0; i<tour.size()-1; i++)
                L += dist[tour[i]][tour[i+1]];
            all_tours.push_back(tour);
            all_lengths.push_back(L);
            if(L < best_length){
                best_length = L; 
                best_tour = tour;
            }
            
        }

        for(int i=0; i<num_nodes; i++) //pheromones evaporation 
            for(int j=0; j<num_nodes; j++)
                pheromones[i][j] *= (1.0 - RHO);

        //update pheromones
        for(int k=0; k<all_tours.size(); k++){
            double delta = Q / all_lengths[k];
            auto &tour = all_tours[k];
            for(int i = 0; i+1 < tour.size(); i++){
                int u = tour[i], v = tour[i+1];
                pheromones[u][v] += delta;
                pheromones[v][u] += delta;
            }
        }
        cout << "Iter "<< iter+1 << ", best = " << best_length << "\n";
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