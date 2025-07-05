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
#include <unordered_set>
#include "geometria.h" 

using namespace std;
using Polygon = vector<Point>;

#define MAX_CAGES 100 
vector<vector<double>> distances;

bitset<MAX_CAGES> create_mask(const vector<int> &numbers)
{
    bitset<MAX_CAGES> mask; // cria um bitset de MAX_CAGES bits, todos inicializados com 0.

    for (int num : numbers)
    {
        mask.set(num);
    }
    return mask;
}

bitset<MAX_CAGES> create_mask_solution(const vector<int> &solution, const vector<vector<bitset<MAX_CAGES>>> &cages_intersections)
{
    bitset<MAX_CAGES> mask; // cria um bitset de MAX_CAGES bits, todos inicializados com 0.

    for (int i = 0; i < solution.size() - 1; ++i)
    {
        mask |= cages_intersections[solution[i]][solution[i + 1]];
    }
    return mask;
}

double tour_length(const vector<int> &tour)
{
    double length = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i)
    {
        length += distances[tour[i]][tour[i + 1]];
    }
    return length;
}

bool is_feasible_tour(const vector<int> &tour, int end_index, int num_cages, const vector<vector<bitset<MAX_CAGES>>> &cages_intersections)
{
    bitset<MAX_CAGES> visited_cages_mask; // controla as jaulas visitadas
    for (size_t i = 0; i < tour.size() - 1; ++i)
    {
        int current_node = tour[i];
        int next_node = tour[i + 1];
        visited_cages_mask |= cages_intersections[current_node][next_node];
    }
    return visited_cages_mask.count() == num_cages; // -2 para ignorar o nó de partida e chegada
}

vector<int> repair_solution( vector<int> tour,int end_index,int num_cages,const vector<vector<bitset<MAX_CAGES>>> &cages_intersections,const vector<Point> &nodes, const vector<int> &node_to_cage_map)
{
   
    bitset<MAX_CAGES> visited = create_mask_solution(tour, cages_intersections);

    for (int cage = 0; cage < num_cages; ++cage) {
        if (visited.test(cage)) continue;

        // coleta todos os nós da jaula
        vector<int> candidates;
        for (int i = 0; i < (int)node_to_cage_map.size(); ++i)
            if (node_to_cage_map[i] == cage)
                candidates.push_back(i);

        double best_delta    = numeric_limits<double>::infinity();
        int    best_pos      = -1;
        int    best_node_ins = -1;

        // busca melhor inserção interna
    
        for (int i = 0; i + 1 < (int)tour.size(); ++i) {
            int u = tour[i], v = tour[i+1];
            auto orig = cages_intersections[u][v];
            for (int k : candidates) {
                if (find(tour.begin(), tour.end(), k) != tour.end()) continue;
                auto neu = cages_intersections[u][k] | cages_intersections[k][v];
                if ((orig & ~neu).any()) continue;

                double delta = distances[u][k] + distances[k][v] - distances[u][v];
                if (delta < best_delta) {
                    best_delta    = delta;
                    best_pos      = i + 1;
                    best_node_ins = k;
                }
            }
        }

        if (best_pos >= 0) {
            // inserção interna
            tour.insert(tour.begin() + best_pos, best_node_ins);
        }
        else {
    
            // escolhe o k que dá menor delta entre tour[size-2] e tour[size-1]
            int  u = tour[tour.size()-2];
            int  v = tour.back();
            double min_delta = numeric_limits<double>::infinity();
            int    pick_k    = -1;
            for (int k : candidates) {
                double delta = distances[u][k] + distances[k][v] - distances[u][v];
                if (delta < min_delta) {
                    min_delta  = delta;
                    pick_k     = k;
                }
            }
            // insere antes do final
            tour.insert(tour.end() - 1, pick_k);
        }

        visited = create_mask_solution(tour, cages_intersections);
    }

    return tour;
}

vector<int> hill_climbing_best_improvement( vector<int> initial_tour,int num_cages,int end_index,const vector<vector<bitset<MAX_CAGES>>> &cages_intersections,const vector<Point> &nodes,const vector<int> &node_to_cage_map)
{
    vector<int> current_tour = initial_tour;
    bitset<MAX_CAGES> visited_cages_mask = create_mask_solution(current_tour, cages_intersections);

    bool improvement_found;
    do
    {
        improvement_found = false;
        double current_length = tour_length(current_tour);

        vector<int> best_valid_neighbor_tour = current_tour;
        double best_valid_neighbor_length = current_length;

        // Gera e avalia a vizinhança 2-opt
        for (size_t i = 1; i < current_tour.size() - 2; ++i)
        {
            for (size_t j = i + 1; j < current_tour.size() - 1; ++j)
            {
                vector<int> neighbor_tour = current_tour;
                reverse(neighbor_tour.begin() + i, neighbor_tour.begin() + j + 1);

                // Verifica se ainda cobre todas as jaulas
                bitset<MAX_CAGES> neighbor_mask = create_mask_solution(neighbor_tour, cages_intersections);

                if (neighbor_mask.count() < num_cages)
                {
                    // Repara o tour para cobrir todas as jaulas
                    neighbor_tour = repair_solution(neighbor_tour, end_index, num_cages, cages_intersections, nodes, node_to_cage_map);

                }

                double neighbor_length = tour_length(neighbor_tour);
                if (neighbor_length < best_valid_neighbor_length)
                {
                    best_valid_neighbor_length = neighbor_length;
                    best_valid_neighbor_tour = neighbor_tour;
                }
            }
        }

        if (best_valid_neighbor_length < current_length)
        {
            current_tour = best_valid_neighbor_tour;
            improvement_found = true;
        }

    } while (improvement_found);

    return current_tour;
}


vector<int> generate_initial_solution(int end_index, int num_nodes, int num_cages, const vector<int> &node_to_cage_map, const vector<vector<bitset<MAX_CAGES>>> &cages_intersections, mt19937 &gen)
{
    vector<int> tour;                             
    tour.push_back(0);                            
    bitset<MAX_CAGES> visited_cages_mask;        
    vector<bool> visited_nodes(num_nodes, false);
    int current_node = 0;

    while (visited_cages_mask.count() < num_cages)
    {
        int current_node = tour.back();
        vector<int> candidates;
        for (int next_node = 1; next_node < end_index; ++next_node)
        {
            if (!visited_cages_mask.test(node_to_cage_map[next_node]))
            {
                candidates.push_back(next_node);
            }
        }

        if (candidates.empty())
        {
            cout << "No more candidates available from node " << current_node << ". Ending tour." << endl;
            break;
        }

        uniform_real_distribution<> dis(0.0, 1.0);
        double r = dis(gen);
        int chosen_index = static_cast<int>(r * candidates.size());
        chosen_index = min(chosen_index, static_cast<int>(candidates.size()) - 1); // Garantir que o índice esteja dentro dos limites
        int chosen_node = candidates[chosen_index];
        visited_cages_mask |= cages_intersections[current_node][chosen_node];

        tour.push_back(chosen_node);
    }

    tour.push_back(end_index);

    return tour;
}

// Função principal do Simulated Annealing
vector<int> simulated_annealing(int num_nodes,int end_index,int num_cages,const vector<int> &node_to_cage_map,const vector<vector<bitset<MAX_CAGES>>> &cages_intersections,const vector<Point> &nodes)
{
    // Parâmetros do SA (valores corrigidos do irace)
    // alpha é a taxa de resfriamento, T_min é a temp. mínima
    double T = 5010.1919;
    double alpha = 0.8851;
    double T_min = 0.6306;

    mt19937 gen(random_device{}());
    uniform_real_distribution<> dis_real(0.0, 1.0);

    cout << "Generating initial solution..." << endl;
    vector<int> current_tour = generate_initial_solution(end_index, num_nodes, num_cages, node_to_cage_map, cages_intersections, gen);
    current_tour = repair_solution(current_tour, end_index, num_cages, cages_intersections, nodes, node_to_cage_map);
    
    double current_cost = tour_length(current_tour);
    vector<int> best_tour = current_tour;
    double best_cost = current_cost;

    cout << "Initial tour cost: " << current_cost << endl;

    while (T > T_min)
    {
        
    
        if (current_tour.size() <= 3) {
            break; // Sai do loop se o tour for muito pequeno para o 2-opt
        }

        vector<int> neighbor = current_tour;

        // Usa uniform_int_distribution para gerar índices de forma segura (as partes de geração aleatória foram feitas com auxilio de IA)
        uniform_int_distribution<> dis_i(1, neighbor.size() - 3);
        int i = dis_i(gen);
        uniform_int_distribution<> dis_j(i + 1, neighbor.size() - 2);
        int j = dis_j(gen);
        
        reverse(neighbor.begin() + i, neighbor.begin() + j + 1);
       

        bitset<MAX_CAGES> mask = create_mask_solution(neighbor, cages_intersections);
        if (mask.count() < num_cages)
        {
            neighbor = repair_solution(neighbor, end_index, num_cages, cages_intersections, nodes, node_to_cage_map);
        }

        double neighbor_cost = tour_length(neighbor);
        double delta = neighbor_cost - current_cost;

     
        if (delta < 0 || (T > 0 && dis_real(gen) < exp(-delta / T)))
        {
            current_tour = neighbor;
            current_cost = neighbor_cost;

            if (current_cost < best_cost)
            {
                best_tour = current_tour;
                best_cost = current_cost;
            }
        }
        T *= alpha;
    }
    
    cout << "Best tour cost found: " << best_cost << endl;
    cout << "Improving final solution with hill climbing..." << endl;
    best_tour = hill_climbing_best_improvement(best_tour, num_cages, end_index, cages_intersections, nodes, node_to_cage_map);

    cout << "Final tour cost: " << tour_length(best_tour) << endl;
    return best_tour;
}
int main()
{
    string name;
    double x, y;
    int num_vertices_safari, num_cages;

    getline(cin, name);
    cin >> x >> y;
    Point start(x, y);
    cin >> x >> y;
    Point end(x, y);
    cin >> num_vertices_safari;
    Polygon safari(num_vertices_safari);
    for (int i = 0; i < num_vertices_safari; i++)
    {
        cin >> x >> y;
        safari[i] = Point(x, y);
    }

    cin >> num_cages;
    vector<Polygon> cages(num_cages);
    for (int i = 0; i < num_cages; i++)
    {
        int num_vertices_cage;
        cin >> num_vertices_cage;
        cages[i] = Polygon(num_vertices_cage);
        for (int j = 0; j < num_vertices_cage; j++)
        {
            cin >> x >> y;
            cages[i][j] = Point(x, y);
        }
    }

    // --- Pré-processamento ---
    vector<Point> nodes;
    nodes.push_back(start);
    vector<int> node_to_cage_map;
    node_to_cage_map.push_back(-1);

    for (int i = 0; i < num_cages; i++)
    {
        for (Point &vertex : cages[i])
        {
            nodes.push_back(vertex);
            node_to_cage_map.push_back(i);
        }
    }
    node_to_cage_map.push_back(-1);

    nodes.push_back(end);
    int end_index = nodes.size() - 1;
    int num_nodes = nodes.size();

    // Construir matriz de distâncias
    distances = vector<vector<double>>(num_nodes, vector<double>(num_nodes));
    for (int i = 0; i < num_nodes; i++)
    {
        for (int j = 0; j < num_nodes; j++)
        {
            distances[i][j] = distance(nodes[i], nodes[j]);
        }
    }

    // **ESSENCIAL**: Construir a matriz de interseções, cages_intersections armazena quais jaulas o segmento de reta entre u e v intercepta
    vector<vector<bitset<MAX_CAGES>>> cages_intersections(num_nodes, vector<bitset<MAX_CAGES>>(num_nodes)); // matrix that stores intersected cages(non visited) by segment
    for (int i = 0; i < num_nodes; i++)
    {
        for (int j = 0; j < num_nodes; j++)
        {
            vector<int> intersections;
            for (int k = 0; k < num_cages; k++)
            {
                if (is_intersection_segment_polygon(nodes[i], nodes[j], cages[k]))
                {
                    intersections.push_back(k);
                }
            }
            cages_intersections[i][j] = create_mask(intersections);
        }
    }
    long double avg_result = 0.0;
    vector<int> best_tour;

    for(int i = 0; i < 5; i++) {
        cout << "========================================================" << endl;
        cout << "Running simulated annealing iteration " << i + 1 << "..." << endl;
        // Executa o Simulated Annealing
        vector<int> current_tour = simulated_annealing(num_nodes,end_index,num_cages,node_to_cage_map,cages_intersections,nodes);
        
        double current_length = tour_length(current_tour);
        avg_result += current_length;
        cout << "Iteration " << i + 1 << " cost: " << current_length << endl;

        if (i == 0 || current_length < tour_length(best_tour)) {
            best_tour = current_tour;
        }
        cout << "========================================================" << endl << endl;

    }
    
    //-----------------------Print the final result and route-----------------------//
    cout << endl << "=== Final results ===" << endl;
    cout << "Approximate optimal distance: "<< tour_length(best_tour) <<"\n";
    cout << "Route: [";
    for(int i=0; i<best_tour.size(); i++){
        nodes[best_tour[i]].show();
        if(i!=best_tour.size()-1) cout << ", ";
    }
    cout << "]" << endl;

    cout << "Average cost over iterations: " << avg_result / 5.0 << endl;



    return 0;
}