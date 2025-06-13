#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <random>
#include <numeric>
#include <algorithm>
#include <limits>
#include "geometria.h"

// Estrutura para manter os dados de uma gaiola
struct Gaiola {
    std::string nome;
    Poligono vertices;
};

// Estrutura para manter todos os dados da instância do problema
struct Instancia {
    std::string nome;
    Ponto ponto_partida;
    Ponto ponto_chegada;
    Poligono poligono_principal;
    std::vector<Gaiola> gaiolas;
};

Instancia carregarInstancia() {
    Instancia inst;

    int n_vertices_poligono, n_gaiolas, n_vertices_gaiola;

    std::cin >> inst.nome;
    std::cin >> inst.ponto_partida.x >> inst.ponto_partida.y;
    std::cin >> inst.ponto_chegada.x >> inst.ponto_chegada.y;

    std::cin >> n_vertices_poligono;
    for (int i = 0; i < n_vertices_poligono; ++i) {
        Ponto p;
        std::cin >> p.x >> p.y;
        inst.poligono_principal.push_back(p);
    }

    std::cin >> n_gaiolas;
    for (int i = 0; i < n_gaiolas; ++i) {
        Gaiola g;
        std::cin >> g.nome >> n_vertices_gaiola;
        for (int j = 0; j < n_vertices_gaiola; ++j) {
            Ponto p;
            std::cin >> p.x >> p.y;
            g.vertices.push_back(p);
        }
        inst.gaiolas.push_back(g);
    }

    return inst;
}


// Classe principal que encapsula o algoritmo ACO
class SolucionadorACO {
public:
    SolucionadorACO(const Instancia& inst, int num_formigas, int num_iteracoes, double alfa, double beta, double rho)
        : m_instancia(inst),
          NUM_FORMIGAS(num_formigas),
          NUM_ITERACOES(num_iteracoes),
          ALFA(alfa),
          BETA(beta),
          RHO(rho),
          Q(100.0) // Fator de deposição de feromônio
    {
        // Mapeia todos os pontos (vértices) para um único vetor para fácil acesso
        // 0: ponto de partida
        // 1 a N: vértices das gaiolas
        // N+1: ponto de chegada
        m_nos.push_back(m_instancia.ponto_partida);
        for (size_t i = 0; i < m_instancia.gaiolas.size(); ++i) {
            m_mapa_gaiola_para_no_inicio[i] = m_nos.size();
            for (const auto& vertice : m_instancia.gaiolas[i].vertices) {
                m_nos.push_back(vertice);
            }
        }
        m_nos.push_back(m_instancia.ponto_chegada);
        m_num_nos = m_nos.size();

        // Inicializa a matriz de feromônios com um valor pequeno e positivo
        m_feromonios.resize(m_num_nos, std::vector<double>(m_num_nos, 0.1));

        // Inicializa o gerador de números aleatórios
        m_gen.seed(std::random_device()());
    }

    void resolver() {
        m_melhor_tour.clear();
        m_melhor_distancia = std::numeric_limits<double>::max();

        for (int iter = 0; iter < NUM_ITERACOES; ++iter) {
            std::vector<std::vector<int>> tours_formigas;
            std::vector<double> distancias_tours;

            // Para cada formiga, construa uma solução (tour)
            for (int i = 0; i < NUM_FORMIGAS; ++i) {
                auto tour = construirTour();
                if (!tour.empty()) {
                    tours_formigas.push_back(tour);
                    distancias_tours.push_back(calcularDistanciaTour(tour));
                }
            }

            // Evaporação de feromônio em todas as arestas
            evaporarFeromonios();
            
            // Deposição de feromônio com base nos tours construídos
            depositarFeromonios(tours_formigas, distancias_tours);

            // Atualiza o melhor tour encontrado até agora
            for (size_t i = 0; i < tours_formigas.size(); ++i) {
                if (distancias_tours[i] < m_melhor_distancia) {
                    m_melhor_distancia = distancias_tours[i];
                    m_melhor_tour = tours_formigas[i];
                }
            }

            std::cout << "Iteração " << iter + 1 << "/" << NUM_ITERACOES 
                      << ", Melhor Distância: " << m_melhor_distancia << std::endl;
        }

        std::cout << "\n--- OTIMIZAÇÃO CONCLUÍDA ---\n";
        std::cout << "Melhor distância encontrada: " << m_melhor_distancia << std::endl;
        std::cout << "Melhor rota (sequência de nós): ";
        for (int no_idx : m_melhor_tour) {
            std::cout << no_idx << " ";
        }
        std::cout << std::endl;
    }

private:
    // --- Membros da classe ---
    const Instancia& m_instancia;
    int NUM_FORMIGAS;
    int NUM_ITERACOES;
    double ALFA;  // Importância do feromônio
    double BETA;  // Importância da heurística (distância)
    double RHO;   // Taxa de evaporação do feromônio
    double Q;     // Constante de deposição de feromônio

    std::vector<Ponto> m_nos;
    std::vector<int> m_mapa_gaiola_para_no_inicio; // Mapeia índice da gaiola para índice do seu primeiro nó
    int m_num_nos;
    
    std::vector<std::vector<double>> m_feromonios;
    
    std::vector<int> m_melhor_tour;
    double m_melhor_distancia;

    std::mt19937 m_gen; // Gerador de números aleatórios


    // --- Métodos principais do ACO ---

    std::vector<int> construirTour() {
        std::vector<int> tour;
        tour.push_back(0); // Começa no nó de partida

        int no_atual_idx = 0;

        // Itera sobre cada gaiola que precisa ser visitada
        for (size_t i = 0; i < m_instancia.gaiolas.size(); ++i) {
            int proxima_gaiola_idx = i;
            
            // Obtém a lista de nós (vértices) candidatos na próxima gaiola
            std::vector<int> candidatos = getCandidatos(proxima_gaiola_idx, no_atual_idx);
            
            if (candidatos.empty()) {
                 // Não há caminho válido, a formiga se perdeu
                return {}; 
            }

            // Calcula a probabilidade de transição para cada nó candidato
            std::vector<double> probabilidades = calcularProbabilidades(no_atual_idx, candidatos);
            
            // Seleciona o próximo nó com base na roleta
            int proximo_no_idx = selecionarProximoNo(candidatos, probabilidades);

            tour.push_back(proximo_no_idx);
            no_atual_idx = proximo_no_idx;
        }

        tour.push_back(m_num_nos - 1); // Adiciona o nó de chegada
        return tour;
    }

    // Retorna os nós candidatos da próxima gaiola que são "visíveis" do nó atual
    std::vector<int> getCandidatos(int gaiola_idx, int no_atual_idx) {
        std::vector<int> candidatos;
        int no_inicio = m_mapa_gaiola_para_no_inicio[gaiola_idx];
        int no_fim = (gaiola_idx + 1 < m_instancia.gaiolas.size()) 
                        ? m_mapa_gaiola_para_no_inicio[gaiola_idx + 1] 
                        : m_num_nos - 1;

        for (int i = no_inicio; i < no_fim; ++i) {
            bool colide = false;
            // Verifica colisão com TODAS as outras gaiolas
            for (size_t j = 0; j < m_instancia.gaiolas.size(); ++j) {
                // Não precisa checar colisão com a gaiola de destino
                if (j == gaiola_idx) continue; 
                
                if (segmentoInterceptaPoligono(m_nos[no_atual_idx], m_nos[i], m_instancia.gaiolas[j].vertices)) {
                    colide = true;
                    break;
                }
            }
            if (!colide) {
                candidatos.push_back(i);
            }
        }
        return candidatos;
    }

    std::vector<double> calcularProbabilidades(int no_atual_idx, const std::vector<int>& candidatos) {
        std::vector<double> probabilidades;
        double soma_total = 0.0;

        for (int candidato_idx : candidatos) {
            double feromonio = m_feromonios[no_atual_idx][candidato_idx];
            double dist = distancia(m_nos[no_atual_idx], m_nos[candidato_idx]);
            double heuristica = 1.0 / dist; // Heurística é o inverso da distância

            double prob = std::pow(feromonio, ALFA) * std::pow(heuristica, BETA);
            probabilidades.push_back(prob);
            soma_total += prob;
        }

        // Normaliza as probabilidades
        if (soma_total > 0) {
            for (double& p : probabilidades) {
                p /= soma_total;
            }
        }
        return probabilidades;
    }

    // Seleção por Roleta
    int selecionarProximoNo(const std::vector<int>& candidatos, const std::vector<double>& probabilidades) {
        std::uniform_real_distribution<> dis(0.0, 1.0);
        double roleta = dis(m_gen);
        double soma_acumulada = 0.0;
        for (size_t i = 0; i < candidatos.size(); ++i) {
            soma_acumulada += probabilidades[i];
            if (roleta <= soma_acumulada) {
                return candidatos[i];
            }
        }
        return candidatos.back(); // Fallback
    }

    void evaporarFeromonios() {
        for (int i = 0; i < m_num_nos; ++i) {
            for (int j = 0; j < m_num_nos; ++j) {
                m_feromonios[i][j] *= (1.0 - RHO);
            }
        }
    }

    void depositarFeromonios(const std::vector<std::vector<int>>& tours, const std::vector<double>& distancias) {
        for (size_t i = 0; i < tours.size(); ++i) {
            const auto& tour = tours[i];
            double delta_feromonio = Q / distancias[i]; // Formigas de caminhos curtos depositam mais feromônio
            for (size_t j = 0; j < tour.size() - 1; ++j) {
                int no1 = tour[j];
                int no2 = tour[j + 1];
                m_feromonios[no1][no2] += delta_feromonio;
                m_feromonios[no2][no1] += delta_feromonio; // Aresta simétrica
            }
        }
    }

    double calcularDistanciaTour(const std::vector<int>& tour) {
        double d = 0.0;
        for (size_t i = 0; i < tour.size() - 1; ++i) {
            d += distancia(m_nos[tour[i]], m_nos[tour[i + 1]]);
        }
        return d;
    }
};


int main() {
    // --- PARÂMETROS DO ACO ---
    // Estes valores geralmente requerem ajuste fino para cada tipo de problema
    int NUM_FORMIGAS = 20;
    int NUM_ITERACOES = 100;
    double ALFA = 1.0;          // Peso do feromônio
    double BETA = 5.0;          // Peso da informação heurística (distância)
    double RHO = 0.5;           // Taxa de evaporação

    // Carrega a instância desejada
    Instancia inst = carregarInstancia();  // Sem caminho!
    std::cout << "Resolvendo instância: " << inst.nome << std::endl;

    // Cria e executa o solucionador
    SolucionadorACO solucionador(inst, NUM_FORMIGAS, NUM_ITERACOES, ALFA, BETA, RHO);
    solucionador.resolver();

    return 0;
}