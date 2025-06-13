import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def ler_instancia_txt(nome_arquivo):
    with open(nome_arquivo, 'r') as f:
        linhas = [l.strip() for l in f if l.strip()]

    nome = linhas[0]
    partida = tuple(map(int, linhas[1].split()))
    chegada = tuple(map(int, linhas[2].split()))

    partes = list(map(int, linhas[3].split()))
    n_vert = partes[0]
    poligono = [tuple(partes[i:i+2]) for i in range(1, 2*n_vert+1, 2)]

    qtd = int(linhas[4])
    gaiolas = []
    idx = 5
    for _ in range(qtd):
        partes = list(map(int, linhas[idx].split()))
        n = partes[0]
        coords = [tuple(partes[i:i+2]) for i in range(1, 2*n+1, 2)]
        gaiolas.append(coords)
        idx += 1

    return nome, partida, chegada, poligono, gaiolas


def desenhar_instancia(nome, partida, chegada, poligono, gaiolas, rota=None, saida_png=None):
    fig, ax = plt.subplots(figsize=(10, 10))

    # Paleta savana
    fundo = '#f4e3b2'; principal = '#c49e52'; gai = '#8b4513'
    cor_part = '#1b5e20'; cor_cheg = '#0d47a1'

    ax.set_facecolor(fundo)
    # principal
    patch = Polygon(poligono, closed=True, facecolor=principal, edgecolor='black', alpha=0.6)
    ax.add_patch(patch)
    # gaiolas
    for coords in gaiolas:
        p = Polygon(coords, closed=True, facecolor=gai, edgecolor='black', alpha=0.7)
        ax.add_patch(p)
    # pontos
    ax.plot(*partida, 'o', markersize=10, color=cor_part, label='Partida')
    ax.plot(*chegada, 'o', markersize=10, color=cor_cheg, label='Chegada')

    # rota
    if rota:
        xs, ys = zip(*rota)
        ax.plot(xs, ys, '-', linewidth=2, color='black', label='Rota')

    ax.set_title(nome, fontsize=14, fontweight='bold', color='#4e342e')
    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', color='gray', alpha=0.3)
    ax.legend()
    plt.autoscale(enable=True, axis='both', tight=True)
    if saida_png:
        plt.savefig(saida_png)
        print(f"Imagem salva em {saida_png}")
    else:
        plt.show()
    plt.close()


if __name__ == '__main__':
    # Exemplo de uso
    filename = 'instances/inst4.txt'  # altere conforme
    nome, partida, chegada, poligono, gaiolas = ler_instancia_txt(filename)
    # Defina aqui sua rota: lista de coordenadas
    rota = [(10, 50), (40, 60), (90, 40), (150, 50), (180, 40), (240, 40), (290, 50)]
    desenhar_instancia(
        nome, partida, chegada, poligono, gaiolas,
        rota=rota,
        saida_png=filename.replace('.txt', '_rota.png')
    )
