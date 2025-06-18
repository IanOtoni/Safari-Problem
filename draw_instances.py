import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def ler_instancia_txt(nome_arquivo):
    with open(nome_arquivo, 'r') as f:
        linhas = [l.strip() for l in f if l.strip()]

    nome = linhas[0]
    # Conversão para float em vez de int
    partida = tuple(map(float, linhas[1].split()))
    chegada = tuple(map(float, linhas[2].split()))

    partes = linhas[3].split()
    n_vert = int(partes[0])
    # Coords como float
    poligono = [ (float(partes[i]), float(partes[i+1]))
                 for i in range(1, 2*n_vert+1, 2) ]

    qtd = int(linhas[4])
    gaiolas = []
    idx = 5
    for _ in range(qtd):
        partes = linhas[idx].split()
        n = int(partes[0])
        coords = [ (float(partes[i]), float(partes[i+1]))
                   for i in range(1, 2*n+1, 2) ]
        gaiolas.append(coords)
        idx += 1

    return nome, partida, chegada, poligono, gaiolas


def desenhar_instancia(nome, partida, chegada, poligono, gaiolas, rota=None, saida_png=None):
    fig, ax = plt.subplots(figsize=(10, 10))

    # Paleta savana
    fundo = '#f4e3b2'; principal = '#c49e52'; gai = '#8b4513'
    cor_part = '#1b5e20'; cor_cheg = '#0d47a1'

    ax.set_facecolor(fundo)
    # área principal
    patch = Polygon(poligono, closed=True, facecolor=principal, edgecolor='black', alpha=0.6)
    ax.add_patch(patch)
    # gaiolas
    for coords in gaiolas:
        p = Polygon(coords, closed=True, facecolor=gai, edgecolor='black', alpha=0.7)
        ax.add_patch(p)
    # pontos de partida e chegada
    ax.plot(partida[0], partida[1], 'o', markersize=10, color=cor_part, label='Partida')
    ax.plot(chegada[0], chegada[1], 'o', markersize=10, color=cor_cheg, label='Chegada')

    # rota (se fornecida)
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
    filename = 'instances/inst8.txt'  # altere conforme necessidade
    nome, partida, chegada, poligono, gaiolas = ler_instancia_txt(filename)
    # Exemplo de rota (pode ser float)
    rota = [(100, 850), (100, 800), (100, 650), (200, 600), (330, 800), (250, 890), (400, 880), (605, 865), (430, 690), (530, 610), (630, 660), (690, 720), (900, 850), (950, 800), (1060, 830), (1130, 760), (1080, 660), (975, 595), (790, 270), (1100, 150)]
    desenhar_instancia(
        nome, partida, chegada, poligono, gaiolas,
        rota=rota,
        saida_png=filename.replace('.txt', '.png')
    )
