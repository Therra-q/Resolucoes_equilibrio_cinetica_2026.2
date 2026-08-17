import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve  # para encontrar as raízes

def exercicio_3():

    T = 290   # K
    R = 8.314 # J mol^-1 K^-1
    alpha = 6000 # J mol^-1 (parâmetro da solução regular)
    w = alpha / (R * T)    

    x_1 = np.linspace(0.001, 0.999, 1000)   # evita log(0)
    x_2 = 1 - x_1
    
    g_mix = x_1 * np.log(x_1) + x_2 * np.log(x_2) + w * x_1 * x_2


    def derivada_g(x):
        # derivada de g(x) = ln(x/(1-x)) + w*(1 - 2x)
        return np.log(x / (1 - x)) + w * (1 - 2 * x)

    # Chute inicial para a raiz entre 0 e 0.5
    chute = 0.01
    x_1_linha = fsolve(derivada_g, chute)[0]       # x1' (fase pobre em 1)
    x_1_duas_linhas = 1 - x_1_linha              # x1'' = 1 - x1'

    # Valor da função nos pontos de tangência (a reta é horizontal)
    y_tangente = x_1_linha * np.log(x_1_linha) + (1 - x_1_linha) * np.log(1 - x_1_linha) + w * x_1_linha * (1 - x_1_linha)

    # Exibir resultados
    print("=== Exercício 3 (d) ===")
    print(f"Temperatura T = {T} K")
    print(f"Parâmetro w = alpha/(RT) = {w:.4f}")
    print(f"Composição da fase 1 (x_1')  = {x_1_linha:.6f}")
    print(f"Composição da fase 2 (x_1'') = {x_1_duas_linhas:.6f}")
    print(f"Valor da tangente comum: y = {y_tangente:.6f}")

    # =========================================================================
    # Construção do gráfico
    # =========================================================================
    plt.figure(figsize=(10, 6))

    plt.plot(x_1, g_mix, 'b-', linewidth=2.5,
             label=r'$\Delta_{\text{mix}}g/RT$')


    plt.axhline(y=y_tangente, color='red', linestyle='--', linewidth=2,
                label='Reta tangente comum (horizontal)')


    plt.axvline(x=x_1_linha, ymin=0, ymax=1, color='green', linestyle=':', linewidth=2)
    plt.axvline(x=x_1_duas_linhas, ymin=0, ymax=1, color='green', linestyle=':', linewidth=2)


    plt.plot(x_1_linha, y_tangente, 'go', markersize=10,
             label=f'$x_1\'$ = {x_1_linha:.4f}')
    plt.plot(x_1_duas_linhas, y_tangente, 'ro', markersize=10,
             label=f'$x_1\'\'$ = {x_1_duas_linhas:.4f}')


    plt.xlabel(r'Fração molar do componente 1 ($x_1$)', fontsize=12)
    plt.ylabel(r'$\Delta_{\text{mix}}g \, / \, RT$', fontsize=12)
    plt.title(f'Solução Regular a T = {T} K  ($w = \\alpha/RT = {w:.3f}$)', fontsize=14)
    plt.legend(loc='best', fontsize=11)
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.xlim(0, 1)

    y_min = np.min(g_mix) - 0.05
    y_max = np.max(g_mix) + 0.05
    plt.ylim(y_min, y_max)

    plt.tight_layout()
    plt.show()

# Executa a função
exercicio_3()
