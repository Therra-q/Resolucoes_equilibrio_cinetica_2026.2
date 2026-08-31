import numpy as np
import matplotlib.pyplot as plt





def pressao_antoine (A,B,C,T):
    return A - (B)/(T+C)


A = 8.07313

B = 1730.63

C = 233.426

# Gerar um intervalo contínuo de pressões (de 0,1 a 500 atm)
T = np.linspace(0, 500, 1000)

# Calcular a fugacidade para cada pressão
P = np.log10(pressao_antoine(A, B, C, T))

# Plotar
plt.figure(figsize=(8, 5))
plt.plot(T, P, linewidth=2, color='#1f77b4', label='Curva líquido-vapor')



# Configurações
plt.xlabel('Temperatura $T$ (K)', fontsize=12)
plt.ylabel('Pressão $P$ (atm)', fontsize=12)
plt.title('Curva líquido-vapor', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.65)
plt.legend()
plt.tight_layout()
plt.show()


