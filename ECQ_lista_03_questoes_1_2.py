import numpy as np
import matplotlib.pyplot as plt

def fugacidade(P,T, b,R=8.314, ): # o valor de b é o valor do coeficiente de vdw do He
    return P*np.exp(b*P/(R*T))

P = np.linspace(0.1, 500, 1000)
T = 298.15
b = 0.0238


# Gerar um intervalo contínuo de pressões (de 0,1 a 500 atm)
P = np.linspace(0.1, 500, 1000)

# Calcular a fugacidade para cada pressão
f = fugacidade(P, T, b)

# Plotar
plt.figure(figsize=(8, 5))
plt.plot(P, f, linewidth=2, color='#1f77b4', label='Modelo de van der Waals (a=0)')



# Configurações
plt.xlabel('Pressão $P$ (atm)', fontsize=12)
plt.ylabel('Fugacidade $f$ (atm)', fontsize=12)
plt.title(f'Fugacidade do hélio a T = {T} K (modelo com a=0)', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.65)
plt.legend()
plt.tight_layout()
plt.show()


###################################################################################

b_amonia = 0.0371
T = 300
P = 101325*10
f_amonia = fugacidade(P, T, b_amonia)
print(f"Valor de f estimado pelo modelo desconsiderando forças atrativas: {f_amonia:.6f}")

gamma = np.exp(b_amonia*P/(8.314*T))



####################################################################################



P2 = np.linspace(0, 550, 500)  # Pressão de 0 a 550 bar
Z = 1 - 0.000612 * P2 + 2.661e-6 * P2**2 - 2.390e-9 * P2**3 - 1.077e-13 * P2**4

plt.figure(figsize=(8, 5))
plt.plot(P2, Z, label='Z(P)', color='darkblue', linewidth=2)
plt.title('Comportamento de Z em função da pressão', fontsize=12)
plt.xlabel('Pressão P (bar)', fontsize=10)
plt.ylabel('Fator de Compressibilidade Z', fontsize=10)
plt.xlim(0, 550)
plt.ylim(0.88, 1.16)
plt.grid(True, linestyle='-.', alpha=0.45)
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()


P3 = np.linspace(0.001, 550, 500)
z_minus_1_over_p = -0.000612 + 2.661e-6 * P3 - 2.390e-9 * P3**2 - 1.077e-13 * P3**3

# Cálculo do ponto onde Z = 1 (ou seja, (Z-1)/P = 0)
coeffs = [-1.077e-13, -2.390e-9, 2.661e-6, -0.000612]
roots = np.roots(coeffs)
P_z1 = [r.real for r in roots if np.isreal(r) and 0 <= r.real <= 550][0]

plt.figure(figsize=(8, 5))
plt.plot(P3, z_minus_1_over_p, color='darkred', linewidth=2,
         label=r'$\frac{Z-1}{P} = -6.12\times 10^{-4} + 2.661\times 10^{-6}P - 2.390\times 10^{-9}P^2 - 1.077\times 10^{-13}P^3$')

plt.axhline(-0.000612, color='navy', linestyle='-.', 
            label=r'$\lim_{P \to 0}\frac{Z-1}{P} = B = -6.12\times 10^{-4}\mathrm{bar}^{-1}$')

plt.plot(P_z1, 0, 'ko', label=f'Z=1 em P = {P_z1:.1f} bar')

plt.title('Análise de $(Z - 1)/P$', fontsize=12)
plt.xlabel('Pressão P (bar)', fontsize=10)
plt.ylabel(r'$\frac{Z-1}{P}\quad (\mathrm{bar}^{-1})$', fontsize=10)
plt.xlim(0, 550)
plt.ylim(-0.00065, 0.00015)
plt.grid(True, linestyle=':', alpha=0.5)
plt.legend(loc='lower right', fontsize=8)
plt.tight_layout()
plt.savefig('figura_3_analise_z.png', dpi=300)
plt.show()
