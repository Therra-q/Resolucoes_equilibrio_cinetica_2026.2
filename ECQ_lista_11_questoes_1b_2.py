import numpy as np
import matplotlib.pyplot as plt
import math  

# =============================================================================
# exercício 1 b
# =============================================================================
def exercicio_1():
    
    x_1 = np.linspace(0, 1, 1000)
    x_2 = 1 - x_1

    # Pressões reais (não ideais)
    P_1 = 120 * x_1 * np.exp(0.20 * x_2**2 + 0.10 * x_2**3)
    P_2 = 140 * x_2 * np.exp(0.35 * x_1**2 - 0.10 * x_1**3)

    # Pressões ideais (lei de Raoult)
    P_1_ideal = 120 * x_1
    P_2_ideal = 140 * x_2

    # Determinação dos parâmetros
    P_1_estrela = 120.0                     # x1 = 1, x2 = 0
    P_2_estrela = 140.0                     # x2 = 1, x1 = 0
    k_H1 = 120 * np.exp(0.20 * 1**2 + 0.10 * 1**3)   # limite x1 -> 0 (x2 -> 1)
    k_H2 = 140 * np.exp(0.35 * 1**2 - 0.10 * 1**3)   # limite x2 -> 0 (x1 -> 1)

    print("=== Exercício 1 ===")
    print(f"P1* = {P_1_estrela:.2f} torr")
    print(f"P2* = {P_2_estrela:.2f} torr")
    print(f"kH,1 = {k_H1:.2f} torr")
    print(f"kH,2 = {k_H2:.2f} torr\n")

    # Gráficos
    plt.figure(figsize=(12, 5))

    # P1 versus x1
    plt.subplot(1, 2, 1)
    plt.plot(x_1, P_1, label=r'$P_1$ real')
    plt.plot(x_1, P_1_ideal, ':', label=r'$P_1$ ideal (Raoult)')
    plt.xlabel(r'$x_1$ (fração molar do componente 1)')
    plt.ylabel(r'Pressão de vapor $P_1$ (torr)')
    plt.legend()
    plt.title('P₁ vs x₁')

    # P2 versus x1
    plt.subplot(1, 2, 2)
    plt.plot(x_1, P_2, label=r'$P_2$ real')
    plt.plot(x_1, P_2_ideal, ':', label=r'$P_2$ ideal (Raoult)')
    plt.xlabel(r'$x_1$ (fração molar do componente 1)')
    plt.ylabel(r'Pressão de vapor $P_2$ (torr)')
    plt.legend()
    plt.title('P₂ vs x₁')

    plt.tight_layout()
    plt.show()

# =============================================================================
# Exercicio 2
# =============================================================================
def exercicio_2():
    # Temperatura em °C (80 a 110 °C) convertida para Kelvin
    T_C = np.linspace(80, 110, 1000)
    T = T_C + 273.15
    
    # Para benzeno: ln(P*) = -3856.56/T + 17.551
    # Para tolueno: ln(P*) = -4514.6/T + 18.397

    A_benz = 17.551    
    A_tolu = 18.397

    P_benz = np.exp(-3856.56 / T + A_benz)
    P_tolu = np.exp(-4514.6  / T + A_tolu)

    # P_total = 760 torr
    x1 = np.zeros_like(T)
    y1 = np.zeros_like(T)
    valido = np.zeros_like(T, dtype=bool)

    for i, (Pb, Pt) in enumerate(zip(P_benz, P_tolu)):
        if Pb != Pt:   # evita divisão por zero
            xi = (760 - Pt) / (Pb - Pt)
            if 0.0 <= xi <= 1.0:
                valido[i] = True
                x1[i] = xi
                y1[i] = xi * Pb / 760.0

    # Filtra apenas os pontos válidos
    T_valid = T[valido]
    x1_valid = x1[valido]
    y1_valid = y1[valido]

    # Gráfico T-x-y
    plt.figure(figsize=(8, 6))
    plt.plot(x1_valid, T_valid, label='Líquido (x₁)')
    plt.plot(y1_valid, T_valid, label='Vapor (y₁)')
    plt.xlabel('Fração molar de benzeno (x₁ ou y₁)')
    plt.ylabel('Temperatura (K)')
    plt.title('Diagrama T-x-y para benzeno–tolueno a 760 torr')
    plt.legend()
    plt.grid(True)
    plt.show()
