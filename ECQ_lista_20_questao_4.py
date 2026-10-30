"""
Lista 20 - Exercício 4(b): cálculos e gráficos.

(iii) Método de Hitchcock: E' = E + (2RT/F) ln m - (2RT/F) ln(10) A sqrt(m)  vs  m
(iv)  Lei limite de Debye-Hückel: E + (2RT/F) ln m  vs  sqrt(m)

Também imprime os valores do Exercício 1.

Requer: matplotlib  (pip install matplotlib)
"""
from math import log, sqrt

import matplotlib.pyplot as plt

# Constantes
R = 8.314462618      # J/(K mol)
F = 96485.33212      # C/mol
T = 298.15           # K
A = 0.509            # constante de Debye-Hückel a 25 °C, (mol/kg)^(-1/2)
k2 = 2 * R * T / F   # 2RT/F em V


def ajuste_linear(x, y):
    """Mínimos quadrados y = a + b x. Retorna (a, b, R2, desvio padrão de a)."""
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((xi - mx) ** 2 for xi in x)
    sxy = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    b = sxy / sxx
    a = my - b * mx
    ssr = sum((yi - (a + b * xi)) ** 2 for xi, yi in zip(x, y))
    sst = sum((yi - my) ** 2 for yi in y)
    sa = sqrt(ssr / (n - 2)) * sqrt(1 / n + mx ** 2 / sxx)
    return a, b, 1 - ssr / sst, sa


# ---------------------------------------------------------------- Exercício 1
E0_1 = 70241 / (2 * F)
c = 0.070
gamma = 10 ** (-A * sqrt(c))
E_1 = E0_1 + k2 * log(gamma * c)
print("Exercício 1")
print(f"  E°  = {E0_1:.5f} V")
print(f"  γ±  = {gamma:.4f}")
print(f"  E   = {E_1:.5f} V\n")

# ---------------------------------------------------------- Exercício 4(b)
m_mmol = [1.6077, 3.0769, 5.0403, 7.6938, 10.9474]          # mmol/kg
E = [0.60080, 0.56825, 0.54366, 0.52267, 0.50532]           # V
m = [x / 1000 for x in m_mmol]                              # mol/kg
raiz_m = [sqrt(x) for x in m]

y_DH = [Ei + k2 * log(mi) for mi, Ei in zip(m, E)]          # E + (2RT/F) ln m
y_H = [yi - k2 * log(10) * A * sqrt(mi) for mi, yi in zip(m, y_DH)]  # E'

print("Exercício 4(b)")
print(f"  2RT/F = {k2:.6f} V ; (2RT/F) ln10 A = {k2 * log(10) * A:.6f} V")
print("  m (mmol/kg)  sqrt(m)   E (V)    E+(2RT/F)ln m   E' (V)")
for row in zip(m_mmol, raiz_m, E, y_DH, y_H):
    print("  {:10.4f}  {:.5f}  {:.5f}  {:.5f}         {:.5f}".format(*row))

aH, bH, R2H, saH = ajuste_linear(m, y_H)
aD, bD, R2D, saD = ajuste_linear(raiz_m, y_DH)
print(f"\n  (iii) Hitchcock: E' = {aH:.6f} {bH:+.5f} m   R² = {R2H:.4f}")
print(f"        E° = {aH:.6f} ± {saH:.6f} V")
print(f"  (iv)  Lei limite: y = {aD:.5f} {bD:+.5f} sqrt(m)   R² = {R2D:.4f}")
print(f"        E° = {aD:.5f} V")

# ------------------------------------------------------- Gráfico (iii)
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(m_mmol, y_H, "o", color="tab:blue", label="dados")
x_fit = [0, 12]
ax.plot(x_fit, [aH + bH * x / 1000 for x in x_fit], "-", color="tab:red",
        label=f"$\\mathcal{{E}}' = {aH:.6f} {bH:+.5f}\\,m$;  $R^2 = {R2H:.4f}$")
ax.set_xlabel("$m$ (mmol/kg)")
ax.set_ylabel("$\\mathcal{E}'$ (V)")
ax.set_title("Método de Hitchcock: $\\mathcal{E}'$ em função de $m$")
ax.set_xlim(0, 12)
ax.ticklabel_format(axis="y", useOffset=False)
ax.grid(True)
ax.legend()
fig.tight_layout()
fig.savefig("grafico_hitchcock.png", dpi=300)

# ------------------------------------------------------- Gráfico (iv)
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(raiz_m, y_DH, "o", color="tab:blue", label="dados")
x_fit = [0, 0.12]
ax.plot(x_fit, [aD + bD * x for x in x_fit], "-", color="tab:red",
        label=f"ajuste: ${aD:.5f} {bD:+.5f}\\sqrt{{m}}$;  $R^2 = {R2D:.4f}$")
ax.plot(x_fit, [aD + k2 * log(10) * A * x for x in x_fit], "--", color="black",
        label=f"inclinação teórica da lei limite ({k2 * log(10) * A:.4f} V)")
ax.set_xlabel("$\\sqrt{m/m^\\circ}$ (adimensional)")
ax.set_ylabel("$\\mathcal{E} + \\frac{2RT}{F}\\ln(m/m^\\circ)$ (V)")
ax.set_title("Extrapolação pela lei limite de Debye–Hückel")
ax.set_xlim(0, 0.12)
ax.ticklabel_format(axis="y", useOffset=False)
ax.grid(True)
ax.legend(loc="lower right")
fig.tight_layout()
fig.savefig("grafico_debye_huckel.png", dpi=300)

plt.show()
