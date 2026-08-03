#questão 3. (a)

P = [23.50, 23.06, 22.86, 22.40, 21.95, 21.48]   
P_estrela = 23.756                                # pressão de referência (solvente puro)
a = []                                            # lista vazia para armazenar as atividades
for j in P:                                       # j recebe cada valor da lista P
    _ = j/P_estrela                               
    a.append(_)

mensagem = "As atividades calculadas na escala pressão, são tem como resultado:"  
torr = " torr"                                    
print(mensagem)
for i in a:
    print(str(i) + torr) 
    
print("\n# ================================================================= #\n")
    
#questao 3. (b)


#\ln a_1(m) = -\frac{\phi(m)\,m}{55{.}509}.

# phi = - 55.509 \cdot \ln{a} 

import numpy as np

def coeficiente_osmotico(atividade):
    return - 55.509 * np.log(atividade)

phi = []

for ativ in a:
    _ = coeficiente_osmotico(ativ)
    phi.append(_)

mensagem_b = "Os coeficientes osmóticos para as respectivas atividades são: "

print(mensagem_b)
for p in phi:
    phi_2 = str(p)
    print(phi_2)
    
print("\n# ================================================================= #\n")


#questao 3. (c):
    


import matplotlib.pyplot as plt


a = np.array([0.7898, 0.7750, 0.7682, 0.7528, 0.7377, 0.7219])
phi = np.array([13.1018, 14.1510, 14.6345, 15.7629, 16.8894, 18.0909])


grau = 1
coef = np.polyfit(a, phi, grau)   
p = np.poly1d(coef)               # função polinomial


phi_pred = p(a)

# R^{2}
ss_res = np.sum((phi - phi_pred) ** 2)          # soma dos quadrados dos resíduos
ss_tot = np.sum((phi - np.mean(phi)) ** 2)      # soma total dos quadrados
r2 = 1 - (ss_res / ss_tot)                      # coeficiente de determinação

#  curva
a_smooth = np.linspace(min(a), max(a), 200)
phi_smooth = p(a_smooth)



# =============================================================================
# Gráfico

plt.figure(figsize=(8, 5))

# Pontos experimentais
plt.plot(a, phi, 'h', color='#FF5A36', markersize=8, label='Dados obtidos')

# Curva de ajuste com R² na legenda
eq_str = f'$ {coef[0]:.4f}x + {coef[1]:.4f}$'
plt.plot(a_smooth, phi_smooth, ':', color='#A100BA', linewidth=1.75,label=f'{eq_str}, $R^2$ = {r2:.4f}')




# Configurações
plt.xlabel('Atividade $a_{1(m)}$', fontsize=12)
plt.ylabel(r'$\phi(m)$', fontsize=12)
plt.title('$\phi$ vs $a$ com ajuste polinomial', fontsize=14)
plt.legend(loc='best', fontsize=11)
plt.grid(True, linestyle=':', alpha=0.65)
plt.tight_layout()
plt.show()

#questao 3 d



# extraolacao dos primeiros  pontos

m = np.array([0.60, 1.60, 2.00, 3.00, 4.00, 5.00])
f = np.array([0.0078, 0.01988, 0.0333, 0.0293, 0.0246, 0.0237])

# Selecionar apenas os dois primeiros pontos (m = 0.6 e m = 1.6)
m_fit = m[:2]
f_fit = f[:2]

# Ajuste linear: f(m) = a * m + b
# Usamos polyfit de grau 1 para obter coeficientes (a, b)
coef = np.polyfit(m_fit, f_fit, 1)
a, b = coef[0], coef[1]

# O valor de f(0) é o intercepto b
f0 = b
print(f"Valor de f(0) estimado pelo modelo linear: {f0:.6f}")

# Criar uma reta para visualização (de m = 0 até m = 5.5)
m_line = np.linspace(0, 5.5, 100)
f_line = a * m_line + b

# Configurar o gráfico
plt.figure(figsize=(8, 5))
plt.scatter(m, f, color='#0B8FDB', label='Dados experimentais')
plt.scatter(m_fit, f_fit, color='#DB0B45', edgecolors='black', s=100,
            label='Pontos usados no ajuste')
plt.plot(m_line, f_line, 'r--', label=f'Reta: f(m) = {a:.4f}·m + {b:.4f}')

# Destacar o ponto f(0)
plt.plot(0, f0, 'go', markersize=10, label=f'f(0) = {f0:.6f}')

# Configurações do gráfico
plt.xlabel('$m$ (mol kg$^{-1}$)')
plt.ylabel('$f(m)$')
plt.title('Ajuste linear usando os dois primeiros pontos')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)
plt.xlim(-0.2, 5.5)
plt.ylim(-0.01, 0.04)

# Exibir o gráfico
plt.show()

################################
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid   # Regra do Trapézio cumulativa
from scipy.interpolate import interp1d

# =============================================================================
# 1. Dados CORRETOS (P_estrela = 23.756 torr, conforme o enunciado)
# =============================================================================
P_estrela = 23.756                     # pressão do solvente puro (torr)
P = np.array([23.50, 23.06, 22.86, 22.40, 21.95, 21.48])
m = np.array([0.6, 1.6, 2.0, 3.0, 4.0, 5.0])

# Cálculo das atividades e coeficientes osmóticos
a1 = P / P_estrela
phi = -55.509 * np.log(a1) / m        # phi(m) = -55.509 * ln(a1) / m

print("Atividades (a1):", a1)
print("Coeficientes osmóticos (phi):", phi)

# =============================================================================
# 2. Ajuste polinomial para obter uma curva suave e o limite em m=0
# =============================================================================
grau = 3
coef = np.polyfit(m, phi, grau)
p_phi = np.poly1d(coef)               # função polinomial para phi(m)

# Limite fundamental: f(0) = lim_{m->0} (phi(m) - 1)/m = phi'(0)
dp_phi = p_phi.deriv()
limite_f_zero = dp_phi(0)             # derivada do polinômio em m=0

# =============================================================================
# 3. Malha densa para integração (de 0 até o máximo de m)
# =============================================================================
m_dense = np.linspace(0, max(m), 1000)   # 1000 pontos de 0 a 5 mol/kg

# Calcula o integrando f(m) = (phi(m) - 1) / m
f_dense = np.zeros_like(m_dense)
f_dense[0] = limite_f_zero               # remove a singularidade em m=0
f_dense[1:] = (p_phi(m_dense[1:]) - 1) / m_dense[1:]

# =============================================================================
# 4. INTEGRAÇÃO PELA REGRA DO TRAPÉZIO (substituindo o Simpson)
# =============================================================================
# cumulative_trapezoid com initial=0 retorna a integral acumulada 
# desde o primeiro ponto, ou seja, F(m) = ∫₀ᵐ f(m') dm'
integral_dense = cumulative_trapezoid(f_dense, x=m_dense, initial=0)

# Interpola a integral para os valores exatos de m fornecidos
integral_func = interp1d(m_dense, integral_dense, kind='cubic')
integral_at_m = integral_func(m)

# =============================================================================
# 5. Cálculo de ln(gamma_2) via Gibbs-Duhem
#    ln(γ₂) = φ(m) - 1 - ∫₀ᵐ (φ(m') - 1)/m' dm'
# =============================================================================
ln_gamma_2 = phi - 1 - integral_at_m

# =============================================================================
# 6. Exibição dos resultados (Tabela)
# =============================================================================
print("\n" + "=" * 60)
print("Resultados do cálculo de ln(γ₂) via Gibbs-Duhem (Método do Trapézio)")
print("=" * 60)
print(f"{'m (mol/kg)':<12} {'φ(m)':<12} {'∫ f dm':<16} {'ln(γ₂)':<12}")
print("-" * 60)
for mi, phii, inti, lng in zip(m, phi, integral_at_m, ln_gamma_2):
    print(f"{mi:<12.2f} {phii:<12.6f} {inti:<16.8f} {lng:<12.8f}")

# =============================================================================
# 7. Gráfico de ln(gamma_2) vs m
# =============================================================================
plt.figure(figsize=(10, 6))

# Curva suave (para o gráfico contínuo)
ln_gamma_2_dense = p_phi(m_dense) - 1 - integral_dense

plt.plot(m_dense, ln_gamma_2_dense, '-.', color='#0B8FDB', linewidth=1.75,
         label='$\ln \gamma_2(m)$ (Trapézio)')
plt.plot(m, ln_gamma_2, 'h', color='#DB0B45', markersize=9,
         label='Dados calculados nos pontos experimentais')

plt.xlabel('Molalidade $m$ (mol kg$^{-1}$)', fontsize=12)
plt.ylabel(r'$\ln \gamma_{2,m}(m)$', fontsize=12)
plt.title('Coeficiente de Atividade do Soluto via Gibbs-Duhem (Regra do Trapézio)', fontsize=14)
plt.legend(loc='best', fontsize=11)
plt.grid(True, linestyle=':', alpha=0.65)
plt.tight_layout()
plt.show()

#questao 3 e

ln_gamma = [-52.683310 , -63.939917 , -66.444040, -71.049928 , -74.323972, -76.798113]
gamma = []

def gamma_exp(ln_gamma):
    return np.e(ln_gamma)

for j in ln_gamma:
    _ = gamma_exp(j)
    gamma.append(_)

mensagem_e = r"Os coeficientes de atividade para as respectivas valores de $ln{\gamma}$ são: "
print(mensagem_e)
for i in gamma:
    g = str(i)
    print(g)
