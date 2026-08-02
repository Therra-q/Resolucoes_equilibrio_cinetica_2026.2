#questão 3. (a)

P = [23.50, 23.06, 22.86, 22.40, 21.95, 21.48]   
P_estrela = 29.756                                # pressão de referência (solvente puro)
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

from scipy.integrate import cumulative_simpson
from scipy.interpolate import interp1d


m = np.array([0.6, 1.6, 2.0, 3.0, 4.0, 5.0])



grau = 3
coef = np.polyfit(m, phi, grau)
p_phi = np.poly1d(coef)          # função polinomial para phi



# Derivada para calcular o limite em m -> 0: f(0) = phi'(0)
dp_phi = p_phi.deriv()
limite_f_zero = dp_phi(0)        # valor de (phi - 1)/m quando m -> 0


#  Malha densa para integração (de 0 até o max(m))

m_dense = np.linspace(0, max(m), 1000)   # 1000 pontos de 0 a 5 mol/kg

# Calcula o integrando f(m) = (phi(m) - 1) / m
f_dense = np.zeros_like(m_dense)
# Para m = 0, usa o limite (remove a singularidade)
f_dense[0] = limite_f_zero
# Para m > 0, calcula normalmente
f_dense[1:] = (p_phi(m_dense[1:]) - 1) / m_dense[1:]

# Integração cumulativa (Simpson) para obter F(m) = ∫₀ᵐ f(m') dm'

integral_dense = cumulative_simpson(f_dense, x=m_dense, initial=0)

# Interpola a integral para os valores exatos de m fornecidos
integral_func = interp1d(m_dense, integral_dense, kind='cubic')
integral_at_m = integral_func(m)


# Cálculo de ln(gamma_2) usando a equação de Gibbs-Duhem
#    \ln(\gamma_2) = \phi(m) - 1 - \int (\phi - 1)/m dm

ln_gamma_2 = phi - 1 - integral_at_m

# =============================================================================
# 5. Exibição dos resultados
# =============================================================================
print("=" * 50)
print("Resultados do cálculo de ln(γ₂) via Gibbs-Duhem")
print("=" * 50)
print(f"{'m (mol/kg)':<12} {'φ(m)':<12} {'∫ f dm':<12} {'ln(γ₂)':<12}")
print("-" * 50)
for mi, phii, inti, lng in zip(m, phi, integral_at_m, ln_gamma_2):
    print(f"{mi:<12.2f} {phii:<12.4f} {inti:<12.6f} {lng:<12.6f}")

# =============================================================================
# 6. Gráfico de ln(gamma_2) vs m
# =============================================================================
plt.figure(figsize=(10, 6))

# Curva suave para visualização (opcional)
ln_gamma_2_dense = p_phi(m_dense) - 1 - integral_dense

plt.plot(m_dense, ln_gamma_2_dense, '-.', color='#0B8FDB', linewidth=1.75,
         label='$\ln \gamma_2(m)$')
plt.plot(m, ln_gamma_2, 'h', color='#DB0B45', markersize=9,
         label='Dados calculados')

plt.xlabel('Molalidade $m$ (mol kg$^{-1}$)', fontsize=12)
plt.ylabel(r'$\ln \gamma_{2,m}(m)$', fontsize=12)
plt.title('Coeficiente de Atividade do Soluto via Gibbs-Duhem', fontsize=14)
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
