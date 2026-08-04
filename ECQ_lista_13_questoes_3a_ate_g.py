#questão 3. (a)

P = [23.50, 23.06, 22.86, 22.40, 21.95, 21.48]   
P_estrela = 23.756                                # pressão de referência (solvente puro)
a = []                                            # lista vazia para armazenar as atividades
for j in P:                                       #
    _ = j/P_estrela                               
    a.append(_)

mensagem = "As atividades calculadas na escala pressão, são tem como resultado:"                                    
print(mensagem)
for i in a:
    print(str(i)) 
    
print("\n# ================================================================= #\n")
    
#questao 3. (b)


#\ln a_1(m) = -\frac{\phi(m)\,m}{55{.}509}.

# phi = - 55.509 \cdot \ln{a} 

import numpy as np

def coeficiente_osmotico(atividade,m):
    return - (55.509/m) * np.log(atividade)

phi = []

m = [0.6,1.6,2.0,3.0,4.0,5.0]

for ativ, molali in zip(a,m):
    _ = coeficiente_osmotico(ativ,molali)
    phi.append(_)

mensagem_b = "Os coeficientes osmóticos para as respectivas atividades são: "

print(mensagem_b)
for p in phi:
    phi_2 = str(p)
    print(phi_2)
    
print("\n# ================================================================= #\n")


#questao 3. (c):
    


import matplotlib.pyplot as plt


a_array = np.array(a)
phi_array = np.array(phi)
m_array = np.array(m)

grau = 5
coef = np.polyfit(m_array, phi_array, grau)   
p = np.poly1d(coef)               # função polinomial


phi_pred = p(m_array)

# R^{2}
ss_res = np.sum((phi_array - phi_pred) ** 2)          # soma dos quadrados dos resíduos
ss_tot = np.sum((phi_array - np.mean(phi_array)) ** 2)      # soma total dos quadrados
r2 = 1 - (ss_res / ss_tot)                      # coeficiente de determinação

#  curva
m_smooth = np.linspace(min(m_array), max(m_array), 200)
phi_smooth = p(m_smooth)



# =============================================================================
# Gráfico

plt.figure(figsize=(8, 5))

# Pontos experimentais
plt.plot(m, phi_array, 'h', color='#FF5A36', markersize=8, label='Dados obtidos')

# Curva de ajuste com R² na legenda
eq_str = f'$ ({coef[0]:.3f})x^5 + ({coef[1]:.3f})x^4 + ({coef[2]:.3f})x^3 + ({coef[3]:.3f})x^2 \n$ $+ ({coef[4]:.3f})x + ({coef[5]:.3f})$'
plt.plot(m_smooth, phi_smooth, ':', color='#A100BA', linewidth=1.75,label=f'{eq_str}, $R^2$ = {r2:.4f}')




# Configurações
plt.xlabel('Molalidade $m$', fontsize=12)
plt.ylabel(r'$\phi(m)$', fontsize=12)
plt.title('$\phi$ vs $m$ com ajuste polinomial', fontsize=14)
plt.legend(loc='best', fontsize=11)
plt.grid(True, linestyle=':', alpha=0.65)
plt.tight_layout()
plt.show()

#questao 3 d
P_array = np.array(P)
a1 = a_array

# Ajuste linear (grau 1) para phi(m)

coef_phi = np.polyfit(m_array, phi_array, 1)
p_phi = np.poly1d(coef_phi)

# Cálculo de f(m) = (phi(m)-1)/m para os pontos dados
f = (phi_array - 1) / m_array

# Extrapolação linear para f(0) usando os dois primeiros pontos
m_fit = m[:2]
f_fit = f[:2]
coef_f0 = np.polyfit(m_fit, f_fit, 1)
a, b = coef_f0[0], coef_f0[1]

f0 = b

print(f"Valor de f(0) estimado pelo modelo linear: {f0:.6f}")

# Agora montamos os arrays com f(0) incluso
m_integ = np.concatenate(([0], m))
f_integ = np.concatenate(([f0], f))

# Regra do trapézio acumulativa

J = np.zeros_like(m)
for i in range(len(m)):
    # integral de 0 até m[i] usando os pontos até i+1 (pois m_integ[0]=0)
    # Podemos usar numpy.trapz mas com cuidado; faremos manual
    integral = 0.0
    for k in range(1, i+2):  # k=1 corresponde ao intervalo [0, m[0]], etc.
        integral += (f_integ[k] + f_integ[k-1]) / 2 * (m_integ[k] - m_integ[k-1])
    J[i] = integral

# (d) ln(gamma_2)
ln_gamma_2 = phi_array - 1 - J


# Exibir tabela comparativa com o PDF
print("m      phi        J          ln(gamma2)   gamma2")
for mi, ph, j, lng in zip(m, phi, J, ln_gamma_2):
    print(f"{mi:5.2f} {ph:10.6f} {j:12.8f} {lng:12.8f}")
###########################################################################

# extraolacao dos primeiros  pontos

# Criar uma reta para visualização (de m = 0 até m = 5.5)
m_line = np.linspace(0, 5.5, 100)
f_line = a * m_line + b

# Configurar o gráfico
plt.figure(figsize=(8, 5))
plt.scatter(m_array, f, color='#0B8FDB')
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
#plt.xlim(-0.2, 5.5)
#plt.ylim(-0.01, 0.04)

# Exibir o gráfico
plt.show()

################################

from scipy.integrate import cumulative_trapezoid   # Regra do Trapézio cumulativa


from scipy.interpolate import CubicSpline







# adicionando o ponto phi = 1 e m = 0
m_array = np.concatenate(([0], m_array))
phi_array = np.concatenate(([1], phi_array))

grau = 6
coef = np.polyfit(m_array, phi_array, grau)
p_phi = np.poly1d(coef)               # função polinomial para phi(m)

# Limite fundamental: f(0) = lim_{m->0} (phi(m) - 1)/m = phi'(0)
limite_f_zero = f0             # derivada do polinômio em m=0


m_array = np.delete(m_array, 0)
phi_array = np.delete(phi_array,0)

# =============================================================================
# 4. INTEGRAÇÃO PELA REGRA DO TRAPÉZIO 
# =============================================================================
# cumulative_trapezoid com initial=0 retorna a integral acumulada 
# desde o primeiro ponto, ou seja, F(m) = ∫₀ᵐ f(m') dm'
f_array = (phi_array -1)/m_array
integral_not_dense = cumulative_trapezoid(f_array, x=m_array, initial=0)

# Interpola a integral para os valores exatos de m fornecidos

integral_at_m = J


# =============================================================================
# 5. Cálculo de ln(gamma_2) via Gibbs-Duhem
#    ln(γ₂) = φ(m) - 1 - ∫₀ᵐ (φ(m') - 1)/m' dm'
# =============================================================================
ln_gamma_2 = (phi_array - 1) - integral_at_m




# =============================================================================
# 6. Exibição dos resultados (Tabela)
# =============================================================================
print("\n" + "=" * 60)
print("Resultados do cálculo de ln(γ₂) via Gibbs-Duhem (Método do Trapézio)")
print("=" * 60)
print(f"{'m (mol/kg)':<12} {'φ(m)':<12} {'∫ f dm':<16} {'ln(γ₂)':<12}")
print("-" * 60)
for mi, phii, inti, lng in zip(m_array, phi_array, integral_at_m, ln_gamma_2):
    print(f"{mi:<12.2f} {phii:<12.6f} {inti:<16.8f} {lng:<12.8f}")

# =============================================================================
# Gráfico de ln(gamma_2) vs m
# =============================================================================
plt.figure(figsize=(10, 6))

# Curva suave (para o gráfico contínuo)
# Cria uma spline que passa exatamente por todos os pontos (m_array, ln_gamma_2)
spline = CubicSpline(m_array, ln_gamma_2, bc_type='natural')  # ou 'clamped'

# Cria uma malha densa para a curva suave (incluindo m=0, se quiser)
m_smooth = np.linspace(0, max(m_array), 500)
ln_gamma_smooth = spline(m_smooth)

# Agora plote
plt.figure(figsize=(8, 5))
plt.plot(m_smooth, ln_gamma_smooth, '-', color='#0B8FDB', linewidth=2, label='Curva suave (spline)')
plt.plot(m_array, ln_gamma_2, 'h', color='#DB0B45', markersize=9, label='Dados calculados')
plt.xlabel('Molalidade $m$ (mol kg$^{-1}$)', fontsize=12)
plt.ylabel(r'$\ln \gamma_{2,m}(m)$', fontsize=12)
plt.title('Coeficiente de Atividade do Soluto via Gibbs-Duhem', fontsize=14)
plt.legend(loc='best')
plt.grid(True, linestyle=':', alpha=0.65)
plt.tight_layout()
plt.show()



#questao 3 e

ln_gamma_2 = phi_array - 1 - integral_at_m

gamma_2 = np.exp(ln_gamma_2) 

mensagem_e = r"Os coeficientes de atividade para as respectivas valores de $ln{\gamma}$ são: "
print(str(mensagem_e))
print(str(gamma_2))
