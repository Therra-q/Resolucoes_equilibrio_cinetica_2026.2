import matplotlib.pyplot as plt
import numpy as np
import statistics

# -----------------------1.(b)
def funcao_S_x(x,N_A = 9, N_B = 3, V_A = 3e-6, V_B = 3e-6):
    S_linha = (N_A*V_A)**(1/3)*x**(1/3) + (N_B*V_B)**(1/3)*(1-x)**(1/3)
    return S_linha



x_vals = np.linspace(0, 1, 400)  # 500 pontos para um gráfico suave

S_vals = funcao_S_x(x_vals)

plt.figure(figsize=(8, 5))
plt.plot(x_vals, S_vals, color='#ad0274', linewidth=2)

plt.xlabel(r'$x$', fontsize=12)
plt.ylabel(r"$S'$", fontsize=12)
plt.title('Entropia do sistema composto em função de $x$', fontsize=14)
plt.grid(True, linestyle='-.', alpha=0.78)
plt.xlim(0, 1)         
plt.tight_layout()




#-------------------------1. (c)

idx_max = np.argmax(S_vals)
x_max_aprox = x_vals[idx_max]
S_max_aprox = S_vals[idx_max]

print(f"Máximo aproximado: x = {x_max_aprox:.6f}, S' = {S_max_aprox:.6f}")
