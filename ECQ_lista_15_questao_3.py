import numpy as np
def m(C, rho, M):
    return (1000*C)/( (1000*rho) - (C * M))

C_NaCl = 0.01
M_NaCl = 58.44
rho_H2O = 1.00
molalidade_NaCl = m(C_NaCl, rho_H2O, rho_H2O)

F_ionica = (1/2)*( (1* molalidade_NaCl)  + (1* molalidade_NaCl)  )

ln_gamma = -1.173*(1*1)*np.sqrt(F_ionica)

gamma = np.exp(ln_gamma)
mensagem = "O logaritmo natural do coeficiente de atividade é: "

mensagem_2 = "E o valor do coeficiente de atividade em si é: "
print("\n" + "=" * 60)
print(mensagem + str(ln_gamma))
print("\n" + mensagem_2 + str(gamma))
print("=" * 60)
