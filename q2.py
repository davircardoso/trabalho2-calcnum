#Leitura de arquivos .mat
import numpy as np
import matplotlib.pyplot as plt
import scipy.io

u_signal = scipy.io.loadmat('u_signal.mat')['u_signal']
y_signal = scipy.io.loadmat('y_signal.mat')['y_signal']
u_signal_test = scipy.io.loadmat('u_signal_test.mat')['u_signal_test']
y_signal_test = scipy.io.loadmat('y_signal_test.mat')['y_signal_test']

def moving_Average(n,N,M,u_signal,y_signal,u_signal_test,y_signal_test):
#Moving Average com n atrasos
    print("______________________________________________________")
    print("\nMoving Average com",n,"atrasos:\n")

#criando a Matriz A
    A = np.zeros((N-n,n))
    for i in range(n,N):
        for j in range(n):
            u_signal_elem = i-j-1
            if(u_signal_elem >= n-1):
                A[i-n][j] = u_signal[u_signal_elem]

#criando a Matriz B
    B = np.zeros((N-n))
    for i in range(N-n-1):
        B[i] = y_signal[i+n-1]

    #resolvendo por minimos quadrados:

    xls = np.dot(np.linalg.inv(np.dot(np.transpose(A),A)),np.dot(np.transpose(A),B))
    y_minQuad = np.dot(A,xls)
    E = np.dot(np.transpose(y_minQuad-B),(y_minQuad-B))

    for i in range(len(xls)):
        print("W"+ str(i) + "=" , xls[i])
    print("\nE = ",E)

    #Plotando o gráfico que compara y_signal em laranja com o resultado obtido com xls em azul

    plt.title("Modelo Fir")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.plot(y_minQuad)
    plt.plot(y_signal)
    plt.legend((y_minQuad,y_signal),("yminquad","ysignal"),loc="upper right")
    plt.show()

    #Testando se o W está correto
    T = np.linalg.lstsq(A,B,rcond=None)
    print("\n\nW por linalg.lstsq dada para teste:\n")
    for i in range(len(T[0])):
        print("W"+ str(i) + "=" , T[0][i])

    #criando a Matriz A_test
    A_test = np.zeros((M-n,n))
    for i in range(n,M):
        for j in range(n):
            u_signal_elem = i-j-1
            if(u_signal_elem >=n-1):
                A_test[i-n][j] = u_signal_test[u_signal_elem]

    #criando a Matriz B_test
    B_test = np.zeros((M-n))
    for i in range(M-n-1):
        B_test[i] = y_signal_test[i+n-1]

    plt.plot(np.dot(A_test,xls))
    plt.plot(y_signal_test)
    plt.show()
    E = np.dot(np.transpose(np.dot(A_test,xls)-B_test),(np.dot(A_test,xls)-B_test))
    print("\nE = ",E)

#chamando a função Moving Average com atraso de 10
moving_Average(10,len(u_signal),len(u_signal_test),u_signal,y_signal,u_signal_test,y_signal_test)
moving_Average(20,len(u_signal),len(u_signal_test),u_signal,y_signal,u_signal_test,y_signal_test)


#criando matriz B, sinal Y de saida do sistema:
N = len(y_signal)
Y = np.zeros((N))
for i in range(N):
    Y[i] = y_signal[i]

#criando matriz A, phi transposto:
d = 4
PHI = np.zeros((N,d))
for i in range(N-1):
    #Definindo uma linha completa de PHI transposto:
    if i>=0: PHI[i][0] = u_signal[i]
    else:PHI[i][0] = 0
    if i-1>=0:PHI[i][1] = u_signal[i-1]
    else:PHI[i][1] = 0
    if i-1>=0:PHI[i][2] = Y[i-1]
    else:PHI[i][2] = 0
    if i-2>=0:PHI[i][3] = Y[i-2]
    else:PHI[i][3] = 0

w = np.dot(np.linalg.inv(np.dot(np.transpose(PHI),PHI)),np.dot(np.transpose(PHI),Y))
y_minQuad_ARX = np.dot(PHI,w)
E = np.dot(np.transpose(y_minQuad_ARX-Y),(y_minQuad_ARX-Y))

for i in range(len(w)):
    print("W"+ str(i) + "=" , w[i])
print("\nE Aw= ",E)

#Plotando o gráfico que compara y_signal em laranja com o resultado obtido com w em azul

plt.plot(y_minQuad_ARX)
plt.plot(Y)

M = len(y_signal_test)
Y_test = np.zeros((M))
for i in range(M):
    Y_test[i] = y_signal_test[i]

#Usando a matriz de teste
PHI_test = np.zeros((M,d))
for i in range(M-1):
#Definindo uma linha completa de PHI transposto:
    if i >= 0:PHI_test[i][0] = u_signal_test[i]
    else:PHI_test[i][0] = 0
    if i-1 >= 0:PHI_test[i][1] = u_signal_test[i-1]
    else:PHI_test[i][1] = 0
    if i-1 >= 0:PHI_test[i][2] = Y_test[i-1]
    else:PHI_test[i][2] = 0
    if i-2 >= 0:PHI_test[i][3] = Y_test[i-2]
    else:PHI_test[i][3] = 0

plt.plot(np.dot(PHI_test,w))
plt.plot(y_signal_test)
E = np.dot(np.transpose(np.dot(PHI_test,w)-Y_test),(np.dot(PHI_test,w)-Y_test))
print("\nE ARX= ",E)


N = len(u_signal)
n = 10
A = np.zeros((N-n,n))
for i in range(n,N):
    for j in range(n):
        u_signal_elem = i-j-1
        if(u_signal_elem >= n-1):
            A[i-n][j] = u_signal[u_signal_elem]