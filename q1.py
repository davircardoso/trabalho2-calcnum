import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np

#1.1 método de Euler

u = 14.19
ca_0 = 5.1
k1 = 14.75
k2 = 14.75
k3 = 2.28

def f_a(ca, _, __):
    fca = u * (ca_0 - ca) - k1*ca - k3*(ca ** 2)
    return fca

def euler(f, ta, tb, i, h):
    x = ta
    y = tb
    k = 0
    list_y = [i]
    list_x = [ta]
    while x < tb:
        y = y + h*f(y, k, h)
        x += h
        k += 1

        list_y.append(y)
        list_x.append(x)

    return list_x, list_y, k

def f_b(cb, k, h):
    ca = euler(f_a, 0, 0.3, 0, h)
    cay = ca[1]
    fcb = -u * cb + k1 * cay[k] - k2 * cb
    return fcb



def rungekutta(f, ta, tb, i, h):
    x = ta
    y = tb
    k = 0
    list_y = [i]
    list_x = [ta]
    while x <= tb:
        z1 = h * f(y, k, h)
        z2 = h * f(y + z1, k, h)
        y = y + (z1 + z2) / 2
        x = x + h
        k = k + 1

        list_y.append(y)
        list_x.append(x)

    return list_x, list_y, k


#1.2 resoluções


#T = 0.05

#euler
v_ca = euler(f_a, 0, 0.3, 0, 0.05)[0]
v_cb = euler(f_b, 0, 0.3, 0, 0.05)[0]
v_ha = euler(f_a, 0, 0.3, 0, 0.05)[1]
v_hb = euler(f_b, 0, 0.3, 0, 0.05)[1]

graf_ca, = plt.plot(v_ca, v_ha)
graf_cb, = plt.plot(v_cb, v_hb)

plt.grid("true")
plt.ylabel("concentração")
plt.xlabel('Iterações')
plt.title('Euler com T =0.05')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper left')
plt.show()

print("Valor Final Método de Euler CA(0.3) = " + str(v_ha[-1]))
print("Valor Final Método de Euler CB(0.3) = " + str(v_hb[-1]))


#rungekutta

v_ca_rk = rungekutta(f_a, 0, 0.3, 0, 0.05)[0]
v_cb_rk = rungekutta(f_b, 0, 0.3, 0, 0.05)[0]
v_ha_rk = rungekutta(f_a, 0, 0.3, 0, 0.05)[1]
v_hb_rk = rungekutta(f_b, 0, 0.3, 0, 0.05)[1]

graf_ca, = plt.plot(v_ca_rk, v_ha_rk)
graf_cb, = plt.plot(v_cb_rk, v_hb_rk)

plt.grid("true")
plt.ylabel("concentração")
plt.xlabel('Iterações')
plt.title('Runge kutta com T =0.05')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper left')
plt.show()

print("Valor Final Método de Runge Kutta CA(0.3) = " + str(v_ha_rk[-1]))
print("Valor Final Método de Runge Kutta CB(0.3) = " + str(v_hb_rk[-1]))

#T = 0.03

#euler

v_ca = euler(f_a, 0, 0.3, 0, 0.03)[0]
v_cb = euler(f_b, 0, 0.3, 0, 0.03)[0]
v_ha = euler(f_a, 0, 0.3, 0, 0.03)[1]
v_hb = euler(f_b, 0, 0.3, 0, 0.03)[1]

graf_ca, = plt.plot(v_ca, v_ha)
graf_cb, = plt.plot(v_cb, v_hb)

plt.grid("true")
plt.ylabel("concentração")
plt.xlabel('Iterações')
plt.title('Euler com T =0.03')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper right')
plt.show()

print("Valor Final Método de Euler CA(0.3) = " + str(v_ha_rk[-1]))
print("Valor Final Método de Euler CB(0.3) = " + str(v_hb_rk[-1]))

#rungekutta

v_ca_rk = rungekutta(f_a, 0, 0.3, 0, 0.03)[0]
v_cb_rk = rungekutta(f_b, 0, 0.3, 0, 0.03)[0]
v_ha_rk = rungekutta(f_a, 0, 0.3, 0, 0.03)[1]
v_hb_rk = rungekutta(f_b, 0, 0.3, 0, 0.03)[1]

print("Valor Final Método de Runge Kutta CA(0.3) = " + str(v_ha_rk[-1]))
print("Valor Final Método de Runge Kutta CB(0.3) = " + str(v_hb_rk[-1]))

graf_ca, = plt.plot(v_ca_rk, v_ha_rk)
graf_cb, = plt.plot(v_cb_rk, v_hb_rk)

plt.grid("true")
plt.ylabel("concentração")
plt.xlabel('Iterações')
plt.title('Runge kutta com T =0.03')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper left')
plt.show()


#T = 0.01

#euler
v_ca = euler(f_a, 0, 0.3, 0, 0.01)[0]
v_cb = euler(f_b, 0, 0.3, 0, 0.01)[0]
v_ha = euler(f_a, 0, 0.3, 0, 0.01)[1]
v_hb = euler(f_b, 0, 0.3, 0, 0.01)[1]

graf_ca, = plt.plot(v_ca, v_ha)
graf_cb, = plt.plot(v_cb, v_hb)

plt.grid("true")
plt.ylabel("concentração")
plt.xlabel('Iterações')
plt.title('Euler com T =0.01')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper right')
plt.show()

print("Valor Final Método de Euler CA(0.3) = " + str(v_ha[-1]))
print("Valor Final Método de Euler CB(0.3) = " + str(v_hb[-1]))


#rungekutta
v_ca_rk = rungekutta(f_a, 0, 0.3, 0, 0.01)[0]
v_cb_rk = rungekutta(f_b, 0, 0.3, 0, 0.01)[0]
v_ha_rk = rungekutta(f_a, 0, 0.3, 0, 0.01)[1]
v_hb_rk = rungekutta(f_b, 0, 0.3, 0, 0.01)[1]

graf_ca, = plt.plot(v_ca_rk, v_ha_rk)
graf_cb, = plt.plot(v_cb_rk, v_hb_rk)

plt.grid("true")
plt.ylabel("concentração")
plt.xlabel('Iterações')
plt.title('Runge kutta com T =0.01')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper left')
plt.show()

print("Valor Final Método de Runge Kutta CA(0.3) = " + str(v_ha_rk[-1]))
print("Valor Final Método de Runge Kutta CB(0.3) = " + str(v_hb_rk[-1]))


#1.3

#Solve ODE

t = np.arange(0.0, 0.3, 0.03)
y0 = [0, 0]
__ = 0
ode = odeint(f_a, y0, t, args = (__,))

#errro em relação à euler

euler_vca = euler(f_a, 0, 0.3, 0, 0.03)[1][-1]
ode1 = odeint(f_a, y0, t, args = (__,))[-1][1]

error = (euler_vca - ode1)/ode1 * 100

print(F"O erro de Euler em relação ao Solver de ODE é de {error} %")

v_ca = euler(f_a, 0, 0.3, 0, 0.03)[0]
v_cb = euler(f_b, 0, 0.3, 0, 0.03)[0]
v_ha = euler(f_a, 0, 0.3, 0, 0.03)[1]
v_hb = euler(f_b, 0, 0.3, 0, 0.03)[1]

graf_ca, = plt.plot(v_ca, v_ha)
graf_solver = plt.plot(t, ode)

plt.title("Gráfico utilizando Solver de ODE x Euler, t = 0.03")
plt.grid(True)
plt.legend((graf_ca, graf_solver), ('Euler', 'Solver de ODE'), loc='upper right')
plt.xlabel('iterações')
plt.ylabel('concentração')
plt.show()

#errro em relação à runge kutta
rk_vca = rungekutta(f_a, 0, 0.3, 0, 0.03)[1][-1]
ode2 = odeint(f_a, y0, t, args = (__,))[-1][1]

error = (ode2 - rk_vca)/ode1 * 100

print(F"O erro de Runge Kutta em relação ao Solver de ODE é de {error} %")

v_ca_rk = rungekutta(f_a, 0, 0.3, 0, 0.03)[0]
v_cb_rk = rungekutta(f_b, 0, 0.3, 0, 0.03)[0]
v_ha_rk = rungekutta(f_a, 0, 0.3, 0, 0.03)[1]
v_hb_rk = rungekutta(f_b, 0, 0.3, 0, 0.03)[1]

graf_ca, = plt.plot(v_ca_rk, v_ha_rk)
graf_solver = plt.plot(t,ode)

plt.title("Gráfico utilizando Solver de Od x Runge Kutta, t = 0.03")
plt.grid(True)
plt.legend((graf_ca, graf_solver), ('Runge Kutta', 'Solver de ODE'), loc='upper right')
plt.xlabel('iterações')
plt.ylabel('concentração')
plt.show()
