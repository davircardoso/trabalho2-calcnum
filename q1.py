import matplotlib.pyplot as plt

#1.1 método de Euler


def f_a(ca ,cb):
    ca_0 = 5.1
    u = 14.19
    k1 = k2 = 14.75
    k3 = 2.28
    fca = u * (ca_0 - ca) - k1*ca - k3*(ca ** 2)
    return fca

def f_b(ca,cb):
    u = 14.19
    k1 = k2 = 14.75
    k3 = 2.28
    fcb = - u * cb + k1 * ca - k2 * cb
    return fcb


def euler(y0, x0, t_max, h, f):
    t_e = []
    y_e = []

    ca = y0
    cb = x0

    while cb < t_max + h:
        t_e.append(ca)
        y_e.append(cb)

        ca += h * f(ca, cb)
        cb += h

    return t_e, y_e


def rungekutta(y0, x0, t_max, h, f):
    t_rk = []
    y_rk = []

    ca_rk = y0
    cb_rk = x0

    while cb_rk < t_max + h:
        t_rk.append(ca_rk)
        y_rk.append(cb_rk)

        f1 = f(cb_rk, ca_rk)
        f2 = f(cb_rk + 0.5*h*f1, ca_rk + 0.5*h)
        ca_rk += h * f2
        cb_rk += h

    return t_rk, y_rk

#1.2 resoluções

f_a(0,0)
f_b(0,0)

#T = 0.05

#euler
v_ca = euler(0, 0, 0.3, 0.05, f_a)[0]
v_cb = euler(0, 0, 0.3, 0.05, f_b)[0]
v_h = euler(0, 0, 0.3, 0.05, f_a)[1]

graf_ca, = plt.plot(v_h, v_ca)
graf_cb, = plt.plot(v_h, v_cb)

plt.xlabel('Iterações')
plt.title('Euler com T =0.05')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper left')
plt.show()

print("Valor Final Método de Euler CA(0.3) = " + str(v_ca[-1]))
print("Valor Final Método de Euler CB(0.3) = " + str(v_cb[-1]))


#rungekutta

v_ca_rk = rungekutta(0, 0, 0.3, 0.05, f_a)[0]
v_cb_rk = rungekutta(0, 0, 0.3, 0.05, f_b)[0]
v_h_rk = rungekutta(0, 0, 0.3, 0.05, f_a)[1]

graf_ca, = plt.plot(v_h_rk, v_ca_rk)
graf_cb, = plt.plot(v_h_rk, v_cb_rk)

plt.xlabel('Iterações')
plt.title('Runge kutta com T =0.05')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper left')
plt.show()

print("Valor Final Método de Runge Kutta CA(0.3) = " + str(v_ca_rk[-1]))
print("Valor Final Método de Runge Kutta CB(0.3) = " + str(v_cb_rk[-1]))

#T = 0.03

#euler

v_ca = euler(0, 0, 0.3, 0.03, f_a)[0]
v_cb = euler(0, 0, 0.3, 0.03, f_b)[0]
v_h = euler(0, 0, 0.3, 0.03, f_a)[1]

graf_ca, = plt.plot(v_h, v_ca)
graf_cb, = plt.plot(v_h, v_cb)

plt.xlabel('Iterações')
plt.title('Euler com T =0.03')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper right')
plt.show()


#rungekutta

v_ca_rk = rungekutta(0, 0, 0.3, 0.03, f_a)[0]
v_cb_rk = rungekutta(0, 0, 0.3, 0.03, f_b)[0]
v_h_rk = rungekutta(0, 0, 0.3, 0.03, f_a)[1]

graf_ca, = plt.plot(v_h_rk, v_ca_rk)
graf_cb, = plt.plot(v_h_rk, v_cb_rk)

plt.xlabel('Iterações')
plt.title('Runge kutta com T =0.03')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper left')
plt.show()


#T = 0.01

#euler
v_ca = euler(0, 0, 0.3, 0.01, f_a)[0]
v_cb = euler(0, 0, 0.3, 0.01, f_b)[0]
v_h = euler(0, 0, 0.3, 0.01, f_a)[1]

graf_ca, = plt.plot(v_h, v_ca)
graf_cb, = plt.plot(v_h, v_cb)

plt.xlabel('Iterações')
plt.title('Euler com T =0.01')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper right')
plt.show()

print("Valor Final Método de Euler CA(0.3) = " + str(v_ca[-1]))
print("Valor Final Método de Euler CB(0.3) = " + str(v_cb[-1]))


#rungekutta

v_ca_rk = rungekutta(0, 0, 0.3, 0.01, f_a)[0]
v_cb_rk = rungekutta(0, 0, 0.3, 0.01, f_b)[0]
v_h_rk = rungekutta(0, 0, 0.3, 0.01, f_a)[1]

graf_ca, = plt.plot(v_h_rk, v_ca_rk)
graf_cb, = plt.plot(v_h_rk, v_cb_rk)

plt.xlabel('Iterações')
plt.title('Runge kutta com T =0.01')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper left')
plt.show()

print("Valor Final Método de Runge Kutta CA(0.3) = " + str(v_ca_rk[-1]))
print("Valor Final Método de Runge Kutta CB(0.3) = " + str(v_cb_rk[-1]))


