import matplotlib.pyplot as plt

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

plt.xlabel('Iterações')
plt.title('Euler com T =0.05')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper left')
plt.show()

print("Valor Final Método de Euler CA(0.3) = " + str(v_ca[-1]))
print("Valor Final Método de Euler CB(0.3) = " + str(v_cb[-1]))


#rungekutta

v_ca_rk = rungekutta(f_a, 0, 0.3, 0, 0.05)[0]
v_cb_rk = rungekutta(f_b, 0, 0.3, 0, 0.05)[0]
v_ha_rk = rungekutta(f_a, 0, 0.3, 0, 0.05)[1]
v_hb_rk = rungekutta(f_b, 0, 0.3, 0, 0.05)[1]

graf_ca, = plt.plot(v_ca_rk, v_ha_rk)
graf_cb, = plt.plot(v_cb_rk, v_hb_rk)

plt.xlabel('Iterações')
plt.title('Runge kutta com T =0.05')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper left')
plt.show()

print("Valor Final Método de Runge Kutta CA(0.3) = " + str(v_ca_rk[-1]))
print("Valor Final Método de Runge Kutta CB(0.3) = " + str(v_cb_rk[-1]))

#T = 0.03

#euler

v_ca = euler(f_a, 0, 0.3, 0, 0.03)[0]
v_cb = euler(f_b, 0, 0.3, 0, 0.03)[0]
v_ha = euler(f_a, 0, 0.3, 0, 0.03)[1]
v_hb = euler(f_b, 0, 0.3, 0, 0.03)[1]

graf_ca, = plt.plot(v_ca, v_ha)
graf_cb, = plt.plot(v_cb, v_hb)


plt.xlabel('Iterações')
plt.title('Euler com T =0.03')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper right')
plt.show()


#rungekutta

v_ca_rk = rungekutta(f_a, 0, 0.3, 0, 0.03)[0]
v_cb_rk = rungekutta(f_b, 0, 0.3, 0, 0.03)[0]
v_ha_rk = rungekutta(f_a, 0, 0.3, 0, 0.03)[1]
v_hb_rk = rungekutta(f_b, 0, 0.3, 0, 0.03)[1]

graf_ca, = plt.plot(v_ca_rk, v_ha_rk)
graf_cb, = plt.plot(v_cb_rk, v_hb_rk)

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

plt.xlabel('Iterações')
plt.title('Euler com T =0.01')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper right')
plt.show()

print("Valor Final Método de Euler CA(0.3) = " + str(v_ca[-1]))
print("Valor Final Método de Euler CB(0.3) = " + str(v_cb[-1]))


#rungekutta
v_ca_rk = rungekutta(f_a, 0, 0.3, 0, 0.01)[0]
v_cb_rk = rungekutta(f_b, 0, 0.3, 0, 0.01)[0]
v_ha_rk = rungekutta(f_a, 0, 0.3, 0, 0.01)[1]
v_hb_rk = rungekutta(f_b, 0, 0.3, 0, 0.01)[1]

graf_ca, = plt.plot(v_ca_rk, v_ha_rk)
graf_cb, = plt.plot(v_cb_rk, v_hb_rk)
plt.xlabel('Iterações')
plt.title('Runge kutta com T =0.01')
plt.legend((graf_ca, graf_cb), ('CA(t)', 'CB(t)'), loc='upper left')
plt.show()

print("Valor Final Método de Runge Kutta CA(0.3) = " + str(v_ca_rk[-1]))
print("Valor Final Método de Runge Kutta CB(0.3) = " + str(v_cb_rk[-1]))


#1.3