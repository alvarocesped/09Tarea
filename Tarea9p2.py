#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import (leastsq, curve_fit)


def hubble_1(H0, D):
    '''
    Modelo 1 propuesto por Hubble en 1929, donde H0 es la constante de Hubble
    en km / s / MPc.
    '''
    v = H0 * D
    return v


def hubble_2(H0, v):
    '''
    Versión modificada del modelo de Hubble.
    '''
    D = v / H0
    return D


def minimizar1(datos_d, H0):
    '''
    Esta función y la siguiente son para la simulación bootstrap
    '''
    return hubble_1(H0, datos_d)


def minimizar2(datos_v, H0):
    return hubble_2(H0, datos_v)


def bootstrap(datos_d, datos_v, H0_inicial):
    '''
    Se realiza un procedimiento bootstrap para que el rango de confianza
    quede por sobre el 95 por ciento para el cálculo de H0 óptimo.
    '''

    D = datos_d
    v = datos_v

    N = len(D)
    N_boot = int(N * np.log10(N)**2) * 100

    np.random.seed(123)
    H0_calculados = np.zeros(N_boot)

    for i in range(N_boot):
        r = np.random.randint(low=0, high=N, size=N)
        D_boot = D[r]
        v_boot = v[r]
        modelo_1_opt, modelo_1_var = curve_fit(minimizar1, D_boot, v_boot,
                                               H0_inicial)
        modelo_2_opt, modelo_2_var = curve_fit(minimizar2, v_boot, D_boot,
                                               H0_inicial)
        H0_promedio = (modelo_1_opt + modelo_2_opt) / 2.
        H0_calculados[i] = H0_promedio

    H0_sort = np.sort(H0_calculados)  # Valores ordenados

    # Confianza de un 95 por ciento
    lim_bajo = H0_sort[int(N_boot * 0.025)]
    lim_alto = H0_sort[int(N_boot * 0.975)]
    return [H0_calculados, lim_bajo, lim_alto]

# Setup

H0_inicial = 3
distancia = np.loadtxt("data/SNIa.dat", usecols=[-1])
velocidad = np.loadtxt("data/SNIa.dat", usecols=[-2])

D_arreglo = np.linspace(min(distancia), max(distancia), 100)
v_arreglo = np.linspace(min(velocidad), max(velocidad), 100)

# Encontrar H0 óptimo sin privilegiar modelo
H0_optimo_mod1, H0_cov_mod1 = curve_fit(minimizar1, distancia, velocidad,
                                        H0_inicial)
H0_optimo_mod2, H0_cov_mod2 = curve_fit(minimizar2, velocidad, distancia,
                                        H0_inicial)

H0_optimo_modprom = (H0_optimo_mod1 + H0_optimo_mod2) / 2.  # H0 optimo prom

print ' H0 modelo 1 = ', H0_optimo_mod1[0], '[km / s / Mpc]'
print ' H0 modelo 2 = ', H0_optimo_mod2[0], '[km / s / Mpc]'
print ' H0 promedio = ', H0_optimo_modprom[0], '[km / s / Mpc]'

# Intervalo de confianza
boot = bootstrap(distancia, velocidad, H0_inicial)
H0 = boot[0]
lim_bajo = boot[1]
lim_alto = boot[2]
print "Intervalo de confianza: [{}, {}]".format(lim_bajo, lim_alto)

plt.figure(1)
plt.clf
plt.plot(distancia, velocidad, '*', color='black', label='Datos utilizados')
plt.plot(D_arreglo, minimizar1(D_arreglo, H0_optimo_mod1), color='green',
         label='Modelo 1 $v = H_0 * D$')
plt.plot(minimizar2(v_arreglo, H0_optimo_mod2), v_arreglo, color='blue',
         label='Modelo 2 $D = v / H_0$')
plt.plot(D_arreglo, minimizar1(D_arreglo, H0_optimo_modprom), color='red',
         label='Promedio de ambos modelos')
plt.title('Velocidad de recesion "Nebulosas" segun la distancia (1929)',
          fontsize=12)
plt.xlabel('Distancia [Mpc]', fontsize=12)
plt.ylabel('Velocidad de recesion [km / s]', fontsize=12)
plt.legend(loc='lower right')
plt.savefig('p2.png')
plt.draw()
plt.show()
