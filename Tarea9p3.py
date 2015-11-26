#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt


def montecarlo(flujo_i, error_i, flujo_z, error_z):
    N_montecarlo = 100000  # En tarea 8 se usó N = 100k
    pendiente = np.zeros(N_montecarlo)
    coef_posicion = np.zeros(N_montecarlo)
    for i in range(N_montecarlo):
        r = np.random.normal(0, 1, size=len(flujo_i))
        muestra_i = flujo_i + error_i * r
        muestra_z = flujo_z + error_z * r
        pendiente[i], coef_posicion[i] = np.polyfit(muestra_i, muestra_z, 1)
    pendiente = np.sort(pendiente)
    coef_posicion = np.sort(coef_posicion)
    lim_bajo_pendiente = pendiente[int(N_montecarlo * 0.025)]
    lim_alto_pendiente = pendiente[int(N_montecarlo * 0.975)]
    lim_bajo_posicion = coef_posicion[int(N_montecarlo * 0.025)]
    lim_alto_posicion = coef_posicion[int(N_montecarlo * 0.975)]
    conf_pendiente = (lim_bajo_pendiente, lim_alto_pendiente)
    conf_posicion = (lim_bajo_posicion, lim_alto_posicion)
    return (conf_pendiente, conf_posicion)


# Setup

np.random.seed(123)
banda_i = np.loadtxt('data/DR9Q.dat', usecols=[80]) * 3.631
err_banda_i = np.loadtxt('data/DR9Q.dat', usecols=[81]) * 3.631
banda_z = np.loadtxt('data/DR9Q.dat', usecols=[82]) * 3.631
err_banda_z = np.loadtxt('data/DR9Q.dat', usecols=[83]) * 3.631

p = np.polyfit(banda_i, banda_z, 1)
y = np.zeros(len(banda_i))
for i in range(len(banda_i)):
    y[i] = p[0] * banda_i[i] + p[1]

confianza = montecarlo(banda_i, err_banda_i, banda_z, err_banda_z)
print " Pendiente y posicion"
print p[0]
print p[1]
print " Confianza pendiente y posicion"
print confianza[0]
print confianza[1]

plt.figure(1)
plt.clf
plt.plot(banda_i, banda_z, 'k*')
plt.plot(banda_i, y, 'g', label='$H_0$')
plt.title('$Flujo \ banda \ I \ vs \ banda \ Z $', fontsize=12)
plt.xlabel("$ Flujo \ Banda  \ i [J \ year \ 10^{-6}]$", fontsize=12)
plt.ylabel("$ Flujo \ Banda  \ z [J \ year \ 10^{-6}]$", fontsize=12)
plt.savefig("p3.png")
plt.show()
plt.draw()
