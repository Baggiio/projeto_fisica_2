# Simulador: Partícula em uma caixa
# em orientacao a objetos
# 1. Caixa 1D – Determinação da função de onda quântica e outros parâmetros
# 2. Caixa 1D – Cálculo dos parâmetros da caixa e partícula, dada a função de onda

# Importando as bibliotecas

from math import *
from scipy import integrate
from matplotlib import pyplot as plt
import numpy as np

# Constantes

h = 4.136e-15  # Constante de Planck
hm = 6.626e-34  # Constante de Planck
c = 3e8  # Velocidade da luz
me = 9.11e-31  # Massa do elétron
mp = 1.67e-27  # Massa do próton
j = 6.242e18  # Constante de conversão de Joule para eV


# Definindo a classe Caixa
class Caixa:

    # Construtor
    def __init__(
        self, L=None, ni=None, nf=None, a=None, b=None, A=None, k=None, xp=None
    ):
        self.L = L
        self.ni = ni
        self.nf = nf
        self.a = a
        self.b = b
        self.A = A
        self.k = k
        self.xp = xp

    # Métodos

    # Função de onda

    def funcaoQuanticaInicial(self, ni, L):
        A = sqrt(2 / L)
        k1 = (ni * pi) / L
        return A, k1

    def funcaoQuanticaFinal(self, nf, L):
        A = sqrt(2 / L)
        k2 = (nf * pi) / L
        return A, k2

    # Energia

    def energiaQuanticaEletron(self, ni, nf, L):
        Eij = (ni**2 * hm**2) / (8 * me * L**2)
        Eiv = Eij * j
        Efj = (nf**2 * hm**2) / (8 * me * L**2)
        Efv = Efj * j
        return Eij, Eiv, Efj, Efv

    def energiaQuanticaProton(self, ni, nf, L):
        Eij = ((hm**2) / (8 * mp * L**2)) * ni**2
        Eiv = Eij * j
        Efj = ((hm**2) / (8 * mp * L**2)) * nf**2
        Efv = Efj * j
        return Eij, Eiv, Efj, Efv

    # Fóton

    def calculoFoton(self, Eiv, Efv):
        E = abs(Efv - Eiv)
        lambada = (h * c) / E
        f = E / h
        return E, lambada, f

    # Velocidade no nivel quantico inicial e final

    def velocidadeEletron(self, Eij, Efj):
        viE = sqrt((2 * Eij) / me)
        vfE = sqrt((2 * Efj) / me)
        return viE, vfE

    def velocidadeProton(self, Eij, Efj):
        viP = sqrt((2 * Eij) / mp)
        vfP = sqrt((2 * Efj) / mp)
        return viP, vfP

    # Comprimento de onda de De Broglie

    def comprimentoDeBroglieEletron(self, viE, vfE):
        lambadai = hm / (me * viE)
        lambadaf = hm / (me * vfE)
        return lambadai, lambadaf

    def comprimentoDeBroglieProton(self, viP, vfP):
        lambadai = hm / (mp * viP)
        lambadaf = hm / (mp * vfP)
        return lambadai, lambadaf

    # |𝜓(𝑥)|**2
    # Probabilidade da integral definida de a até b de |𝜓(𝑥)|**2

    # def funcaoEletron(self, A, k, x):
    # funcao = pow(A * sin(k * x),2)
    # return funcao

    def probabilidadeIntegralNi(self, A, k1, a, b):
        integrali = integrate.quad((lambda x: pow(A * sin(k1 * x), 2)), a, b)
        percent = integrali[0] * 100
        return percent

    def probabilidadeIntegralNf(self, A, k2, a, b):
        integralf = integrate.quad(lambda x: pow(A * sin(k2 * x), 2), a, b)
        percent = integralf[0] * 100
        return percent

    # Largura da caixa

    def larguraCaixa(self, A):
        L = 2 / (A**2)
        return L

    # Numero Quantico da particula

    def numeroQuantico(self, L, k):
        n = round((k * L) / pi)
        return n

    # Propabilidade de encontarar a particula na posicao xp

    def probabilidadeParticula(self, L, n, xp):
        xp = xp * L
        p = (2 / L) * pow((sin((n * pi * xp) / L)), 2)
        return p

    def salvaGraficoFuncaoOnda(self, A, k, n, L):
        def f(x):
            return A * np.sin(k * x)

        fig, ax = plt.subplots(1)
        ax.set(
            title=("Gráfico da Função de Onda da Partícula no nível %d" % n),
            xlabel="x (Â)",
            ylabel=("Ψ%d" % n),
        )
        ax.axhline(y=0, color="k", linestyle="-", linewidth=0.5)
        x = np.linspace(0, L, 100)
        y = f(x)
        ax.plot(x, y)
        plt.savefig(("função_n%d.png" % n), dpi=200)
        plt.clf()

    def salvaGraficoProbabilidade(self, A, k, n, L):
        def f(x):
            return (A * np.sin(k * x))**2

        fig, ax = plt.subplots(1)
        ax.set(
            title=("Gráfico da Probabilidade da Partícula no nível %d" % n),
            xlabel="x (Â)",
            ylabel=("|Ψ%d|²" % n),
        )
        ax.axhline(y=0, color="k", linestyle="-", linewidth=0.5)
        x = np.linspace(0, L, 100)
        y = f(x)
        ax.plot(x, y)
        plt.savefig(("probabilidade_n%d.png" % n), dpi=200)
        plt.clf()


while True:
    print("Partícula em uma caixa\n")
    print("1. Caixa 1D – Determinação da função de onda quântica e outros parâmetros")
    print(
        "2. Caixa 1D – Cálculo dos parâmetros da caixa e partícula, dada a função de onda"
    )
    while True:
        opc_principal = input("Insira a opcao desejada: ")
        if opc_principal == "1" or opc_principal == "2":
            break
        else:
            print("Opcao invalida!")
    if opc_principal == "1":
        print("\n1- Confinando um elétron")
        print("2- Confinando um próton")
        while True:
            opc_principal1 = input("Insira a opcao desejada: ")
            if opc_principal1 == "1" or opc_principal1 == "2":
                break
            else:
                print("Opcao invalida!")
        if opc_principal1 == "1":
            L = float(input("\nLargura da caixa (L) em metros: "))
            ni = int(input("Número quântico inicial (ni): "))
            nf = int(input("Número quântico final (nf): "))
            a = float(input("Limite inferior da integral (a) em metros: "))
            b = float(input("Limite superior da integral (b) em metros: "))

            caixa = Caixa(L, ni, nf, a, b)
            A, k1 = caixa.funcaoQuanticaInicial(ni, L)
            A, k2 = caixa.funcaoQuanticaFinal(nf, L)
            Eij, Eiv, Efj, Efv = caixa.energiaQuanticaEletron(ni, nf, L)
            E, lambada, f = caixa.calculoFoton(Eiv, Efv)  # arrumar
            viE, vfE = caixa.velocidadeEletron(Eij, Efj)
            lambadai, lambdaf = caixa.comprimentoDeBroglieEletron(viE, vfE)
            integrali = caixa.probabilidadeIntegralNi(A, k1, a, b)
            integralf = caixa.probabilidadeIntegralNf(A, k2, a, b)
            caixa.salvaGraficoFuncaoOnda(A, k1, ni, L)
            caixa.salvaGraficoFuncaoOnda(A, k2, nf, L)
            caixa.salvaGraficoProbabilidade(A, k1, ni, L)
            caixa.salvaGraficoProbabilidade(A, k2, nf, L)
            print("\nFunção de onda inicial: {:.3e} * sin({:.3e} * x)".format(A, k1))
            print("Função de onda final: {:.3e} * sin({:.3e} * x)".format(A, k2))
            print("Energia inicial joules: {:.3e} j".format(Eij))
            print("Energia inicial eV: {:.3e} eV".format(Eiv))
            print("Energia final joules: {:.3e} j".format(Efj))
            print("Energia final eV: {:.3e} eV".format(Efv))
            print("Energia do fóton: {:.3e} eV".format(E))  # arrumar
            print("Frequência do fóton: {:.3e} Hz".format(f))  # arrumar
            print("Comprimento de onda do fóton: {:.3e} m".format(lambada))  # arrumar
            print("Velocidade inicial do elétron: {:.3e} m/s".format(viE))  # arrumar
            print("Velocidade final do elétron: {:.3e} m/s".format(vfE))
            print("Comprimento de onda inicial do elétron: {:.3e} m".format(lambadai))
            print("Comprimento de onda final do elétron: {:.3e} m".format(lambdaf))
            print(
                "Probabilidade de encontrar a partícula na posição x inicial: %.3f%%"
                % (integrali)
            )
            print(
                "Probabilidade de encontrar a partícula na posição x final: %.3f%%\n"
                % (integralf)
            )

        elif opc_principal1 == "2":
            L = float(input("Largura da caixa (L) em metros: "))
            ni = int(input("Número quântico inicial (ni): "))
            nf = int(input("Número quântico final (nf): "))
            a = float(input("Limite inferior da integral (a) em metros: "))
            b = float(input("Limite superior da integral (b) em metros: "))

            caixa = Caixa(L, ni, nf, a, b)
            A, k1 = caixa.funcaoQuanticaInicial(ni, L)
            A, k2 = caixa.funcaoQuanticaFinal(nf, L)
            Eij, Eiv, Efj, Efv = caixa.energiaQuanticaProton(ni, nf, L)
            E, lambada, f = caixa.calculoFoton(Eiv, Efv)  # arrumar
            viP, vfP = caixa.velocidadeProton(Eij, Efj)  # arrumar
            lambadai, lambdaf = caixa.comprimentoDeBroglieProton(viP, vfP)  # arrumar
            integrali = caixa.probabilidadeIntegralNi(A, k1, a, b)
            integralf = caixa.probabilidadeIntegralNf(A, k2, a, b)
            print("\nFunção de onda inicial: {:.3e} * sin({:.3e} * x)".format(A, k1))
            print("Função de onda final: {:.3e} * sin({:.3e} * x)".format(A, k2))
            print("Energia inicial joules: {:.3e} j".format(Eij))
            print("Energia inicial eV: {:.3e} eV".format(Eiv))
            print("Energia final joules: {:.3e} j".format(Efj))
            print("Energia final eV: {:.3e} eV".format(Efv))
            print("Energia do fóton: {:.3e} eV".format(E))  # arrumar
            print("Frequência do fóton: {:.3e} Hz".format(f))  # arrumar
            print("Comprimento de onda do fóton: {:.3e} m".format(lambada))  # arrumar
            print("Velocidade inicial do próton: {:.3e} m/s".format(viP))  # arrumar
            print("Velocidade final do próton: {:.3e} m/s".format(vfP))
            print("Comprimento de onda inicial do próton: {:.3e} m".format(lambadai))
            print("Comprimento de onda final do próton: {:.3e} m".format(lambdaf))
            print(
                "Probabilidade de encontrar a partícula na posição x inicial: %.3f%%".format(
                    integrali
                )
            )
            print(
                "Probabilidade de encontrar a partícula na posição x final: %.3f%%\n".format(
                    integralf
                )
            )

    elif opc_principal == "2":
        print("\n1- Confinando um elétron")
        print("2- Confinando um próton")
        while True:
            opc_principal2 = input("Insira a opcao desejada: ")
            if opc_principal2 == "1" or opc_principal2 == "2":
                break
            else:
                print("Opcao invalida!")
        if opc_principal2 == "1":
            A = float(input("\nAmplitude (A): "))
            k = float(input("K: "))
            xp = float(input("xp (m): "))

            caixa = Caixa(A, k, xp)
            L = caixa.larguraCaixa(A)
            n = caixa.numeroQuantico(k, L)
            p = caixa.probabilidadeParticula(L, n, xp)
            print("\nLargura da caixa: {:.3e} m".format(L))
            print("Número quântico: %.1f" % n)
            print("P(x) = {:.3e} .dx".format(p))

        elif opc_principal2 == "2":
            A = float(input("\nAmplitude (A): "))
            k = float(input("K: "))
            xp = float(input("xp (m): "))

            caixa = Caixa(A, k, xp)
            L = caixa.larguraCaixa(A)
            n = caixa.numeroQuantico(k, L)
            p = caixa.probabilidadeParticula(L, n, xp)
            print("\nLargura da caixa: {:.3e} m".format(L))
            print("Número quântico: %.1f" % n)
            print("P(x) = {:.3e} .dx".format(p))
