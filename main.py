# Simulador: Partícula em uma caixa
# em orientacao a objetos
# 1. Caixa 1D – Determinação da função de onda quântica e outros parâmetros
# 2. Caixa 1D – Cálculo dos parâmetros da caixa e partícula, dada a função de onda

# Importando as bibliotecas

from math import *
from scipy import integrate
from matplotlib import pyplot as plt
import numpy as np
from tkinter import *

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
    
    def plotaGraficos(self, A, k1, k2, n1, n2, L):

        def f1(x):
            return A * np.sin(k1 * x)
        
        def f2(x):
            return A * np.sin(k2 * x)
        
        def g1(x):
            return (A * np.sin(k1 * x))**2
        
        def g2(x):
            return (A * np.sin(k2 * x))**2
        
        x = np.linspace(0, L, 100)

        try:
            plt.close()
        except:
            pass
        
        fig, ax = plt.subplots(2,2, figsize=(12, 8))
        fig.tight_layout(pad=5.0)

        ax[0,0].set(
            title=("Gráfico da Função de Onda da Partícula no nível %d" % n1),
            xlabel="x (Â)",
            ylabel=("Ψ%d" % n1),
        )
        ax[0,0].axhline(y=0, color="k", linestyle="-", linewidth=0.5)

        yf1 = f1(x)
        ax[0,0].plot(x, yf1)

        ax[0,1].set(
            title=("Gráfico da Função de Onda da Partícula no nível %d" % n2),
            xlabel="x (Â)",
            ylabel=("Ψ%d" % n2),
        )
        ax[0,1].axhline(y=0, color="k", linestyle="-", linewidth=0.5)

        yf2 = f2(x)
        ax[0,1].plot(x, yf2)

        ax[1,0].set(
            title=("Gráfico da Probabilidade da Partícula no nível %d" % n1),
            xlabel="x (Â)",
            ylabel=("|Ψ%d|²" % n1),
        )
        ax[1,0].axhline(y=0, color="k", linestyle="-", linewidth=0.5)

        yg1 = g1(x)
        ax[1,0].plot(x, yg1)

        ax[1,1].set(
            title=("Gráfico da Probabilidade da Partícula no nível %d" % n2),
            xlabel="x (Â)",
            ylabel=("|Ψ%d|²" % n2),
        )
        ax[1,1].axhline(y=0, color="k", linestyle="-", linewidth=0.5)

        yg2 = g2(x)
        ax[1,1].plot(x, yg2)

class GUI:
    def __init__(self):

        # cria a janela
        self.janela = Tk()
        self.janela.title("Simulador: Partícula em uma caixa")
        self.janela.geometry("800x800")
        self.janela.configure(background="white")
        self.janela.resizable(width=False, height=False)

        # label título
        self.titulo = Label(self.janela, text="Simulador: Partícula em uma caixa", bg="white", fg="black", font="Calibri 25 bold", anchor=CENTER)
        self.titulo.place(relx=0.5, rely=0.1, relwidth=0.9, relheight=0.1, anchor=CENTER)

        # mensagem de introdução
        self.mensagem = Label(self.janela, text="Aqui o Ruan irá escrever um resumo sobre a matéria, o que é partícula numa caixa,\n etc etc etc.", bg="white", fg="black", font="Calibri 14", anchor=CENTER)
        self.mensagem.place(relx=0.5, rely=0.2, relwidth=0.9, relheight=0.1, anchor=CENTER)

        # botão da primeira opção de simulação
        self.botao1 = Button(self.janela, text="Caixa 1D – Determinação da função de onda quântica e outros parâmetros", bg="SeaGreen1", fg="black", font="Calibri 14", anchor=CENTER, command=self.opcao1, borderwidth=2, relief="solid", activebackground="SeaGreen3")
        self.botao1.place(relx=0.5, rely=0.65, relwidth=0.8, relheight=0.1, anchor=CENTER)

        # botão da segunda opção de simulação
        self.botao2 = Button(self.janela, text="Caixa 1D – Cálculo dos parâmetros da caixa e partícula, dada a função de onda", bg="SkyBlue1", fg="black", font="Calibri 14", anchor=CENTER, command=self.opcao2, borderwidth=2, relief="solid", activebackground="SkyBlue3")
        self.botao2.place(relx=0.5, rely=0.80, relwidth=0.8, relheight=0.1, anchor=CENTER)

        # botão de sair
        self.botaoSair = Button(self.janela, text="Sair", bg="white", fg="black", font="Calibri 12", anchor=CENTER, command=self.janela.quit, borderwidth=2, relief="solid", activebackground="red3")
        self.botaoSair.place(relx=0.925, rely=0.95, relwidth=0.1, relheight=0.05, anchor=CENTER)

    def opcao1(self):
        # cria a janela da opcao1
        self.janelaOpcao1 = Toplevel()
        self.janelaOpcao1.title("Caixa 1D – Determinação da função de onda quântica e outros parâmetros")
        self.janelaOpcao1.geometry("800x800")
        self.janelaOpcao1.configure(background="white")
        self.janelaOpcao1.resizable(width=False, height=False)

        # label título
        self.titulo = Label(self.janelaOpcao1, text="Caixa 1D – Determinação da função de onda quântica e outros parâmetros", bg="white", fg="black", font="Calibri 14 bold", anchor=CENTER)
        self.titulo.place(relx=0.5, rely=0.05, relwidth=0.9, relheight=0.1, anchor=CENTER)

        # radiobutton de seleção para elétron ou próton
        self.varOpcao = StringVar()
        self.labelParticula = Label(self.janelaOpcao1, text="Selecione a partícula:", bg="white", fg="black", font="Calibri 14", anchor=W)
        self.labelParticula.place(relx=0.02, rely=0.12, relwidth=0.3, relheight=0.05, anchor=W)
        self.radioEletron = Radiobutton(self.janelaOpcao1, text="Elétron", bg="white", fg="black", font="Calibri 14", anchor=W, variable=self.varOpcao, value="eletron")
        self.radioEletron.place(relx=0.02, rely=0.18, relwidth=0.1, relheight=0.03, anchor=W)
        self.radioProton = Radiobutton(self.janelaOpcao1, text="Próton", bg="white", fg="black", font="Calibri 14", anchor=W, variable=self.varOpcao, value="proton")
        self.radioProton.place(relx=0.02, rely=0.23, relwidth=0.1, relheight=0.03, anchor=W)
        self.varOpcao.set("eletron")

        # entrada das informações
        self.labelL = Label(self.janelaOpcao1, text="Largura da caixa (L): ", bg="white", fg="black", font="Calibri 14", anchor=W)
        self.labelL.place(relx=0.02, rely=0.32, relwidth=0.4, relheight=0.1, anchor=W)
        self.entradaL = Entry(self.janelaOpcao1, bg="white", fg="black", font="Calibri 14", justify=CENTER, borderwidth=2, relief="solid")
        self.entradaL.place(relx=0.33, rely=0.32, relwidth=0.12, relheight=0.05, anchor=CENTER)
        self.labelMetro = Label(self.janelaOpcao1, text="m", bg="white", fg="black", font="Calibri 14", anchor=CENTER)
        self.labelMetro.place(relx=0.41, rely=0.32, relwidth=0.02, relheight=0.05, anchor=CENTER)

        self.labelNi = Label(self.janelaOpcao1, text="n inicial da partícula (ni): ", bg="white", fg="black", font="Calibri 14", anchor=W)
        self.labelNi.place(relx=0.02, rely=0.4, relwidth=0.4, relheight=0.08, anchor=W)
        self.entradaNi = Entry(self.janelaOpcao1, bg="white", fg="black", font="Calibri 14", justify=CENTER, borderwidth=2, relief="solid")
        self.entradaNi.place(relx=0.33, rely=0.4, relwidth=0.12, relheight=0.05, anchor=CENTER)

        self.labelNf = Label(self.janelaOpcao1, text="n final da partícula (nf): ", bg="white", fg="black", font="Calibri 14", anchor=W)
        self.labelNf.place(relx=0.02, rely=0.48, relwidth=0.4, relheight=0.08, anchor=W)
        self.entradaNf = Entry(self.janelaOpcao1, bg="white", fg="black", font="Calibri 14", justify=CENTER, borderwidth=2, relief="solid")
        self.entradaNf.place(relx=0.33, rely=0.48, relwidth=0.12, relheight=0.05, anchor=CENTER)

        self.labelProbabilidade = Label(self.janelaOpcao1, text="Dados para probabilidade:\nP(a <= x <= b) ", bg="white", fg="black", font="Calibri 14", anchor=W)
        self.labelProbabilidade.place(relx=0.02, rely=0.6, relwidth=0.4, relheight=0.1, anchor=W)

        self.labelA = Label(self.janelaOpcao1, text="a: ", bg="white", fg="black", font="Calibri 14", anchor=W)
        self.labelA.place(relx=0.02, rely=0.7, relwidth=0.4, relheight=0.05, anchor=W)
        self.entradaA = Entry(self.janelaOpcao1, bg="white", fg="black", font="Calibri 14", justify=CENTER, borderwidth=2, relief="solid")
        self.entradaA.place(relx=0.12, rely=0.7, relwidth=0.12, relheight=0.05, anchor=CENTER)
        self.labelMetroA = Label(self.janelaOpcao1, text="m", bg="white", fg="black", font="Calibri 14", anchor=CENTER)
        self.labelMetroA.place(relx=0.20, rely=0.7, relwidth=0.02, relheight=0.05, anchor=CENTER)

        self.labelB = Label(self.janelaOpcao1, text="b: ", bg="white", fg="black", font="Calibri 14", anchor=W)
        self.labelB.place(relx=0.02, rely=0.78, relwidth=0.4, relheight=0.05, anchor=W)
        self.entradaB = Entry(self.janelaOpcao1, bg="white", fg="black", font="Calibri 14", justify=CENTER, borderwidth=2, relief="solid")
        self.entradaB.place(relx=0.12, rely=0.78, relwidth=0.12, relheight=0.05, anchor=CENTER)
        self.labelMetroB = Label(self.janelaOpcao1, text="m", bg="white", fg="black", font="Calibri 14", anchor=CENTER)
        self.labelMetroB.place(relx=0.20, rely=0.78, relwidth=0.02, relheight=0.05, anchor=CENTER)

        # botao calcular
        self.botaoCalcular = Button(self.janelaOpcao1, text="Calcular", bg="SeaGreen1", fg="black", font="Calibri 14 bold", anchor=CENTER, command=self.calcularFuncaoOndaOutrosParametros, borderwidth=2, relief="solid", activebackground="SeaGreen3")
        self.botaoCalcular.place(relx=0.2, rely=0.9, relwidth=0.17, relheight=0.08, anchor=CENTER)

        # exibição dos resultados
        self.labelFuncaoOndaInicial = Label(self.janelaOpcao1, text="Função de onda inicial: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelFuncaoOndaInicial.place(relx=0.75, rely=0.1, relwidth=0.4, relheight=0.025, anchor=CENTER)
        self.funcaoOndaInicial = Text(self.janelaOpcao1, height=1, borderwidth=0, font="Calibri 14 bold")
        self.funcaoOndaInicial.tag_configure("center", justify="center")
        self.funcaoOndaInicial.insert(1.0, "ψi(x) = A · sin(k · x)")
        self.funcaoOndaInicial.configure(state="disabled")
        self.funcaoOndaInicial.tag_add("center", "1.0", "end")
        self.funcaoOndaInicial.place(relx=0.75, rely=0.13, relwidth=0.4, relheight=0.04, anchor=CENTER)

        self.labelFuncaoOndaFinal = Label(self.janelaOpcao1, text="Função de onda final: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelFuncaoOndaFinal.place(relx=0.75, rely=0.16, relwidth=0.4, relheight=0.025, anchor=CENTER)
        self.funcaoOndaFinal = Text(self.janelaOpcao1, height=1, borderwidth=0, font="Calibri 14 bold")
        self.funcaoOndaFinal.tag_configure("center", justify="center")
        self.funcaoOndaFinal.insert(1.0, "ψf(x) = A· sin(k · x)")
        self.funcaoOndaFinal.configure(state="disabled")
        self.funcaoOndaFinal.tag_add("center", "1.0", "end")
        self.funcaoOndaFinal.place(relx=0.75, rely=0.19, relwidth=0.4, relheight=0.04, anchor=CENTER)

        self.labelEnergiaInicial = Label(self.janelaOpcao1, text="Energia inicial: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelEnergiaInicial.place(relx=0.75, rely=0.22, relwidth=0.4, relheight=0.025, anchor=CENTER)
        self.energiaInicial = Text(self.janelaOpcao1, height=1, borderwidth=0, font="Calibri 14 bold")
        self.energiaInicial.tag_configure("center", justify="center")
        self.energiaInicial.insert(1.0, "Ei = 0.0 J = 0.0 eV")
        self.energiaInicial.configure(state="disabled")
        self.energiaInicial.tag_add("center", "1.0", "end")
        self.energiaInicial.place(relx=0.75, rely=0.25, relwidth=0.4, relheight=0.04, anchor=CENTER)

        self.labelEnergiaFinal = Label(self.janelaOpcao1, text="Energia final: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelEnergiaFinal.place(relx=0.75, rely=0.28, relwidth=0.4, relheight=0.025, anchor=CENTER)
        self.energiaFinal = Text(self.janelaOpcao1, height=1, borderwidth=0, font="Calibri 14 bold")
        self.energiaFinal.tag_configure("center", justify="center")
        self.energiaFinal.insert(1.0, "Ef = 0.0 J = 0.0 eV")
        self.energiaFinal.configure(state="disabled")
        self.energiaFinal.tag_add("center", "1.0", "end")
        self.energiaFinal.place(relx=0.75, rely=0.31, relwidth=0.4, relheight=0.04, anchor=CENTER)

        self.labelEnergiaFoton = Label(self.janelaOpcao1, text="Energia do fóton: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelEnergiaFoton.place(relx=0.75, rely=0.34, relwidth=0.4, relheight=0.025, anchor=CENTER)
        self.energiaFoton = Text(self.janelaOpcao1, height=1, borderwidth=0, font="Calibri 14 bold")
        self.energiaFoton.tag_configure("center", justify="center")
        self.energiaFoton.insert(1.0, "Eγ = 0.0 J = 0.0 eV")
        self.energiaFoton.configure(state="disabled")
        self.energiaFoton.tag_add("center", "1.0", "end")
        self.energiaFoton.place(relx=0.75, rely=0.37, relwidth=0.4, relheight=0.04, anchor=CENTER)

        self.labelComprimentoOndaFoton = Label(self.janelaOpcao1, text="Comprimento de onda do fóton: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelComprimentoOndaFoton.place(relx=0.75, rely=0.40, relwidth=0.4, relheight=0.025, anchor=CENTER)
        self.comprimentoOndaFoton = Text(self.janelaOpcao1, height=1, borderwidth=0, font="Calibri 14 bold")
        self.comprimentoOndaFoton.tag_configure("center", justify="center")
        self.comprimentoOndaFoton.insert(1.0, "λγ = 0.0 m")
        self.comprimentoOndaFoton.configure(state="disabled")
        self.comprimentoOndaFoton.tag_add("center", "1.0", "end")
        self.comprimentoOndaFoton.place(relx=0.75, rely=0.43, relwidth=0.4, relheight=0.04, anchor=CENTER)

        self.labelFrequenciaFoton = Label(self.janelaOpcao1, text="Frequência do fóton: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelFrequenciaFoton.place(relx=0.75, rely=0.46, relwidth=0.4, relheight=0.025, anchor=CENTER)
        self.frequenciaFoton = Text(self.janelaOpcao1, height=1, borderwidth=0, font="Calibri 14 bold")
        self.frequenciaFoton.tag_configure("center", justify="center")
        self.frequenciaFoton.insert(1.0, "νγ = 0.0 Hz")
        self.frequenciaFoton.configure(state="disabled")
        self.frequenciaFoton.tag_add("center", "1.0", "end")
        self.frequenciaFoton.place(relx=0.75, rely=0.49, relwidth=0.4, relheight=0.04, anchor=CENTER)

        self.labelVelocidadeInicial = Label(self.janelaOpcao1, text="Velocidade inicial da partícula: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelVelocidadeInicial.place(relx=0.75, rely=0.52, relwidth=0.4, relheight=0.025, anchor=CENTER)
        self.velocidadeInicial = Text(self.janelaOpcao1, height=1, borderwidth=0, font="Calibri 14 bold")
        self.velocidadeInicial.tag_configure("center", justify="center")
        self.velocidadeInicial.insert(1.0, "Vi = 0.0 m/s")
        self.velocidadeInicial.configure(state="disabled")
        self.velocidadeInicial.tag_add("center", "1.0", "end")
        self.velocidadeInicial.place(relx=0.75, rely=0.55, relwidth=0.4, relheight=0.04, anchor=CENTER)

        self.labelVelocidadeFinal = Label(self.janelaOpcao1, text="Velocidade final da partícula: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelVelocidadeFinal.place(relx=0.75, rely=0.58, relwidth=0.4, relheight=0.025, anchor=CENTER)
        self.velocidadeFinal = Text(self.janelaOpcao1, height=1, borderwidth=0, font="Calibri 14 bold")
        self.velocidadeFinal.tag_configure("center", justify="center")
        self.velocidadeFinal.insert(1.0, "Vf = 0.0 m/s")
        self.velocidadeFinal.configure(state="disabled")
        self.velocidadeFinal.tag_add("center", "1.0", "end")
        self.velocidadeFinal.place(relx=0.75, rely=0.61, relwidth=0.4, relheight=0.04, anchor=CENTER)

        self.labelComprimentoOndaInicialParticula = Label(self.janelaOpcao1, text="Comprimento de onda da partícula no nível inicial: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelComprimentoOndaInicialParticula.place(relx=0.75, rely=0.64, relwidth=0.5, relheight=0.025, anchor=CENTER)
        self.comprimentoOndaInicialParticula = Text(self.janelaOpcao1, height=1, borderwidth=0, font="Calibri 14 bold")
        self.comprimentoOndaInicialParticula.tag_configure("center", justify="center")
        self.comprimentoOndaInicialParticula.insert(1.0, "λi = 0.0 m")
        self.comprimentoOndaInicialParticula.configure(state="disabled")
        self.comprimentoOndaInicialParticula.tag_add("center", "1.0", "end")
        self.comprimentoOndaInicialParticula.place(relx=0.75, rely=0.67, relwidth=0.4, relheight=0.04, anchor=CENTER)

        self.labelComprimentoOndaFinalParticula = Label(self.janelaOpcao1, text="Comprimento de onda da partícula no nível final: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelComprimentoOndaFinalParticula.place(relx=0.75, rely=0.70, relwidth=0.4, relheight=0.025, anchor=CENTER)
        self.comprimentoOndaFinalParticula = Text(self.janelaOpcao1, height=1, borderwidth=0, font="Calibri 14 bold")
        self.comprimentoOndaFinalParticula.tag_configure("center", justify="center")
        self.comprimentoOndaFinalParticula.insert(1.0, "λf = 0.0 m")
        self.comprimentoOndaFinalParticula.configure(state="disabled")
        self.comprimentoOndaFinalParticula.tag_add("center", "1.0", "end")
        self.comprimentoOndaFinalParticula.place(relx=0.75, rely=0.73, relwidth=0.4, relheight=0.04, anchor=CENTER)

        self.labelProbabilidadeParticulaInicial = Label(self.janelaOpcao1, text="P(a <= x <= b) no nível inicial: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelProbabilidadeParticulaInicial.place(relx=0.75, rely=0.76, relwidth=0.5, relheight=0.025, anchor=CENTER)
        self.probabilidadeParticulaInicial = Text(self.janelaOpcao1, height=1, borderwidth=0, font="Calibri 14 bold")
        self.probabilidadeParticulaInicial.tag_configure("center", justify="center")
        self.probabilidadeParticulaInicial.insert(1.0, "Pi = 0.0 %")
        self.probabilidadeParticulaInicial.configure(state="disabled")
        self.probabilidadeParticulaInicial.tag_add("center", "1.0", "end")
        self.probabilidadeParticulaInicial.place(relx=0.75, rely=0.79, relwidth=0.4, relheight=0.04, anchor=CENTER)

        self.labelProbabilidadeParticulaFinal = Label(self.janelaOpcao1, text="P(a <= x <= b) no nível final: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelProbabilidadeParticulaFinal.place(relx=0.75, rely=0.82, relwidth=0.5, relheight=0.025, anchor=CENTER)
        self.probabilidadeParticulaFinal = Text(self.janelaOpcao1, height=1, borderwidth=0, font="Calibri 14 bold")
        self.probabilidadeParticulaFinal.tag_configure("center", justify="center")
        self.probabilidadeParticulaFinal.insert(1.0, "Pf = 0.0 %")
        self.probabilidadeParticulaFinal.configure(state="disabled")
        self.probabilidadeParticulaFinal.tag_add("center", "1.0", "end")
        self.probabilidadeParticulaFinal.place(relx=0.75, rely=0.85, relwidth=0.4, relheight=0.04, anchor=CENTER)

        # botao que abre janela dos gráficos
        self.botaoGraficos = Button(self.janelaOpcao1, text="Gráficos", bg="SeaGreen1", fg="black", font="Calibri 14 bold", command=self.abrirGraficos, borderwidth=2, relief="solid", activebackground="SeaGreen3")
        self.botaoGraficos.place(relx=0.75, rely=0.93, relwidth=0.2, relheight=0.05, anchor=CENTER)

    def calcularFuncaoOndaOutrosParametros(self):
        try:
            L = float(self.entradaL.get())
            ni = int(self.entradaNi.get())
            nf = int(self.entradaNf.get())
            a = float(self.entradaA.get())
            b = float(self.entradaB.get())
            caixa = Caixa(L, ni, nf, a, b)
        except:
            print("Insira informações válidas!")

        try:
            A, k1 = caixa.funcaoQuanticaInicial(ni, L)
            A, k2 = caixa.funcaoQuanticaFinal(nf, L)
            integrali = caixa.probabilidadeIntegralNi(A, k1, a, b)
            integralf = caixa.probabilidadeIntegralNf(A, k2, a, b)
        except:
            print("Erro ao calcular!")

        if self.varOpcao.get() == "eletron":
            try:
                Eij, Eiv, Efj, Efv = caixa.energiaQuanticaEletron(ni, nf, L)
                E, lambada, f = caixa.calculoFoton(Eiv, Efv)
                viE, vfE = caixa.velocidadeEletron(Eij, Efj)
                lambadai, lambdaf = caixa.comprimentoDeBroglieEletron(viE, vfE)

                self.energiaInicial.configure(state="normal")
                self.energiaInicial.delete(1.0, END)
                self.energiaInicial.insert(1.0, "E{} = {:.3e} j = {:.3e} ev".format(ni, Eij, Eiv))
                self.energiaInicial.tag_add("center", "1.0", "end")
                self.energiaInicial.configure(state="disabled")

                self.energiaFinal.configure(state="normal")
                self.energiaFinal.delete(1.0, END)
                self.energiaFinal.insert(1.0, "E{} = {:.3e} j = {:.3e} ev".format(nf, Efj, Efv))
                self.energiaFinal.tag_add("center", "1.0", "end")
                self.energiaFinal.configure(state="disabled")

                self.energiaFoton.configure(state="normal")
                self.energiaFoton.delete(1.0, END)
                self.energiaFoton.insert(1.0, "Ef = {:.3e} ev".format(E))
                self.energiaFoton.tag_add("center", "1.0", "end")
                self.energiaFoton.configure(state="disabled")

                self.frequenciaFoton.configure(state="normal")
                self.frequenciaFoton.delete(1.0, END)
                self.frequenciaFoton.insert(1.0, "f = {:.3e} Hz".format(f))
                self.frequenciaFoton.tag_add("center", "1.0", "end")
                self.frequenciaFoton.configure(state="disabled")

                self.comprimentoOndaFoton.configure(state="normal")
                self.comprimentoOndaFoton.delete(1.0, END)
                self.comprimentoOndaFoton.insert(1.0, "λf = {:.3e} m".format(lambada))
                self.comprimentoOndaFoton.tag_add("center", "1.0", "end")
                self.comprimentoOndaFoton.configure(state="disabled")

                self.velocidadeInicial.configure(state="normal")
                self.velocidadeInicial.delete(1.0, END)
                self.velocidadeInicial.insert(1.0, "V{} = {:.3e} m/s".format(ni, viE))
                self.velocidadeInicial.tag_add("center", "1.0", "end")
                self.velocidadeInicial.configure(state="disabled")

                self.velocidadeFinal.configure(state="normal")
                self.velocidadeFinal.delete(1.0, END)
                self.velocidadeFinal.insert(1.0, "V{} = {:.3e} m/s".format(nf, vfE))
                self.velocidadeFinal.tag_add("center", "1.0", "end")
                self.velocidadeFinal.configure(state="disabled")

                self.comprimentoOndaInicialParticula.configure(state="normal")
                self.comprimentoOndaInicialParticula.delete(1.0, END)
                self.comprimentoOndaInicialParticula.insert(1.0, "λ{} = {:.3e} m".format(ni, lambadai))
                self.comprimentoOndaInicialParticula.tag_add("center", "1.0", "end")
                self.comprimentoOndaInicialParticula.configure(state="disabled")

                self.comprimentoOndaFinalParticula.configure(state="normal")
                self.comprimentoOndaFinalParticula.delete(1.0, END)
                self.comprimentoOndaFinalParticula.insert(1.0, "λ{} = {:.3e} m".format(nf, lambdaf))
                self.comprimentoOndaFinalParticula.tag_add("center", "1.0", "end")
                self.comprimentoOndaFinalParticula.configure(state="disabled")
            except:
                print("Erro ao calcular!")
        else:
            try:
                Eij, Eiv, Efj, Efv = caixa.energiaQuanticaProton(ni, nf, L)
                E, lambada, f = caixa.calculoFoton(Eiv, Efv)
                viP, vfP = caixa.velocidadeProton(Eij, Efj)
                lambadai, lambdaf = caixa.comprimentoDeBroglieProton(viP, vfP)

                self.energiaInicial.configure(state="normal")
                self.energiaInicial.delete(1.0, END)
                self.energiaInicial.insert(1.0, "E{} = {:.3e} j = {:.3e} ev".format(ni, Eij, Eiv))
                self.energiaInicial.tag_add("center", "1.0", "end")
                self.energiaInicial.configure(state="disabled")

                self.energiaFinal.configure(state="normal")
                self.energiaFinal.delete(1.0, END)
                self.energiaFinal.insert(1.0, "E{} = {:.3e} j = {:.3e} ev".format(nf, Efj, Efv))
                self.energiaFinal.tag_add("center", "1.0", "end")
                self.energiaFinal.configure(state="disabled")

                self.energiaFoton.configure(state="normal")
                self.energiaFoton.delete(1.0, END)
                self.energiaFoton.insert(1.0, "Ef = {:.3e} ev".format(E))
                self.energiaFoton.tag_add("center", "1.0", "end")
                self.energiaFoton.configure(state="disabled")

                self.frequenciaFoton.configure(state="normal")
                self.frequenciaFoton.delete(1.0, END)
                self.frequenciaFoton.insert(1.0, "f = {:.3e} Hz".format(f))
                self.frequenciaFoton.tag_add("center", "1.0", "end")
                self.frequenciaFoton.configure(state="disabled")

                self.comprimentoOndaFoton.configure(state="normal")
                self.comprimentoOndaFoton.delete(1.0, END)
                self.comprimentoOndaFoton.insert(1.0, "λf = {:.3e} m".format(lambada))
                self.comprimentoOndaFoton.tag_add("center", "1.0", "end")
                self.comprimentoOndaFoton.configure(state="disabled")

                self.velocidadeInicial.configure(state="normal")
                self.velocidadeInicial.delete(1.0, END)
                self.velocidadeInicial.insert(1.0, "V{} = {:.3e} m/s".format(ni, viP))
                self.velocidadeInicial.tag_add("center", "1.0", "end")
                self.velocidadeInicial.configure(state="disabled")

                self.velocidadeFinal.configure(state="normal")
                self.velocidadeFinal.delete(1.0, END)
                self.velocidadeFinal.insert(1.0, "V{} = {:.3e} m/s".format(nf, vfP))
                self.velocidadeFinal.tag_add("center", "1.0", "end")
                self.velocidadeFinal.configure(state="disabled")

                self.comprimentoOndaInicialParticula.configure(state="normal")
                self.comprimentoOndaInicialParticula.delete(1.0, END)
                self.comprimentoOndaInicialParticula.insert(1.0, "λ{} = {:.3e} m".format(ni, lambadai))
                self.comprimentoOndaInicialParticula.tag_add("center", "1.0", "end")
                self.comprimentoOndaInicialParticula.configure(state="disabled")

                self.comprimentoOndaFinalParticula.configure(state="normal")
                self.comprimentoOndaFinalParticula.delete(1.0, END)
                self.comprimentoOndaFinalParticula.insert(1.0, "λ{} = {:.3e} m".format(nf, lambdaf))
                self.comprimentoOndaFinalParticula.tag_add("center", "1.0", "end")
                self.comprimentoOndaFinalParticula.configure(state="disabled")
            except:
                print("Erro ao calcular!")

        try:
            self.funcaoOndaInicial.configure(state="normal")
            self.funcaoOndaInicial.delete(1.0, END)
            self.funcaoOndaInicial.insert(1.0, "ψ{}(x) = {:.3e} · sin({:.3e} · x)".format(ni, A, k1))
            self.funcaoOndaInicial.tag_add("center", "1.0", "end")
            self.funcaoOndaInicial.configure(state="disabled")
            
            self.funcaoOndaFinal.configure(state="normal")
            self.funcaoOndaFinal.delete(1.0, END)
            self.funcaoOndaFinal.insert(1.0, "ψ{}(x) = {:.3e} · sin({:.3e} · x)".format(nf, A, k2))
            self.funcaoOndaFinal.tag_add("center", "1.0", "end")
            self.funcaoOndaFinal.configure(state="disabled")

            self.probabilidadeParticulaInicial.configure(state="normal")
            self.probabilidadeParticulaInicial.delete(1.0, END)
            self.probabilidadeParticulaInicial.insert(1.0, "P%d = %.3f %%" % (ni, integrali))
            self.probabilidadeParticulaInicial.tag_add("center", "1.0", "end")
            self.probabilidadeParticulaInicial.configure(state="disabled")

            self.probabilidadeParticulaFinal.configure(state="normal")
            self.probabilidadeParticulaFinal.delete(1.0, END)
            self.probabilidadeParticulaFinal.insert(1.0, "P%d = %.3f %%" % (nf, integralf))
            self.probabilidadeParticulaFinal.tag_add("center", "1.0", "end")
            self.probabilidadeParticulaFinal.configure(state="disabled")
        except:
            print("Erro ao calcular!")

        caixa.plotaGraficos(A, k1, k2, ni, nf, L)
 
    def abrirGraficos(self):
        plt.show()

    def opcao2(self):
        # cria a janela da opcao2
        self.janelaOpcao2 = Toplevel(self.janela)
        self.janelaOpcao2.title("Caixa 1D – Cálculo dos parâmetros da caixa e partícula, dada a função de onda")
        self.janelaOpcao2.geometry("800x400")
        self.janelaOpcao2.configure(background="white")
        self.janelaOpcao2.resizable(width=False, height=False)

        self.titulo2 = Label(self.janelaOpcao2, text="Caixa 1D – Cálculo dos parâmetros da caixa e partícula, dada a função de onda", bg="white", fg="black", font="Calibri 16 bold", anchor=CENTER)
        self.titulo2.place(relx=0.5, rely=0.05, relwidth=0.9, relheight=0.1, anchor=CENTER)

        self.labelA = Label(self.janelaOpcao2, text="A: ", bg="white", fg="black", font="Calibri 14", anchor=W)
        self.labelA.place(relx=0.02, rely=0.3, relwidth=0.4, relheight=0.1, anchor=W)
        self.entradaA = Entry(self.janelaOpcao2, bg="white", fg="black", font="Calibri 14", justify=CENTER, borderwidth=2, relief="solid")
        self.entradaA.place(relx=0.12, rely=0.3, relwidth=0.12, relheight=0.1, anchor=CENTER)
        self.labelMetroAA = Label(self.janelaOpcao2, text="m", bg="white", fg="black", font="Calibri 14", anchor=CENTER)
        self.labelMetroAA.place(relx=0.2, rely=0.3, relwidth=0.02, relheight=0.1, anchor=CENTER)

        self.labelK = Label(self.janelaOpcao2, text="K: ", bg="white", fg="black", font="Calibri 14", anchor=W)
        self.labelK.place(relx=0.02, rely=0.5, relwidth=0.4, relheight=0.1, anchor=W)
        self.entradaK = Entry(self.janelaOpcao2, bg="white", fg="black", font="Calibri 14", justify=CENTER, borderwidth=2, relief="solid")
        self.entradaK.place(relx=0.12, rely=0.5, relwidth=0.12, relheight=0.1, anchor=CENTER)

        self.labelXP = Label(self.janelaOpcao2, text="xp: ", bg="white", fg="black", font="Calibri 14", anchor=W)
        self.labelXP.place(relx=0.02, rely=0.7, relwidth=0.4, relheight=0.1, anchor=W)
        self.entradaXP = Entry(self.janelaOpcao2, bg="white", fg="black", font="Calibri 14", justify=CENTER, borderwidth=2, relief="solid")
        self.entradaXP.place(relx=0.12, rely=0.7, relwidth=0.12, relheight=0.1, anchor=CENTER)

        self.labelLarguraCaixa = Label(self.janelaOpcao2, text="Largura da caixa: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelLarguraCaixa.place(relx=0.75, rely=0.25, relwidth=0.4, relheight=0.05, anchor=CENTER)
        self.larguraCaixa = Text(self.janelaOpcao2, height=1, borderwidth=0, font="Calibri 14 bold")
        self.larguraCaixa.tag_configure("center", justify="center")
        self.larguraCaixa.insert(1.0, "0.0 m")
        self.larguraCaixa.configure(state="disabled")
        self.larguraCaixa.tag_add("center", "1.0", "end")
        self.larguraCaixa.place(relx=0.75, rely=0.30, relwidth=0.4, relheight=0.065, anchor=CENTER)

        self.labelNivelQuanticoParticula = Label(self.janelaOpcao2, text="Nível quântico da partícula: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelNivelQuanticoParticula.place(relx=0.75, rely=0.45, relwidth=0.4, relheight=0.05, anchor=CENTER)
        self.nivelQuanticoParticula = Text(self.janelaOpcao2, height=1, borderwidth=0, font="Calibri 14 bold")
        self.nivelQuanticoParticula.tag_configure("center", justify="center")
        self.nivelQuanticoParticula.insert(1.0, "0")
        self.nivelQuanticoParticula.configure(state="disabled")
        self.nivelQuanticoParticula.tag_add("center", "1.0", "end")
        self.nivelQuanticoParticula.place(relx=0.75, rely=0.5, relwidth=0.4, relheight=0.065, anchor=CENTER)

        self.labelProbabilidadePosicaoXP = Label(self.janelaOpcao2, text="Probabilidade de encontrar a partícula em xp: ", bg="white", fg="black", font="Calibri 12", anchor=CENTER)
        self.labelProbabilidadePosicaoXP.place(relx=0.75, rely=0.65, relwidth=0.4, relheight=0.05, anchor=CENTER)
        self.probabilidadePosicaoXP = Text(self.janelaOpcao2, height=1, borderwidth=0, font="Calibri 14 bold")
        self.probabilidadePosicaoXP.tag_configure("center", justify="center")
        self.probabilidadePosicaoXP.insert(1.0, "0.0")
        self.probabilidadePosicaoXP.configure(state="disabled")
        self.probabilidadePosicaoXP.tag_add("center", "1.0", "end")
        self.probabilidadePosicaoXP.place(relx=0.75, rely=0.7, relwidth=0.4, relheight=0.065, anchor=CENTER)

        self.botaoCalcular2 = Button(self.janelaOpcao2, text="Calcular", bg="SkyBlue1", fg="black", font="Calibri 14 bold", anchor=CENTER, command=self.calcularParametrosCaixa, borderwidth=2, relief="solid", activebackground="SkyBlue3")
        self.botaoCalcular2.place(relx=0.2, rely=0.9, relwidth=0.17, relheight=0.08, anchor=CENTER)

    def calcularParametrosCaixa(self):
        try:
            A = float(self.entradaA.get())
            k = float(self.entradaK.get())
            xp = float(self.entradaXP.get())
            caixa = Caixa(A, k, xp)
        except:
            print("Insira informações válidas!")

        try:
            L = caixa.larguraCaixa(A)
            n = caixa.numeroQuantico(k, L)
            p = caixa.probabilidadeParticula(L, n, xp)
        except:
            print("Erro ao calcular!")

        try:
            self.larguraCaixa.configure(state="normal")
            self.larguraCaixa.delete(1.0, END)
            self.larguraCaixa.insert(1.0, "L = {:.3e} m".format(L))
            self.larguraCaixa.tag_add("center", "1.0", "end")
            self.larguraCaixa.configure(state="disabled")

            self.nivelQuanticoParticula.configure(state="normal")
            self.nivelQuanticoParticula.delete(1.0, END)
            self.nivelQuanticoParticula.insert(1.0, ("n = %d" % (n)))
            self.nivelQuanticoParticula.tag_add("center", "1.0", "end")
            self.nivelQuanticoParticula.configure(state="disabled")

            self.probabilidadePosicaoXP.configure(state="normal")
            self.probabilidadePosicaoXP.delete(1.0, END)
            self.probabilidadePosicaoXP.insert(1.0, "P(x) = {:.3e} · dx".format(p))
            self.probabilidadePosicaoXP.tag_add("center", "1.0", "end")
            self.probabilidadePosicaoXP.configure(state="disabled")
        except:
            print("Erro ao mostrar resultados!")

if __name__ == "__main__":
    gui = GUI()
    gui.janela.mainloop()