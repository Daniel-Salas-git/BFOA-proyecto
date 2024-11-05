from fastaReader import fastaReader
import random as rd
import numpy as np
import copy
from evaluadorBlosum import evaluadorBlosum

class bacteria:
    
    def __init__(self, path):
        self.matrix = fastaReader(path)
        self.blosumScore = 0
        self.fitness = 0
        self.interaction = 0
        self.NFE = 0
        
    def showGenome(self):
        for seq in self.matrix.seqs:
            print(seq)

    def clonar(self, path):
        newBacteria = bacteria(path)
        newBacteria.matrix.seqs = np.array(copy.deepcopy(self.matrix.seqs))
        return newBacteria

    def tumboNado(self, numGaps):
        self.cuadra()
        matrixCopy = copy.deepcopy(self.matrix.seqs).tolist()
        gapRandomNumber = rd.randint(0, numGaps)  # número de gaps a insertar
        for i in range(gapRandomNumber):  # ciclo de gaps 
            seqnum = rd.randint(0, len(matrixCopy) - 1)  # selecciono secuencia
            pos = rd.randint(0, len(matrixCopy[0]) - 1)  # ajuste aquí
            part1 = matrixCopy[seqnum][:pos]
            part2 = matrixCopy[seqnum][pos:]
            temp = "-".join([part1, part2])  # inserto gap
            matrixCopy[seqnum] = temp
        self.matrix.seqs = np.array(matrixCopy)  # convierto a numpy array
        self.cuadra()
        self.limpiaColumnas()

    def cuadra(self):
        """Rellena con gaps las secuencias más cortas"""
        seq = self.matrix.seqs
        maxLen = len(max(seq, key=len))
        for i in range(len(seq)):
            if len(seq[i]) < maxLen:
                seq[i] = seq[i] + "-" * (maxLen - len(seq[i]))
        self.matrix.seqs = np.array(seq)

    def gapColumn(self, col):
        for i in range(len(self.matrix.seqs)):
            if self.matrix.seqs[i][col] != "-":
                return False
        return True

    def limpiaColumnas(self):
        i = 0
        while i < len(self.matrix.seqs[0]):
            if self.gapColumn(i):
                self.deleteColumn(i)  # ajuste en el nombre
            else:
                i += 1
        
    def deleteColumn(self, pos):
        for i in range(len(self.matrix.seqs)):
            self.matrix.seqs[i] = self.matrix.seqs[i][:pos] + self.matrix.seqs[i][pos + 1:]

    def getColumn(self, col):
        return [self.matrix.seqs[i][col] for i in range(len(self.matrix.seqs))]

    def autoEvalua(self):   
        evaluador = evaluadorBlosum()
        score = 0
        for i in range(len(self.matrix.seqs[0])):
            column = self.getColumn(i)
            gapCount = column.count("-")
            column = [x for x in column if x != "-"]
            pares = self.obtener_pares_unicos(column)
            for par in pares:
                score += evaluador.getScore(par[0], par[1])
            score -= gapCount * 2
        self.blosumScore = score
        self.NFE += 1

    def obtener_pares_unicos(self, columna):
        return list({tuple(sorted([columna[i], columna[j]])) for i in range(len(columna)) for j in range(i + 1, len(columna))})
    
    #------------------------------------------------------------------
    
    def mejorarFitness(self):
        """Método para aumentar el fitness mediante variaciones estratégicas."""
        # Multiplicar el fitness por un factor que depende de su puntuación actual
        if self.blosumScore > 0:
            factor = 1 + (self.blosumScore / 100)  # Factor de mejora basado en la puntuación actual
            self.fitness *= factor  # Aumentar el fitness
        else:
            self.fitness += 1  # Aumentar en un valor fijo si la puntuación es negativa

        # Opcionalmente, aplicar un pequeño ajuste aleatorio
        ajuste_aleatorio = rd.uniform(-0.5, 0.5)
        self.fitness += ajuste_aleatorio  # Ajuste aleatorio
