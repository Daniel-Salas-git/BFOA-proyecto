import math as mt
import random
from bacteria import bacteria

class chemiotaxis:
    def __init__(self):
        self.parcialNFE = 0  

    def compute_cell_interaction(self, bacteria, poblacion, d, w):
        total = 30
        for other in poblacion:
            diff = (bacteria.blosumScore - other.blosumScore) ** 2.0
            total += d * mt.exp(w * diff)
        return total

    def attract_repel(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
        attract = self.compute_cell_interaction(bacteria, poblacion, -d_attr, -w_attr)
        repel = self.compute_cell_interaction(bacteria, poblacion, h_rep, -w_rep)
        return attract + repel
    
    def chemio(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
        bacteria.interaction = self.attract_repel(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
        bacteria.fitness = bacteria.blosumScore + bacteria.interaction

    def doChemioTaxis(self, poblacion, d_attr, w_attr, h_rep, w_rep):
        self.parcialNFE = 0
        for bacteria in poblacion:
            self.chemio(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
            self.parcialNFE += bacteria.NFE
            bacteria.NFE = 0

    def eliminarClonar(self, path, poblacion):
        """Elimina el 50% de las bacterias con menos fitness"""
        poblacion.sort(key=lambda x: x.fitness)
        half_size = len(poblacion) // 2
        del poblacion[:half_size]
          
        clones = self.clonacion(path, poblacion)
        poblacion.extend(clones)

    
    def clonacion(self, path, poblacion):
        poblacionClones = []
        best = max(poblacion, key=lambda x: x.fitness)
        for bacteria in poblacion:
            newBacteria = bacteria.clonar(path)
            mutacion = int((best.fitness - bacteria.fitness) / 10)  # mutación en función de la diferencia de fitness
            newBacteria.tumboNado(mutacion)
            newBacteria.autoEvalua()
            poblacionClones.append(newBacteria)
        return poblacionClones

    def randomBacteria(self, path):
        bact = bacteria(path)
        bact.tumboNado(random.randint(1, 10))
        return bact 
   
    def insertRamdomBacterias(self, path, num, poblacion):
        for _ in range(num):
            poblacion.append(self.randomBacteria(path))
        # Elimina la bacteria con menos fitness
        poblacion.sort(key=lambda x: x.fitness)
        del poblacion[0]  # Solo elimina el de menor fitness después de agregar
