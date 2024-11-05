from bacteria import bacteria as Bacteria
from chemiotaxis import chemiotaxis
import numpy
import matplotlib.pyplot as plt  # Para graficar

poblacion = []
path = "multiFasta.fasta"
numeroDeBacterias = 5
numRandomBacteria = 1
iteraciones = 30
tumbo = 1  # número de gaps a insertar
nado = 3
chemio = chemiotaxis()
bestBacteria = None  # Inicializa como None
tempBacteria = Bacteria(path)  # bacteria temporal para validaciones
original = Bacteria(path)  # bacteria original sin gaps
globalNFE = 0  # número de evaluaciones de la función objetivo

dAttr = 0.1
wAttr = 0.2
hRep = dAttr
wRep = 10

def clonaBest(bestBacteria, best):
    bestBacteria.matrix.seqs = numpy.array(best.matrix.seqs)
    bestBacteria.blosumScore = best.blosumScore
    bestBacteria.fitness = best.fitness
    bestBacteria.interaction = best.interaction

def validaSecuencias(path, bestBacteria):
    tempBacteria.matrix.seqs = numpy.array(bestBacteria.matrix.seqs)
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-", "")
    
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return

for i in range(numeroDeBacterias):  # población inicial
    poblacion.append(Bacteria(path))

for _ in range(iteraciones):  # número de iteraciones  
    for bacteria_instance in poblacion:
        bacteria_instance.tumboNado(tumbo)
#--------------------        
        bacteria_instance.mejorarFitness()  # Mejora el fitness tras la evaluación
#--------------------
    chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)
    globalNFE += chemio.parcialNFE 
    best = max(poblacion, key=lambda x: x.fitness)
    
    if (bestBacteria is None) or (best.fitness > bestBacteria.fitness):
        bestBacteria = Bacteria(path)  # Crea una nueva instancia
        clonaBest(bestBacteria, best)

    print("\n--- Resultados de la Iteración ---")
    print(f"Interacción: {bestBacteria.interaction:.2f}")
    print(f"Fitness: {bestBacteria.fitness:.2f}")
    print(f"NFE: {globalNFE}")
    print(f"Tamaño de la población: {len(poblacion)}")
    print("-------------------------------\n")
    
    chemio.eliminarClonar(path, poblacion)
    chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)  # inserta bacterias aleatorias

bestBacteria.showGenome()
validaSecuencias(path, bestBacteria)
#-----------------------------------------------
def experimentar(self, path, dAttr, wAttr, hRep, wRep, num_iteraciones=30):
    fitness_results = []
    nfe_results = []

    for _ in range(num_iteraciones):
        self.run(path, dAttr, wAttr, hRep, wRep)
        fitness_results.append(self.veryBest.fitness)
        nfe_results.append(self.globalNFE)

    # Llama al método para graficar los resultados
    self.graficar_resultados(fitness_results, nfe_results)
    
def graficar_resultados(self, fitness_results, nfe_results):
    plt.figure(figsize=(10, 5))

    # Gráfica de Fitness
    plt.subplot(1, 2, 1)
    plt.plot(fitness_results, marker='o', color='b')
    plt.title("Evolución del Fitness")
    plt.xlabel("Iteración")
    plt.ylabel("Fitness")

    # Gráfica de NFE
    plt.subplot(1, 2, 2)
    plt.plot(nfe_results, marker='o', color='r')
    plt.title("Evolución del NFE")
    plt.xlabel("Iteración")
    plt.ylabel("NFE")

    plt.tight_layout()
    plt.show()