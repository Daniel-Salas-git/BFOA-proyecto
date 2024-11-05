import blosum as bl

class evaluadorBlosum:
    
    def __init__(self):
        self.matrix = bl.BLOSUM(62)
        
    def showMatrix(self):
        print(self.matrix)
        
    def getScore(self, A, B):
        if A in self.matrix and B in self.matrix[A]:
            score = self.matrix[A][B]
            return score
        else:
            raise ValueError(f"One or both residues '{A}' or '{B}' are not in the BLOSUM matrix.")
