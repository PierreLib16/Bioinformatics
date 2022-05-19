'''
Created on 19 May 2022

@author: pierreliboureau
'''
from root.nested.RandomizedMotifSearch import Profile, Profile_Most, Score, Consensus


def GreedyMotifSearch(dna, k, t):
    best_motifs = []
    for strand in dna:
        best_motifs.append(strand[:k])
    for i in range(1, len(dna[0])-k):
        motifs = [dna[0][i:i+k]]
        for j in range(1, t):
            profile = Profile(motifs)
            motifs.append(Profile_Most(profile, k, dna[j]))
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs[:]
            
    return best_motifs, Score(best_motifs)

        
        
infile = '/Users/pierreliboureau/Downloads/dataset_160_9.txt'
with open(infile, 'r') as file:
    k, t = file.readline().split()
    k = int(k)
    t = int(t)
    dna = file.readline().split()
    
best_motifs, best_score = GreedyMotifSearch(dna, k, t)
print(' '.join(best_motifs))
print(best_score)