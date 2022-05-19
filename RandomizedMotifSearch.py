'''
Created on 18 May 2022

@author: pierreliboureau
'''
from random import randint

def RandomizedMotifSearch(dna, k, t):
    motifs = []
    for strand in dna:
        l = len(strand)-k
        rand = randint(0, l)
        motifs.append(strand[rand:rand+k])
    best_motifs = motifs[:]
    i = 0
    while i == 0:
        profile = Profile(motifs)
        motifs = Motifs(profile, dna)
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs [:]
        else: 
            return best_motifs, Score(best_motifs)

def Profile(motifs):
    profiles = []
    for i in range(len(motifs[0])):
        counts = [1,1,1,1]
        tot_count = 4
        for j in range(len(motifs)):
            motif = list(motifs[j])
            if motif[i] == 'A': counts[0] +=1
            elif motif[i] == 'C': counts[1] +=1
            elif motif[i] == 'G': counts[2] +=1
            elif motif[i] == 'T': counts[3] +=1
            tot_count += 1
        profile = [counts[0]/tot_count, counts[1]/tot_count, counts[2]/tot_count, counts[3]/tot_count]
        profiles.append(profile)
    return profiles

def Motifs(profile, dna):
    motifs = []
    k = len(profile)
    for strand in dna:
        motifs.append(Profile_Most(profile, k, strand))
    return motifs


def Profile_Most(profile, k, strand):
    max_prob = 0
    for i in range(len(strand)-k+1):
        pattern = strand[i:i+k]
        prob_pattern = 1
        for j in range(len(pattern)):
            if pattern[j] == 'A':
                prob_pattern = prob_pattern * profile[j][0]
            elif pattern[j] == 'C':
                prob_pattern = prob_pattern * profile[j][1] 
            elif pattern[j] == 'G':
                prob_pattern = prob_pattern * profile[j][2] 
            elif pattern[j] == 'T':
                prob_pattern = prob_pattern * profile[j][3]
        if prob_pattern > max_prob:
            max_prob = prob_pattern
            best_pattern = pattern[:]
    return best_pattern


def Score(motifs):
    consensus = Consensus(motifs)
    score = 0
    for motif in motifs:
        score += HammingDistance(motif, consensus)
    return score

def Consensus(motifs):
    consensus = ''
    for i in range(len(motifs[0])):
        counts = [0,0,0,0]
        for j in range(len(motifs)):
            if motifs[j][i] == 'A': counts[0] += 1
            elif motifs[j][i] == 'C': counts[1] += 1
            elif motifs[j][i] == 'G': counts[2] += 1
            elif motifs[j][i] == 'T': counts[3] += 1
        max_value = max(counts)
        max_index = counts.index(max_value)
        if max_index == 0: consensus += 'A'
        elif max_index == 1: consensus += 'C'
        elif max_index == 2: consensus += 'G'
        elif max_index == 3: consensus += 'T'
    return consensus

def HammingDistance(gen1, gen2):
    counter = 0
    for i in range(0,len(gen1)):
            if gen1[i] != gen2[i]:
                    counter +=1
    return counter
             
        
def multiRandom(dna, k, t, iterations):
    i = 0
    min_score = 1000000
    while i < iterations:
        motifs, score = RandomizedMotifSearch(dna, k, t)
        if score < min_score:
            min_score = score
            best_motifs = motifs
        i += 1
    return best_motifs, min_score
                
    
#infile = ('/Users/pierreliboureau/Downloads/dataset_161_5.txt')  
#with open(infile, 'r') as file:
#    k, t = file.readline().split()
#    k = int(k)
#    t = int(t)
#    dna = file.readline().split()
    

#best_motif, best_score = multiRandom(dna, k, t, 5000)
#print(' '.join(best_motif))
#print(best_score)
#with open('output.txt', 'w') as f:
#    print('\r\n'.join(best_motif), file=f, flush=True)