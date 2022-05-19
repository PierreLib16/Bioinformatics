'''
Created on 19 May 2022

@author: pierreliboureau
'''

import random

def GibbsSampler(dna, k, t, N):
    best_motifs = []
    for strand in dna:
        n = random.randint(0, len(strand) - k)
        best_motifs.append(strand[n:n + k])
    for j in range(N):
        i = random.randint(0, t-1)
        r_motifs = best_motifs[:i]+best_motifs[i+1:]
        profile = Profile(r_motifs)
        new_motif = WeightedSelect(profile, dna[i], k)
        motifs = best_motifs[:i] + [new_motif] + best_motifs[i+1:]
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs[:]
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

def WeightedSelect(profile, strand, k):
    all_probs = {}
    for i in range(len(strand)-k):
        motif = strand[i:i+k]
        if motif in all_probs.keys():
            continue
        prob = 1
        for j in range(len(motif)):
            if motif[j] == 'A':
                prob = prob * profile[j][0]
            elif motif[j] == 'C':
                prob = prob * profile[j][1]
            elif motif[j] == 'G':
                prob = prob * profile[j][2]
            elif motif[j] == 'T':
                prob = prob * profile[j][3]
        all_probs[motif] = prob
    sum_probs = sum(all_probs.values())
    distrib = [0]
    list_values = list(all_probs.values())
    list_motifs = list(all_probs.keys())
    for i in range(len(list_values)):
        distrib.append(list_values[i]/sum_probs+distrib[-1])
    n = random.random()
    for i in range(len(distrib)-1):
        if n > distrib[i] and n < distrib[i+1]:
            return list_motifs[i]

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



def MultiRun(dna, k, t, N, iterations):
    i = 0
    min_score = 1000000
    while i < iterations:
        motifs, score = GibbsSampler(dna, k, t, N)
        if score < min_score:
            min_score = score
            best_motifs = motifs
        i += 1
    return best_motifs, min_score

infile = ('/Users/pierreliboureau/Downloads/dataset_163_4.txt')
with open(infile, 'r') as file:
    k, t, N = file.readline().split()
    k = int(k)
    t = int(t)
    N = int(N)
    dna = file.readline().split()
    
best_motifs, best_score = MultiRun(dna, k, t, N, 20)
print(best_score)
print(' '.join(best_motifs))