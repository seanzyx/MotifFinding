import random
import numpy as np
import math

##get input
def get_input(datasetNumber):
    fa_in = open("data set "+str(datasetNumber)+"/sequences.fa","r")
    sequences = []# define a list to store the sequences
    num_lines = 0
    for line in fa_in.readlines():# read FASTA format line by line
        line = line.rstrip()# trim the line
        num_lines += 1
        if num_lines % 2 == 0:
            sequences.append(line)
    fa_in.close()
    num_seq = int(num_lines/2)
    length_in = open("data set "+str(datasetNumber)+"motiflength.txt","r")
    motif_len = int(length_in.readline().strip())
    length_in.close()
    return sequences, motif_len

def get_pwm(updated_motifs, motif_len, num_seq):
    counts = [[0 for i in range(motif_len)] for j in range(4)]
    for i in range(motif_len):
        cur = [k[i] for k in updated_motifs]
        counts[0][i] = cur.count('A')
        counts[1][i] = cur.count('T')
        counts[2][i] = cur.count('C')
        counts[3][i] = cur.count('G')
    counts = np.array(counts)
    probs = counts/(num_seq-1)
    ##print(pro_matrix)
    return probs

def generate_prob(probs,temp):
    dic = {'A':0,'T':1,'C':2,'G':3}
    p = 1
    for i in range(len(temp)):
        p = p * probs[dic[temp[i]]][i]
    return p

def motifs_score(updated_motifs,background_p, motif_len, num_seq):
    #motif_count_matrix
    counts = [[0 for i in range(motif_len)] for j in range(4)]
    for i in range(motif_len):
        cur = [k[i] for k in updated_motifs]
        counts[0][i] = cur.count('A')
        counts[1][i] = cur.count('T')
        counts[2][i] = cur.count('C')
        counts[3][i] = cur.count('G')
    counts = np.array(counts)
    #motif_probability_matrix
    probs = counts/num_seq
    #print(count_matrix)
    score = 0
    for i in range(motif_len):
        for j in range(4):
            if counts[j][i] == 0:
                score += 0
            else:
                score += counts[j][i]*math.log2(probs[j][i]/background_p[j])
    return score





background = [0.25, 0.25, 0.25, 0.25]

def findBestSites(datasetNumber, drops):
    ##get input
    sequences, motif_len = get_input(datasetNumber)

    currentBestSites = []
    currentBestScore = float('-inf')
    for d in range(drops):
        currentSites, currentScore = oneDrop(10, sequences, motif_len)
        if currentScore > currentBestScore:
            currentBestScore = currentScore
            currentBestSites = currentSites
    return currentBestSites, currentBestScore


def oneDrop(iters, sequences, motif_len):
    ## generate random locations in the sequences, and get the motifs
    num_seq = len(sequences)
    sites = [random.randint(0, (len(sequences[0])-motif_len)) for i in range(num_seq)]
    random_motifs = []
    for i in range(num_seq):
        random_motifs.append(list(sequences[i][sites[i]:sites[i]+motif_len]))
    score = float('-inf')

    for it in range(iters):
        ## remove one of the sequences randomly
        hide_index = random.randint(0,num_seq-1)
        updated_motifs = random_motifs.copy()
        cur_sites = sites.copy()
        updated_motifs.pop(hide_index)

        ## get PWM from the rest of the motifs, generate probability distribution
        probs = get_pwm(updated_motifs, motif_len, num_seq)
        prob = []
        hidden_seq = sequences[hide_index]
        for i in  range(len(sequences[0])-motif_len+1):
            temp_motif = hidden_seq[i:i+motif_len]
            prob.append(generate_prob(probs,temp_motif))
        ##use distribution to generate motif site

        # updated_site = np.random.choice(len(sequences[0])-motif_len+1,1,prob)[0] ##do we use random

        updated_site = np.argmax(prob)

        cur_sites[hide_index] = updated_site
        updated_motifs.insert(hide_index, list(sequences[hide_index][cur_sites[hide_index]:cur_sites[hide_index]+motif_len]))

        ##calculate F for updated motifs
        cur_score = motifs_score(updated_motifs, background, motif_len, num_seq)
        # print(cur_score, score)
        ##update if score is higher
        if cur_score > score:
            sites = cur_sites
            random_motifs = updated_motifs
            score = cur_score
    return sites, score

for i in range(70):
    sites, score = findBestSites(i+1, 100)
    print(i+1, score)
    print(sites)
