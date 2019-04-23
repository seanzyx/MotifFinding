import random
import numpy as np
import math

##get input
def get_input(datasetNumber):
    fa_in = open("data set "+str(datasetNumber)+"/sequences.fa","r")
    fa_Seq = []# define a list to store the sequences
    fa_Num = 0
    for line in fa_in.readlines():# read FASTA format line by line
        line = line.rstrip()# trim the line
        fa_Num = fa_Num + 1
        if fa_Num % 2 == 0:
            fa_Seq.append(line)
    fa_in.close()
    num_seq = int(fa_Num/2)
    length_in = open("data set 1/motiflength.txt","r")
    motif_len = int(length_in.readline().strip())
    length_in.close()
    return fa_Seq, motif_len

def get_pwm(updated_motifs, motif_len, num_seq):
    counts = [[0 for i in range(motif_len)] for j in range(4)]
    for i in range(motif_len):
        for j in range(num_seq-1):
            if updated_motifs[j][i] == 'A':
                counts[0][i] += 1
            elif updated_motifs[j][i] == 'T':
                counts[1][i] += 1
            elif updated_motifs[j][i] == 'C':
                counts[2][i] += 1
            elif updated_motifs[j][i] == 'G':
                counts[3][i] += 1
    counts = np.array(counts) 
    probs = counts/(num_seq-1)
    ##print(pro_matrix)
    return probs

def generate_prob(probs,temp): #k_mer 'CAAATCCC'
    '''
    Calculate the probability of current k_mer
    '''
    dic = {'A':0,'T':1,'C':2,'G':3}
    p = 1
    for i in range(len(temp)):
        p = p * probs[dic[temp[i]]][i]
    return p

def motifs_score(temp_k_mers,background_p, motif_len, num_seq):
    dic = {'A':0,'T':1,'C':2,'G':3}
    sigma = ['A','T','C','G']
    
    #motif_count_matrix
    count_matrix = [[0 for i in range(motif_len)] for j in range(4)]
    for i in range(motif_len):
        for j in range(num_seq):
            if temp_k_mers[j][i] == 'A':
                count_matrix[0][i] += 1
            elif temp_k_mers[j][i] == 'T':
                count_matrix[1][i] += 1
            elif temp_k_mers[j][i] == 'C':
                count_matrix[2][i] += 1
            elif temp_k_mers[j][i] == 'G':
                count_matrix[3][i] += 1
    count_matrix = np.array(count_matrix)
    #motif_probability_matrix
    pro_matrix = count_matrix/num_seq
    #print(count_matrix)
    F = 0
    for i in range(len(temp_k_mers[0])):
        for j in range(4):
            if count_matrix[j][i] == 0:
                F += 0
            else:
                F += count_matrix[j][i]*math.log2(pro_matrix[j][i]/background_p[j])
    return F





background = [0.25, 0.25, 0.25, 0.25]

def findBestSites(datasetNumber, drops):
    ##get input
    fa_Seq, motif_len = get_input(datasetNumber)

    currentBestSites = []
    currentBestScore = float('-inf')
    for d in range(drops):
        currentSites, currentScore = oneDrop(10, fa_Seq, motif_len)
        if currentScore > currentBestScore:
            currentBestScore = currentScore
            currentBestSites = currentSites
    return currentBestSites, currentBestScore


def oneDrop(iters, fa_Seq, motif_len):  
    ## generate random locations in the sequences, and get the motifs
    num_seq = len(fa_Seq)
    sites = [random.randint(0, (len(fa_Seq[0])-motif_len)) for i in range(num_seq)] 
    random_motifs = [] 
    for i in range(num_seq):
        random_motifs.append(list(fa_Seq[i][sites[i]:sites[i]+motif_len]))
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
        hidden_seq = fa_Seq[hide_index]
        for i in  range(len(fa_Seq[0])-motif_len+1):
            temp_motif = hidden_seq[i:i+motif_len]
            prob.append(generate_prob(probs,temp_motif))
        ##use distribution to generate motif site

        # updated_site = np.random.choice(len(fa_Seq[0])-motif_len+1,1,prob)[0] ##do we use random

        updated_site = np.argmax(prob)

        cur_sites[hide_index] = updated_site
        updated_motifs.insert(hide_index, list(fa_Seq[hide_index][cur_sites[hide_index]:cur_sites[hide_index]+motif_len]))
        
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
