import PWM
import numpy as np
import os  


def generateOneRandomSequence(SL):
	sequence = [0]*SL
	for i in range(SL):
		sequence[i] = alphabet[np.random.randint(4)]
	return sequence

def generateRandomSequences(SC, SL):
	sequences = [0]*SC
	for i in range(SC):
		sequences[i] = generateOneRandomSequence(SL)
	return sequences

def generateOneMotif(ML, pwm):
	motif = [0]*ML
	for i in range(ML):
		motif[i] = np.random.choice(alphabet, p=pwm[i])
	return motif

def generateMotifs(ML, pwm):
	motifs = [0]*SC
	for i in range(SC):
		motifs[i] = generateOneMotif(ML, pwm)
	return motifs

def plantMotifs(sequences, motifs):
	positions = []
	for i in range(SC):
		position = np.random.randint(SL - ML)
		positions.append(position)
		sequences[i][position:position+ML] = motifs[i]
	return sequences, positions

path = os.getcwd()



alphabet = ['A', 'T', 'C', 'G']
#ICPC, ML, SL, SC
parameters = [
	[2, 	8, 	500, 	10],
	[1, 	8, 	500, 	10],
	[1.5, 	8, 	500, 	10],
	[2, 	6, 	500, 	10],
	[2, 	7, 	500, 	10],
	[2, 	8, 	500, 	5],
	[2, 	8, 	500, 	20]
]

for p in range(len(parameters)):
	[ICPC, ML, SL, SC] = parameters[p]
	for k in range(10):
		index = p*10 + k + 1
		datadir = path+'/data set ' + str(index) + '/'
		try:
			os.mkdir(datadir)
		except:
			print('data set folder already there')

		sequences = generateRandomSequences(SC, SL)
		pwm = PWM.getPWM(ICPC, ML)
		motifs = generateMotifs(ML, pwm)
		plantedSequnces, positions = plantMotifs(sequences, motifs)

		sequencesFile = open(datadir+'sequences.fa', 'w')
		for i in range(len(plantedSequnces)):
			sequencesFile.write('>seq' + str(i) + '\n')
			sequencesFile.write(''.join(plantedSequnces[i]) + '\n')
		sequencesFile.close()


		#position of sites
		sitesFile = open(datadir+'sites.txt', 'w')
		for i in range(len(positions)):
			sitesFile.write(str(positions[i]) + '\n')
		sitesFile.close()


		motifFile = open(datadir+'motif.txt', 'w')
		motifFile.write('>MOTIF' + str(index) + ' ' + str(ML) + '\n')
		for i in range(len(pwm)):
			motifFile.write(' '.join(str(j) for j in pwm[i]) + '\n')
		motifFile.write('<\n')
		motifFile.close()

		motifLengthFile = open(datadir+'motiflength.txt', 'w')
		motifLengthFile.write(str(ML))
		motifLengthFile.close()
