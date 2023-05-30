from src import readSequenceFile, getParameterDetails, pcaRegressionAlgorithm, motifsAlgorithm, processResults, writeFile,pca_training,reg_training,dataframe_woraround, cross_correlation
import datetime
# import src.pca
import matplotlib.pyplot as plt
# import seaborn as sns
import numpy as np
import gc

start = datetime.datetime.now()
filepath_tss = "F:\\Meethi Folder\\INTERNSHIPS\\IIT Delhi\modelknn\\train_data\\tss_20.txt"
try:
    f =open(filepath_tss)
except NameError:
    filepath_tss = str(input("Please enter the input sequence file."))

# filepath_no_tss = "F:\\Meethi Folder\\INTERNSHIPS\\IIT Delhi\\modelknn\\train_data\\cds.txt"
# df = dataframe_woraround.readingCSV()
# cross_correlation.corr_with_output(df)
# print(df.head)
# cross_correlation.corr_removal(df)
# try:
#     f = open(writeFilePath)
# except NameError:
#     writeFilePath = str(input("Please enter the output file path."))
#
# arr = filepath_tss.split("/")
# len = len(arr)
# writeFileName = arr[-1]
# path = writeFilePath+"/"+writeFileName

#POSITIVE DATASET
sequence_map_tss = readSequenceFile.readSequenceFile(filepath_tss)
parameter_map = {}
parameter_map = getParameterDetails.iterateSequences(sequence_map_tss)
# print(parameter_map_tss['normalized_params_map'].keys())
# print(sequence_map_tss.keys())
gc.collect()

#NEGATIVE DATASET
# sequence_map_nt = readSequenceFile.readSequenceFile(filepath_no_tss)
# parameter_map_nt = {}
# parameter_map_nt = getParameterDetails.iterateSequences(sequence_map_nt)
# print(parameter_map_nt['normalized_params_map'].keys())
# print(sequence_map_nt.keys())


#plotting to check if at 500, there is a dip CORRESPONDING TO EACH PARAMETRE(to indicate TSS)

# lineplt_ex = parameter_map_tss['normalized_params_map'][0]['a']
# a = [0 for i in range(976)]
# b = [0 for i in range(976)]
# c = [0 for i in range(976)]
#
#
# for i in range(19):
#     a=np.add(a,parameter_map_tss['normalized_params_map'][i]['a'])
#     b = np.add(b, parameter_map_tss['normalized_params_map'][i]['b'])
#     c = np.add(c, parameter_map_tss['normalized_params_map'][i]['g'])
#
# a = [i/19 for i in a]
# b = [i/19 for i in b]
# c = [i/19 for i in c]
# plt.plot(a,label = 'a',color = 'b')
# plt.plot(b,label = 'b',color='r')
# plt.plot(c,label = 'g',color='g')
# # plt.xlim(0,1000)
# # plt.ylim(0,1000)
#
# # plt.xticks(ran)
# plt.show()
#
# import sys
# sys.exit()


#NOW SEND THIS NORMALISED DATA TO TRAINING MODULE AND TEST WITH SAME ORGANISM SEQUENCE
pca_vectors_tss = pca_training.pca_training(parameter_map["normalized_params_map"])

# log_regression = reg_training.log_reg(pca_vectors_tss,pca_vector_nt)
#
# end = datetime.datetime.now()
#
# print(end-start)
# import sys
# sys.exit()

start = datetime.datetime.now()

map_pca = pcaRegressionAlgorithm.iterateSequences(parameter_map)
end = datetime.datetime.now()
# print(end-start)
# import sys
# sys.exit()
# print(parameter_map['combined_params_map'][0]['structuralIncreasing_params'])
# import sys
# sys.exit()
start = datetime.datetime.now()
map_motif = motifsAlgorithm.iterateSequences(parameter_map)
end = datetime.datetime.now()
print(end-start)

# print("fine till here")

final_map = {}
for seq in map_pca.keys():
    map_pca_1 = map_pca[seq]
    map_motif_1 = map_motif[seq]
    final_map[seq] = {}
    for start in map_pca_1.keys():
        final_map[seq][start] = []
        map_pca_2 = map_pca_1[start]
        map_motif_2 = map_motif_1[start]
        final_map[seq][start].append(map_pca_2[0] and map_motif_2[0])
        final_map[seq][start].append(map_pca_2[1] and map_motif_2[1])

for seq in final_map.keys():
    for start in (final_map[seq]).keys():
        # print(final_map[1][569])
        if (final_map[seq][start][0] == 1):
            final_map[seq][start] = 1
        else:
            final_map[seq][start] = 0


writeFile.writeResults("F:\\Meethi Folder\\INTERNSHIPS\\IIT Delhi\modelknn\\train_data\\",final_map)

# print (final_map[2])
import sys
sys.exit()

result_map = processResults.predictSequencewiseTss(final_map)
# print (result_map.keys())
pos_map = processResults.getTssPositions(result_map)
print(pos_map)
import sys
sys.exit()
writeFile.writeResults(path,final_map)
