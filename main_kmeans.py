from src import readSequenceFile, getParameterDetails, pcaRegressionAlgorithm, motifsAlgorithm, processResults, writeFile,pca_training,reg_training,dataframe_woraround, cross_correlation
import datetime
import pandas as pd

start = datetime.datetime.now()
filepath_tss = "F:\\Meethi Folder\\INTERNSHIPS\\IIT Delhi\modelknn\\train_data\\tss_20.txt"
try:
    f =open(filepath_tss)
except NameError:
    filepath_tss = str(input("Please enter the input sequence file."))


sequence_map_tss = readSequenceFile.readSequenceFile(filepath_tss)
parameter_map = {}
parameter_map = getParameterDetails.iterateSequences(sequence_map_tss)

parameter_map = pd.DataFrame(parameter_map['normalized_params_map'])
parameter_map.to_csv("parameter_map.csv")

import sys
sys.exit()