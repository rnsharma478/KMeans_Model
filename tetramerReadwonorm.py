from src import readSequenceFile, getParameterDetails, pcaRegressionAlgorithm, motifsAlgorithm, processResults, writeFile,pca_training,reg_training,dataframe_woraround, cross_correlation
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

tetramer = 'C:/Users/User/SEProm-master/train_data/TSS_Seq/Synechocystis sp. PCC6803'
paramVals = 'C:/Users/User/SEProm-master/train_data/new_tetramer.xlsx'
tetramer_cds = 'C:/Users/User/SEProm-master/train_data/CDS_Seq/Synechocystis_sp_PCC6803_CDS'

strInc = []
strDec = []

def dataCleaning(df):
    # df = df.drop(['Tetramer'], axis =1)
    duplicates = df.duplicated(subset = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r'], keep = False)
    df2 = df[~duplicates]
    df1 = df[duplicates]
    df1 = (df1.groupby(df1.columns.tolist())
           .apply(lambda x: tuple(x.index))
           .reset_index(name='Tetramer'))
    df1=df1.set_index('Tetramer')
    df = pd.concat([df1,df2])
    df = df.T
    return df

paramValues = pd.read_excel(paramVals, sheet_name='Sheet2', index_col=0)
paramValues = dataCleaning(paramValues)
#print(paramValues)





def calculateParameters(sequence_map_per_seq,paramValues):
    param_map = {'a':[],'b':[],	'c':[],	'd':[],	'e':[],	'f':[], 'g':[],'h':[],'i':[],'j':[],'k':[],'l':[],'m':[],'n':[],'o':[],'p':[],'q':[],'r':[]}
    # shift = slide= rise = tilt = roll = twist =0
    #sequence_map_per_seq = sequence_map_per_seq[2:-1]
    no_of_bases = len(sequence_map_per_seq)
    list_motifs = []
    for m in range(no_of_bases-3):
        list_motifs.append(sequence_map_per_seq[m:m+4])

# to continue from here

    for motif in list_motifs:
        if motif in paramValues.columns:
            param_map['a'].append(paramValues[motif]['a'])
            param_map['b'].append(paramValues[motif]['b'])
            param_map['c'].append(paramValues[motif]['c'])
            param_map['d'].append(paramValues[motif]['d'])
            param_map['e'].append(paramValues[motif]['e'])
            param_map['f'].append(paramValues[motif]['f'])
            param_map['g'].append(paramValues[motif]['g'])
            param_map['h'].append(paramValues[motif]['h'])
            param_map['i'].append(paramValues[motif]['i'])
            param_map['j'].append(paramValues[motif]['j'])
            param_map['k'].append(paramValues[motif]['k'])
            param_map['l'].append(paramValues[motif]['l'])
            param_map['m'].append(paramValues[motif]['m'])
            param_map['n'].append(paramValues[motif]['n'])
            param_map['o'].append(paramValues[motif]['o'])
            param_map['p'].append(paramValues[motif]['p'])
            param_map['q'].append(paramValues[motif]['q'])
            param_map['r'].append(paramValues[motif]['r'])
        else:
            for j in range(len(paramValues.columns)):
                if motif in paramValues.columns[j]:
                    # print(("found inside"))
                    param_map['a'].append(paramValues[paramValues.columns[j]]['a'])
                    param_map['b'].append(paramValues[paramValues.columns[j]]['b'])
                    param_map['c'].append(paramValues[paramValues.columns[j]]['c'])
                    param_map['d'].append(paramValues[paramValues.columns[j]]['d'])
                    param_map['e'].append(paramValues[paramValues.columns[j]]['e'])
                    param_map['f'].append(paramValues[paramValues.columns[j]]['f'])
                    param_map['g'].append(paramValues[paramValues.columns[j]]['g'])
                    param_map['h'].append(paramValues[paramValues.columns[j]]['h'])
                    param_map['i'].append(paramValues[paramValues.columns[j]]['i'])
                    param_map['j'].append(paramValues[paramValues.columns[j]]['j'])
                    param_map['k'].append(paramValues[paramValues.columns[j]]['k'])
                    param_map['l'].append(paramValues[paramValues.columns[j]]['l'])
                    param_map['m'].append(paramValues[paramValues.columns[j]]['m'])
                    param_map['n'].append(paramValues[paramValues.columns[j]]['n'])
                    param_map['o'].append(paramValues[paramValues.columns[j]]['o'])
                    param_map['p'].append(paramValues[paramValues.columns[j]]['p'])
                    param_map['q'].append(paramValues[paramValues.columns[j]]['q'])
                    param_map['r'].append(paramValues[paramValues.columns[j]]['r'])
    return calculateMovingAverages(param_map)

def calculateMovingAverages(param_map):
    moving_win_size = 25
    moving_param_map = {}
    for k, v in param_map.items():
        arr = v
        moving_param_map[k] = []
        for i in range(0, len(arr) - moving_win_size + 1):
            sum = 0
            for j in range(i, i + moving_win_size):
                sum += arr[j]
            avg = sum / moving_win_size
            moving_param_map[k].append(avg)
    return moving_param_map

'''
def normalizeMovingAverages(moving_param_map):
    normalized_map = {}
    for k in moving_param_map.keys():
        arr = moving_param_map[k]
        # maxArr = max(arr)
        # minArr = min(arr)
        rang = max(arr)-min(arr)
        normalized_map[k] = []
        for i in arr:
            norm_val = (i-min(arr))/rang
            normalized_map[k].append(norm_val)
    return normalized_map
'''


def iterateSequences(sequence_map):
    parameters = {
        'values_map': {},
        'moving_averages_map': {},
        'normalized_params_map': {},
        'combined_params_map': {
            'structuralIncreasing_params': {},
            'structuralDecreasing_params': {}
        }
    }
    for key in sequence_map.keys():
        parameters['normalized_params_map'][key] = calculateParameters(sequence_map[key],paramValues)
        # print(parameters['normalized_params_map'][key]['l'])

    return parameters


#main
sequence_map_tss = readSequenceFile.readSequenceFile(tetramer)
sequence_map_cds = readSequenceFile.readSequenceFile(tetramer_cds)
tetra_seq_map = iterateSequences(sequence_map_tss)
#tetra_seq_map.to_csv('tetra_seq_map.csv')
tetra_cds_seq_map = iterateSequences(sequence_map_cds)
#tetra_cds_seq_map.to_csv('tetra_cds_seq_map.csv')



#plotting to check if at 500, there is a dip CORRESPONDING TO EACH PARAMETRE(to indicate TSS)
a = [0 for i in range(974)]
b = [0 for i in range(974)]
c = [0 for i in range(974)]
d = [0 for i in range(974)]
e = [0 for i in range(974)]
f = [0 for i in range(974)]
g = [0 for i in range(974)]
h = [0 for i in range(974)]
i = [0 for z in range(974)]
j = [0 for i in range(974)]
k = [0 for i in range(974)]
l = [0 for i in range(974)]
m = [0 for i in range(974)]
n = [0 for i in range(974)]
o = [0 for i in range(974)]
p = [0 for i in range(974)]
q = [0 for i in range(974)]
r = [0 for i in range(974)]


a_cds = [0 for i in range(974)]
b_cds = [0 for i in range(974)]
c_cds = [0 for i in range(974)]
d_cds = [0 for i in range(974)]
e_cds = [0 for i in range(974)]
f_cds = [0 for i in range(974)]
g_cds = [0 for i in range(974)]
h_cds = [0 for i in range(974)]
i_cds = [0 for z in range(974)]
j_cds = [0 for i in range(974)]
k_cds = [0 for i in range(974)]
l_cds = [0 for i in range(974)]
m_cds = [0 for i in range(974)]
n_cds = [0 for i in range(974)]
o_cds = [0 for i in range(974)]
p_cds = [0 for i in range(974)]
q_cds = [0 for i in range(974)]
r_cds = [0 for i in range(974)]


for z in range(len(sequence_map_tss)):
    a = np.add(a, tetra_seq_map['normalized_params_map'][z]['a'])

    b = np.add(b, tetra_seq_map['normalized_params_map'][z]['b'])

    c = np.add(c, tetra_seq_map['normalized_params_map'][z]['c'])

    d = np.add(d, tetra_seq_map['normalized_params_map'][z]['d'])

    e = np.add(e, tetra_seq_map['normalized_params_map'][z]['e'])

    f = np.add(f, tetra_seq_map['normalized_params_map'][z]['f'])

    g = np.add(g, tetra_seq_map['normalized_params_map'][z]['g'])

    h = np.add(h, tetra_seq_map['normalized_params_map'][z]['h'])

    i = np.add(i, tetra_seq_map['normalized_params_map'][z]['i'])

    j = np.add(j, tetra_seq_map['normalized_params_map'][z]['j'])

    k = np.add(k, tetra_seq_map['normalized_params_map'][z]['k'])

    l = np.add(l,tetra_seq_map['normalized_params_map'][z]['l'])

    m = np.add(m,tetra_seq_map['normalized_params_map'][z]['m'])

    n = np.add(n,tetra_seq_map['normalized_params_map'][z]['n'])

    o = np.add(o,tetra_seq_map['normalized_params_map'][z]['o'])

    p = np.add(p, tetra_seq_map['normalized_params_map'][z]['p'])

    q = np.add(q, tetra_seq_map['normalized_params_map'][z]['q'])

    r = np.add(r, tetra_seq_map['normalized_params_map'][z]['r'])


for z in range(len(sequence_map_cds)):
    a_cds = np.add(a_cds, tetra_cds_seq_map['normalized_params_map'][z]['a'])
    b_cds = np.add(b_cds, tetra_cds_seq_map['normalized_params_map'][z]['b'])
    c_cds = np.add(c_cds, tetra_cds_seq_map['normalized_params_map'][z]['c'])
    d_cds = np.add(d_cds, tetra_cds_seq_map['normalized_params_map'][z]['d'])
    e_cds = np.add(e_cds, tetra_cds_seq_map['normalized_params_map'][z]['e'])
    f_cds = np.add(f_cds, tetra_cds_seq_map['normalized_params_map'][z]['f'])
    g_cds = np.add(g_cds, tetra_cds_seq_map['normalized_params_map'][z]['g'])
    h_cds = np.add(h_cds, tetra_cds_seq_map['normalized_params_map'][z]['h'])
    i_cds = np.add(i_cds, tetra_cds_seq_map['normalized_params_map'][z]['i'])
    j_cds = np.add(j_cds, tetra_cds_seq_map['normalized_params_map'][z]['j'])
    k_cds = np.add(k_cds, tetra_cds_seq_map['normalized_params_map'][z]['k'])
    l_cds = np.add(l_cds,tetra_cds_seq_map['normalized_params_map'][z]['l'])
    m_cds = np.add(m_cds,tetra_cds_seq_map['normalized_params_map'][z]['m'])
    n_cds = np.add(n_cds,tetra_cds_seq_map['normalized_params_map'][z]['n'])
    o_cds = np.add(o_cds,tetra_cds_seq_map['normalized_params_map'][z]['o'])
    p_cds = np.add(p_cds,tetra_cds_seq_map['normalized_params_map'][z]['p'])
    q_cds = np.add(q_cds,tetra_cds_seq_map['normalized_params_map'][z]['q'])
    r_cds = np.add(r_cds,tetra_cds_seq_map['normalized_params_map'][z]['r'])

a = [z/len(sequence_map_tss) for z in a]
b = [z/len(sequence_map_tss) for z in b]
c = [z/len(sequence_map_tss) for z in c]
d = [z/len(sequence_map_tss) for z in d]
e = [z/len(sequence_map_tss) for z in e]
f = [z/len(sequence_map_tss) for z in f]
g = [z/len(sequence_map_tss) for z in g]
h = [z/len(sequence_map_tss) for z in h]
i = [z/len(sequence_map_tss) for z in i]
j = [z/len(sequence_map_tss) for z in j]
k = [z/len(sequence_map_tss) for z in k]
l = [z/len(sequence_map_tss) for z in l]
m = [z/len(sequence_map_tss) for z in m]
n = [z/len(sequence_map_tss) for z in n]
o = [z/len(sequence_map_tss) for z in o]
p = [z/len(sequence_map_tss) for z in p]
q = [z/len(sequence_map_tss) for z in q]
r = [z/len(sequence_map_tss) for z in r]

tss_df = pd.DataFrame(
    {'a': a,
     'b': b,
     'c': c,
     'd': d,
     'e': e,
     'f': f,
     'g': g,
     'h': h,
     'i': i,
     'j': j,
     'k': k,
     'l': l,
     'm': m,
     'n': n,
     'o': o,
     'p': p,
     'q': q,
     'r': r})

a_cds = [z/len(sequence_map_cds) for z in a_cds]
b_cds = [z/len(sequence_map_cds) for z in b_cds]
c_cds = [z/len(sequence_map_cds) for z in c_cds]
d_cds = [z/len(sequence_map_cds) for z in d_cds]
e_cds = [z/len(sequence_map_cds) for z in e_cds]
f_cds = [z/len(sequence_map_cds) for z in f_cds]
g_cds = [z/len(sequence_map_cds) for z in g_cds]
h_cds = [z/len(sequence_map_cds) for z in h_cds]
i_cds = [z/len(sequence_map_cds) for z in i_cds]
j_cds = [z/len(sequence_map_cds) for z in j_cds]
k_cds = [z/len(sequence_map_cds) for z in k_cds]
l_cds = [z/len(sequence_map_cds) for z in l_cds]
m_cds = [z/len(sequence_map_cds) for z in m_cds]
n_cds = [z/len(sequence_map_cds) for z in n_cds]
o_cds = [z/len(sequence_map_cds) for z in o_cds]
p_cds = [z/len(sequence_map_cds) for z in p_cds]
q_cds = [z/len(sequence_map_cds) for z in q_cds]
r_cds = [z/len(sequence_map_cds) for z in r_cds]

cds_df = pd.DataFrame(
    {'a_cds': a_cds,
     'b_cds': b_cds,
     'c_cds': c_cds,
     'd_cds': d_cds,
     'e_cds': e_cds,
     'f_cds': f_cds,
     'g_cds': g_cds,
     'h_cds': h_cds,
     'i_cds': i_cds,
     'j_cds': j_cds,
     'k_cds': k_cds,
     'l_cds': l_cds,
     'm_cds': m_cds,
     'n_cds': n_cds,
     'o_cds': o_cds,
     'p_cds': p_cds,
     'q_cds': q_cds,
     'r_cds': r_cds})


tss_df.to_csv("tss_df_tetra_syne_wonorm.csv")
cds_df.to_csv("cds_df_tetra_syne_wonorm.csv")


fig, axs = plt.subplots(18)

axs[0].plot(a,label = 'a')
axs[0].plot(a_cds)
axs[1].plot(b,label = 'b')
axs[1].plot(b_cds)
axs[2].plot(c,label = 'c')
axs[2].plot(c_cds)
axs[3].plot(d,label = 'd')
axs[3].plot(d_cds)
axs[4].plot(e,label = 'e')
axs[4].plot(e_cds)
axs[5].plot(f,label = 'f')
axs[5].plot(f_cds)
axs[6].plot(g,label = 'g')
axs[6].plot(g_cds)
axs[7].plot(h,label = 'h')
axs[7].plot(h_cds)
axs[8].plot(i,label = 'i')
axs[8].plot(i_cds)
axs[9].plot(j,label = 'j')
axs[9].plot(j_cds)
axs[10].plot(k,label = 'k')
axs[10].plot(k_cds)
axs[11].plot(l,label = 'l')
axs[11].plot(l_cds)
axs[12].plot(m,label = 'm')
axs[12].plot(m_cds)
axs[13].plot(n,label = 'n')
axs[13].plot(n_cds)
axs[14].plot(o,label = 'o')
axs[14].plot(o_cds)
axs[15].plot(p,label = 'p')
axs[15].plot(p_cds)
axs[16].plot(q,label = 'q')
axs[16].plot(q_cds)
axs[17].plot(r,label = 'r')
axs[17].plot(r_cds)

for i in range(0,18):
    axs[i].legend()
    # axs[i].xticks(np.arange(-500, 500, 100))
    # plt.yticks(np.arange(0, 1, 0.2))

plt.show()

import sys
sys.exit()
