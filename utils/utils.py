import numpy as np
import os, six, sys, subprocess
import tensorflow as tf
import random
from tqdm import tqdm
import pandas as pd
from pathlib import Path

# ------------- one hot encoding of RNA sequences -----------------#
def one_hot(seq):
    RNN_seq = seq
    BASES = 'AUCG'
    bases = np.array([base for base in BASES])
    feat = np.concatenate(
        [[(bases == base.upper()).astype(int)] if str(base).upper() in BASES else np.array([[-1] * len(BASES)]) for base
         in RNN_seq])

    return feat


def z_mask(seq_len):
    mask = np.ones((seq_len, seq_len))
    return np.triu(mask, 2)

def l_mask(inp, seq_len):
    mask = np.ones((seq_len, seq_len))
    return np.triu(mask, 1)

def get_data(seq, rna_id, args):
    seq_len = len(seq)
    one_hot_feat = one_hot(seq)

    with open(os.path.splitext(args.inputs)[0] + '.pssm') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None).values
    seq = ['U' if k == 'T' else k for k in temp[:, 0]]
    profile = temp[:, 1:5].astype(float)
    off_set = np.zeros((len(seq), profile.shape[1])) + 0.3
    for k, K in enumerate(seq):
        try:
            off_set[k, BASES.index(K)] = 8.7
        except:
            pass
    profile += off_set
    profile /= np.sum(profile, axis=1, keepdims=True)
    profile = -np.log(profile)

    profile_one_hot = np.concatenate([profile, one_hot_feat], axis=1)

############ load base-pair prob form linearpartition ##############################
    try:
        with open(os.path.splitext(args.inputs)[0] + '.prob', 'r') as f:
            prob = pd.read_csv(f, delimiter=None, delim_whitespace=True, header=None, skiprows=[0]).values
        bp_prob_lp =  np.zeros((len(seq), len(seq)))
        for i in prob:
            bp_prob_lp[int(i[0])-1, int(i[1])-1] = i[2]
        bp_prob_lp = bp_prob_lp + np.transpose(bp_prob_lp)
    except:
        print("linearpartition output missing",rna_id)
        bp_prob_lp =  np.zeros((len(seq), len(seq)))

############ load dca obtained from gremlin ##############################
    try:
        with open(os.path.splitext(args.inputs)[0] + '.dca') as f:
            temp4 = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0], usecols=[0,1,2]).values
        #print(temp4.shape)
        dca_output = np.zeros((len(seq), len(seq)))
        for k in temp4:
            if abs(int(k[0]) - int(k[1])) < 4:
                dca_output[int(k[0]), int(k[1])] = 0.01*k[2]
            else:
                dca_output[int(k[0]), int(k[1])] = k[2]
        dca_output = dca_output + np.transpose(dca_output)
    except:
        print("dca missing", rna_id)
        dca_output = np.zeros((len(seq), len(seq)))

    zero_mask = z_mask(seq_len)[None, :, :, None]
    label_mask = l_mask(profile_one_hot, seq_len)
    temp = profile_one_hot[None, :, :]
    temp = np.tile(temp, (temp.shape[1], 1, 1))
    feature = np.concatenate([temp, np.transpose(temp, [1, 0, 2])], 2)
    feature = np.concatenate([feature, np.expand_dims(dca_output, axis=2)], axis=2)
    feature = np.concatenate([feature, np.expand_dims(bp_prob_lp, axis=2)], axis=2)

    assert feature.shape==(seq_len,seq_len, 18)

    return seq_len, [i for i in (feature.astype(float)).flatten()], [i for i in zero_mask.flatten()], [i for i in label_mask.flatten()], [i for i in label_mask.flatten()]

def _int64_feature(value):
    if not isinstance(value, list) and not isinstance(value, np.ndarray):
        value = [value]

    return tf.train.Feature(int64_list=tf.train.Int64List(value=value))


def _float_feature(value):
    if not isinstance(value, list) and not isinstance(value, np.ndarray):
        value = [value]

    return tf.train.Feature(float_list=tf.train.FloatList(value=value))


def _bytes_feature(value):
    """Wrapper for inserting bytes features into Example proto."""
    if isinstance(value, six.string_types):
        value = six.binary_type(value, encoding='utf-8')
    return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))

def create_tfr_files(args):

    print('\nPreparing tfr records file for SPOT-RNA2:')
    path_tfrecords = os.path.splitext(args.inputs)[0] + ".tfrecords"
    with open(args.inputs) as file:
        input_data = [line.strip() for line in file.read().splitlines() if line.strip()]

    count = int(len(input_data)/2)

    ids = [input_data[2*i][1:].strip() for i in range(count)]
    
    with tf.io.TFRecordWriter(path_tfrecords) as writer:
        for i in tqdm(range(len(ids))):
            name     = input_data[2*i].replace(">", "") 
            sequence = input_data[2*i+1].replace(" ", "").replace("T", "U").upper()
            #print(sequence[-1])
            
            #print(len(sequence), name)                
            seq_len, feature, zero_mask, label_mask, true_label = get_data(sequence, name, args)

            example = tf.train.Example(features=tf.train.Features(feature={'rna_name': _bytes_feature(name),
                                                                           'seq_len': _int64_feature(seq_len),
                                                                           'feature': _float_feature(feature),
                                                                           'zero_mask': _float_feature(zero_mask),
                                                                           'label_mask': _float_feature(label_mask),
                                                                           'true_label': _float_feature(true_label)}))

            writer.write(example.SerializeToString())

    writer.close()

# ----------------------- hair pin loop assumption i - j < 2 --------------------------------#
def hair_pin_assumption(pred_pairs):
    pred_pairs_all = [i[:2] for i in pred_pairs]
    bad_pairs = []
    for i in pred_pairs_all:
        if abs(i[0] - i[1]) < 3:
            bad_pairs.append(i)
    return bad_pairs

def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, str):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

def type_pairs(pairs, sequence):
    sequence = [i.upper() for i in sequence]
    # seq_pairs = [[sequence[i[0]],sequence[i[1]]] for i in pairs]

    AU_pair = []
    GC_pair = []
    GU_pair = []
    other_pairs = []
    for i in pairs:
        if [sequence[i[0]],sequence[i[1]]] in [["A","U"], ["U","A"]]:
            AU_pair.append(i)
        elif [sequence[i[0]],sequence[i[1]]] in [["G","C"], ["C","G"]]:
            GC_pair.append(i)
        elif [sequence[i[0]],sequence[i[1]]] in [["G","U"], ["U","G"]]:
            GU_pair.append(i)
        else:
            other_pairs.append(i)
    watson_pairs_t = AU_pair + GC_pair
    wobble_pairs_t = GU_pair
    other_pairs_t = other_pairs
        # print(watson_pairs_t, wobble_pairs_t, other_pairs_t)
    return watson_pairs_t, wobble_pairs_t, other_pairs_t

# ----------------------- find multiplets pairs--------------------------------#
def multiplets_pairs(pred_pairs):

    pred_pair = [i[:2] for i in pred_pairs]
    temp_list = flatten(pred_pair)
    temp_list.sort()
    new_list = sorted(set(temp_list))
    dup_list = []
    for i in range(len(new_list)):
        if (temp_list.count(new_list[i]) > 1):
            dup_list.append(new_list[i])

    dub_pairs = []
    for e in pred_pair:
        if e[0] in dup_list:
            dub_pairs.append(e)
        elif e[1] in dup_list:
            dub_pairs.append(e)

    temp3 = []
    for i in dup_list:
        temp4 = []
        for k in dub_pairs:
            if i in k:
                temp4.append(k)
        temp3.append(temp4)
        
    return temp3

def multiplets_free_bp(pred_pairs, y_pred):
    L = len(pred_pairs)
    multiplets_bp = multiplets_pairs(pred_pairs)
    save_multiplets = []
    while len(multiplets_bp) > 0:
        remove_pairs = []
        for i in multiplets_bp:
            save_prob = []
            for j in i:
                save_prob.append(y_pred[j[0], j[1]])
            remove_pairs.append(i[save_prob.index(min(save_prob))])
            save_multiplets.append(i[save_prob.index(min(save_prob))])
        pred_pairs = [k for k in pred_pairs if k not in remove_pairs]
        multiplets_bp = multiplets_pairs(pred_pairs)
    save_multiplets = [list(x) for x in set(tuple(x) for x in save_multiplets)]
    assert L == len(pred_pairs)+len(save_multiplets)
    #print(L, len(pred_pairs), save_multiplets)
    return pred_pairs, save_multiplets
        
def output_mask(seq, NC=True):
    if NC:
        include_pairs = ['AU', 'UA', 'GC', 'CG', 'GU', 'UG', 'CC', 'GG', 'AG', 'CA', 'AC', 'UU', 'AA', 'CU', 'GA', 'UC']
    else:
        include_pairs = ['AU', 'UA', 'GC', 'CG', 'GU', 'UG']
    mask = np.zeros((len(seq), len(seq)))
    for i, I in enumerate(seq):
        for j, J in enumerate(seq):
            if str(I) + str(J) in include_pairs:
                mask[i, j] = 1
    return mask

def ct_file_output(pairs, seq, id, save_result_path):

    col1 = np.arange(1, len(seq) + 1, 1)
    col2 = np.array([i for i in seq])
    col3 = np.arange(0, len(seq), 1)
    col4 = np.append(np.delete(col1, 0), [0])
    col5 = np.zeros(len(seq), dtype=int)

    for i, I in enumerate(pairs):
        col5[I[0]] = int(I[1]) + 1
        col5[I[1]] = int(I[0]) + 1
    col6 = np.arange(1, len(seq) + 1, 1)
    temp = np.vstack((np.char.mod('%d', col1), col2, np.char.mod('%d', col3), np.char.mod('%d', col4),
                      np.char.mod('%d', col5), np.char.mod('%d', col6))).T

    np.savetxt(os.path.join(save_result_path, str(id))+'.ct', (temp), delimiter='\t\t', fmt="%s", header=str(len(seq)) + '\t\t' + str(id) + '\t\t' + 'SPOT-RNA2 output\n' , comments='')

    return

def bpseq_file_output(pairs, seq, id, save_result_path):

    col1 = np.arange(1, len(seq) + 1, 1)
    col2 = np.array([i for i in seq])
    #col3 = np.arange(0, len(seq), 1)
    #col4 = np.append(np.delete(col1, 0), [0])
    col5 = np.zeros(len(seq), dtype=int)

    for i, I in enumerate(pairs):
        col5[I[0]] = int(I[1]) + 1
        col5[I[1]] = int(I[0]) + 1
    #col6 = np.arange(1, len(seq) + 1, 1)
    temp = np.vstack((np.char.mod('%d', col1), col2, np.char.mod('%d', col5))).T
    #os.chdir(save_result_path)
    #print(os.path.join(save_result_path, str(id[0:-1]))+'.spotrna')
    np.savetxt(os.path.join(save_result_path, str(id))+'.bpseq', (temp), delimiter=' ', fmt="%s", header='#' + str(id) , comments='')

    return

def lone_pair(pairs):
    lone_pairs = []
    pairs.sort()
    for i, I in enumerate(pairs):
        if ([I[0] - 1, I[1] + 1] not in pairs) and ([I[0] + 1, I[1] - 1] not in pairs):
            lone_pairs.append(I)

    return lone_pairs

def prob_to_secondary_structure(ensemble_outputs, label_mask, seq, name, args):
    #save_result_path = 'outputs'
    Threshold = 0.795
    label_mask = np.triu(np.ones((len(seq), len(seq))),1)
    test_output = ensemble_outputs
    mask = output_mask(seq)
    inds = np.where(label_mask == 1)
    y_pred = np.zeros(label_mask.shape)
    for i in range(test_output.shape[0]):
        y_pred[inds[0][i], inds[1][i]] = test_output[i]
    #y_pred = np.multiply(y_pred, mask)

    tri_inds = np.triu_indices(y_pred.shape[0], k=1)

    out_pred = y_pred[tri_inds]
    outputs = out_pred[:, None]
    seq_pairs = [[tri_inds[0][j], tri_inds[1][j], ''.join([seq[tri_inds[0][j]], seq[tri_inds[1][j]]])] for j in
                 range(tri_inds[0].shape[0])]

    outputs_T = np.greater_equal(outputs, Threshold)
    pred_pairs = [i for I, i in enumerate(seq_pairs) if outputs_T[I]]
    pred_pairs = [i[:2] for i in pred_pairs]
    pred_pairs, save_multiplets = multiplets_free_bp(pred_pairs, y_pred)
    
    if args.outputs=='outputs/':
        output_path = os.path.join(Path(os.path.dirname(os.path.realpath(__file__))).parent, args.outputs)
    else:
        output_path = args.outputs

    ct_file_output(pred_pairs, seq, name, output_path)
    bpseq_file_output(pred_pairs, seq, name, output_path)
    np.savetxt(output_path + '/'+ name +'.prob', y_pred, delimiter='\t')
    
    if args.motifs:
        try:
            os.chdir(args.outputs)
            p = subprocess.Popen(['perl', os.path.join(Path(os.path.dirname(os.path.realpath(__file__))).parent, 'utils/bpRNA.pl'), os.path.join(args.outputs, name + '.bpseq')])
        except:
            print('\nUnable to run bpRNA script;\nplease refer to "https://github.com/hendrixlab/bpRNA/" for system requirments to use bpRNA')
        #os.chdir('../')

    return
