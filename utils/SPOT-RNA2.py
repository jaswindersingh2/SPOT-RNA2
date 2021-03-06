import tensorflow as tf
import numpy as np
import os
from tqdm import tqdm
import argparse
from utils import create_tfr_files, prob_to_secondary_structure
import time
start = time.time()
from argparse import RawTextHelpFormatter
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--inputs', default='inputs/single_seq.fasta', type=str, help='Path to input file in fasta format, accept multiple sequences as well in fasta format; default = ''inputs/2zzm-1-B.fasta''\n', metavar='')
parser.add_argument('--outputs',default='outputs/', type=str, help='Path to output files; SPOT-RNA outputs at least three files .ct, .bpseq, and .prob files; default = ''outputs/\n', metavar='')
parser.add_argument('--gpu', default=1, type=int, help='To run on GPU, specifiy GPU number. If only one GPU in computer specifiy 0; default = -1 (no GPU)\n', metavar='')
parser.add_argument('--plots',default=False, type=bool, help='Set this to "True" to get the 2D plots of predicted secondary structure by SPOT-RNA; default = False\n', metavar='')
parser.add_argument('--motifs',default=False, type=bool, help='Set this to "True" to get the motifs of predicted secondary structure by SPOT-RNA; default = False\n', metavar='')
#parser.add_argument('--NC',default=True, type=bool, help='Set this to "False" to predict only canonical pairs; default = True\n', metavar='')
args = parser.parse_args()

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

base_path = os.path.dirname(os.path.realpath(__file__))

create_tfr_files(args)

with open(args.inputs) as file:
    input_data = [line.strip() for line in file.read().splitlines() if line.strip()]

count = int(len(input_data)/2)

ids = [input_data[2*i].replace(">", "") for i in range(count)]
sequences = {}
for i,I in enumerate(ids):
    sequences[I] = input_data[2*i+1].replace(" ", "").replace("T", "U").upper()

os.environ["CUDA_VISIBLE_DEVICES"]= str(args.gpu)
#os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
NUM_MODELS = 4

test_loc = [os.path.splitext(args.inputs)[0] + ".tfrecords"]

outputs = {}
mask = {}
def sigmoid(x):
    return 1/(1+np.exp(-np.array(x, dtype=np.float128)))

#for MODEL in range(NUM_MODELS):
for MODEL in [0, 1, 2, 3]:
#for MODEL in [0, 1, 2, 3]:
    print(MODEL)
    config = tf.ConfigProto()
    #config.gpu_options.allow_growth = True
    config.allow_soft_placement=True
    config.log_device_placement=False
    print('\nPredicting for SPOT-RNA2 model '+str(MODEL))
    with tf.Session(config=config) as sess:
        saver = tf.train.import_meta_graph(os.path.join(base_path, 'models_ckps'+'/model_'+str(MODEL)+'.meta'))
        saver.restore(sess, os.path.join(base_path, 'models_ckps'+'/model_'+str(MODEL)))
        graph = tf.get_default_graph()
        init_test =  graph.get_operation_by_name('make_initializer_1')
        tmp_out = graph.get_tensor_by_name('output_FC/fully_connected/BiasAdd:0')
        name_tensor = graph.get_tensor_by_name('tensors_1/component_0:0')
        RNA_name = graph.get_tensor_by_name('IteratorGetNext:0')
        label_mask = graph.get_tensor_by_name('IteratorGetNext:4')
        sess.run([init_test], feed_dict={name_tensor:test_loc})
        
        pbar = tqdm(total = count)
        for rna in ids:
            out = sess.run([tmp_out,RNA_name,label_mask],feed_dict={'dropout:0':1})
            out[1] = rna

            mask[out[1]] = out[2]
            
            if MODEL == 0:
                outputs[out[1]] = [sigmoid(out[0])]
            else:
                outputs[out[1]].append(sigmoid(out[0]))
            pbar.update(1)
        pbar.close()
    tf.reset_default_graph()


RNA_ids = [i for i in list(outputs.keys())]
ensemble_outputs = {}

print('\nPost Processing and Saving Output')
for i in RNA_ids:
    #print(i, mask[i].shape, len(sequences[i]))
    ensemble_outputs[i] = np.mean(outputs[i],0)
    prob_to_secondary_structure(ensemble_outputs[i], mask[i], sequences[i], i, args)

print('\nFinished!')
end = time.time()
print('\nProcesssing Time {} seconds'.format(end - start))
