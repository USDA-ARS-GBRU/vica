#!/usr/bin/env python3
'''tfrecord_maker.py: a module to create tf record files'''
import tensorflow as tf
import subprocess
import itertools
import tempfile

def external_sort(infile, outfile, sep, size=1073741824, key=1):
    options = ['sort', '--field-separator=' + sep, '--buffer-size=' + size , '--key=' + key, '>', '--output=' + outfile ]
    subprocess.run(options)

def external_join(infile1, infile2, outfile, sep):
    options = ['join', '-t ' + sep, '-1 1', '-1 2', infile1, infile2, '> ' + outfile ]
    subprocess.run(options)

# def parse_mapping(file):
#         '''parse two colum file in the format: label <tab> id'''
#         mdict = {}
#         with open(file, 'r') as f:
#             for line in f:
#                 ll = line.strip()split()
#                 mdict[ll[1]] = ll[0]
#         return mdict


def combine_libsvm_csv_to_tfrecords(label, libsvmfile, csvfile,  tfrecordfile):
    '''convert a sorted libsvm file of sparse data and a sorted csv file into a tfrecord file'''
    writer = tf.python_io.TFRecordWriter(tfrecordfile)
    with open(libsvmfile) as file1, open(csvfile) as file2:
        for line1, line2 in izip(file1, file2):
            ll1 = line1.split(' ')
            ll2 = line2.split(',')
            assert ll1[0] == ll2[0], "Lines in the svmlib file and the csv file do not match"
            densefeatures = [float(i) for i in ll2[1:]]
            values = []
            ids =  []
            for fea in ll1[1:]:
                id, value = fea.split(":")
                ids.append(int(id))
                values.append(float(value))
            example = tf.train.Example(features=tf.train.Features(feature={
                "label":
                    tf.train.Feature(float_list=tf.train.FloatList(value=[label])),
                "densefeatures":
                    tf.train.Feature(float_list=tf.train.FloatList(value=densefeatures)),
                "minhashids":
                    tf.train.Feature(int64_list=tf.train.Int64List(value=ids)),
                "minhashvalues":
                    tf.train.Feature(float_list=tf.train.FloatList(value=values))
                }))
            writer.write(example.SerializeToString())

            writer.close()
            print("Successfully converted {} and {} to {}".format(libsvmfile, csvfile,
                                                                  tfrecordfile))


def convert_to_tfrecords(kmerfile, codonfile, minhashfile, tfrecordfile, label):
    dtemp = tempfile.mkdtemp()
    ksorted = os.path.join(dtemp, kmer_sorted.csv)
    csorted = os.path.join(dtemp, codon_sorted.csv)
    densefile = os.path.join(dtemp, dense.csv)
    external_sort(infile=kmerfile, oufile=ksorted, sep=",")
    external_sort(infile=codonfile, oufile=csorted, sep=",")
    external_join(infile1=ksorted, infile2=csorted, outfile=densefile, sep=",")
    combine_libsvm_csv_to_tfrecords(label=label, libsvmfile=minhashfile,
                                    csvfile=densefile,  tfrecordfile=tfrecordfile)

def main():

    parser = argparse.ArgumentParser(description='A script to generate a tfrecord file from feature files')
    parser.add_argument('--kmerin', help="A csv of kmer frequencies")
    parser.add_argument('--codonin', help="A csv of codon frequencies")
    parser.add_argument('--minhashin', help="A libsvm file of taxonids and minhash matches")
    parser.add_argument('--outfile', help= "A tfrecord file")
    parser.add_argument('--label', help= "An interger label of the class")
    args = parser.parse_args()

    convert_to_tfrecords(kmerfile=args.kmerin, codonfile=args.codonin,
                         minhashfile=args.minhashin, tfrecordfile=args.utfile,
                         label=args.label)

if __name__ == '__main__':
    main()
