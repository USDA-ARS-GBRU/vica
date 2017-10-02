#!/usr/bin/env python3
'''tfrecord_maker.py: a module to create tf record files'''

def csv_to_tfrecords(input_filename, output_filename, delimiter=","):
    '''convert csv file to a dense tfrecord file, label should be in first column freatures should be in all other columns'''
    print("Start to convert {} to {}".format(input_filename, output_filename))
    writer = tf.python_io.TFRecordWriter(output_filename)

    for line in open(input_filename, "r"):
        data = line.split(delimiter)
        label = float(data[0])
        features = [float(i) for i in data[1:]]

    example = tf.train.Example(features=tf.train.Features(feature={
        "label":
        tf.train.Feature(float_list=tf.train.FloatList(value=[label])),
        "features":
        tf.train.Feature(float_list=tf.train.FloatList(value=features)),
        }))
    writer.write(example.SerializeToString())

    writer.close()
    print("Successfully convert {} to {}".format(input_filename,
                                               output_filename))
def libsvm_to_tfrecords(input_filename, output_filename):
    print("Start to convert {} to {}".format(input_filename, output_filename))
    writer = tf.python_io.TFRecordWriter(output_filename)

    for line in open(input_filename, "r"):
        data = line.split(" ")
        label = float(data[0])
        ids = []
        values = []
        for fea in data[1:]:
            id, value = fea.split(":")
            ids.append(int(id))
            values.append(float(value))

    # Write each example one by one
    example = tf.train.Example(features=tf.train.Features(feature={
        "label":
        tf.train.Feature(float_list=tf.train.FloatList(value=[label])),
        "ids": tf.train.Feature(int64_list=tf.train.Int64List(value=ids)),
        "values": tf.train.Feature(float_list=tf.train.FloatList(value=values))
    }))

    writer.write(example.SerializeToString())
    writer.close()
    print("Successfully convert {} to {}".format(input_filename,
                                               output_filename))
def parse_sketch_to_tfrecord(sketchfile, tfrecordfile, labeldict):
    writer = tf.python_io.TFRecordWriter(tfrecordfile)
    try:
        with open(sketchfile, 'r') as f:
            label = None
            taxidlist = []
            matchlist = []
            for line in f:
                if line.strip() == '':
                    # write tf record examples
                    pass
                    # 3 reset lists and label
                    taxidlist = []
                    matchlist = []
                    label = None
                elif line.startswith("Query:"):
                    ll = line.strip().split("\t")
                    key1 = ll[0].split(":")[1].split()[0]
                    label = labeldict[key1]
                elif line.startswith("WKID"):
                    next
                else:
                    ll2 = line.strip().split("\t")
                    taxidlist.append(int(ll2[5]))
                    matchlist.append(int(ll2[3]))
    except IOError:
        print("could not parse of the sketch file %s" % (file))
