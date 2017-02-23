# ⓒ Copyright IBM Corp. 2017
import numpy as np
import sys
from graphframes import *

# from utils import *

encoding_map = {'A': 0, 'C': 1, 'T': 2, 'G': 3}  # chosen to make rc(x) == XOR(x)
decoding_lst = ['A', 'C', 'T', 'G']


def encode(x):                  
    code = 0
    for ch in x:
        code *= 4
        code += encoding_map[ch]
    return code
    
def get_kmers(seq, k):
    return [seq[i:i+k] for i in range(len(seq)-k+1)]  

def get_minimizer(kmer, l):
    assert  l <= len(kmer)
    return min(get_kmers(kmer, l))

def getKmerToNextCharCounts(x,k):
    """ given a string x and k val, gets legal kmers in the read, then
        filters kmers containing N chars, translates kmers to (K,V) pairs
        where K is the kmer, and V is a list with only the next character set to one
        e.g., [0,1,0,0] if the (k+1)th character is a C
        Finally, encodes K as its 2-bit version inside each (K,V) pair
    """
    allKmers = get_kmers(x, k+1)   
    filtered = filter(lambda y: "N" not in y and "n" not in y, allKmers)
    res = []
    for a in filtered:
        b = np.array([0,0,0,0], dtype='uint16')
        b[encoding_map[a[-1]]] = 1
        res.append((a[:-1], b)) #(encode(a[:-1]),b))
    return res

def my_filter(a):
    """ counts the number of non-zeros in a numpy array
        a[0] is expected to be a kmer and a[1] is a numpy array
        this function is needed because spark filter operations require 
        function arguments that take one argument only
    """
    return np.sum(a[1] != 0) > 1

def build_partial_junctions_set():
    """ based on 
        http://tech.magnetic.com/2016/01/bloom-filter-assisted-joins-using-pyspark.html
        but I replaced BFs with sets for now; I don't know if the closure stuff is 
        needed or helpful...
    """
    def _build_partial_junc_set(junctions):
        junc_set = set([])
        for junc in junctions:
            junc_set.add(junc[0])
        yield (None, junc_set)

    return _build_partial_junc_set

def merge_sets(set1, set2): # not sure if this function is needed
    return set1.union(set2)

def filter_reads_by_junctions(partition):
    """ given junctions broadcast, 
        out of all reads in each partition
        yields only reads that include some junction
    """
    juncs = juncs_broadcast.value

    for read_name, seq in partition:
        allKmers = get_kmers(seq,k)
        for kmer in allKmers:
            if kmer in juncs:
                yield (read_name, seq) 


# def build_anchors_graph(reads_RDD, juncs):


def launch_spark_job():
    from pyspark import SparkContext, SparkConf
    from pyspark.sql import SQLContext
    from pyspark.sql.functions import concat, col, lit

    readFile = sys.argv[1]
    k = int(sys.argv[2])
    num_partitions = int(sys.argv[3])
    conf = SparkConf().setAppName("reads Loader"+str(num_partitions))
    sc = SparkContext(conf=conf)
    sc.addPyFile("utils.py")
    sc.setCheckpointDir("hdfs://doop-mng1.haifa.ibm.com:8020/projects/Store_Analytics/SparkCheckPoints")
    import utils
    # from utils import map_read_to_anchors_list, convert_anchors_list_to_seq_edges
    readLines = (
        sc.newAPIHadoopFile(
            readFile,
            'org.apache.hadoop.mapreduce.lib.input.TextInputFormat',
            'org.apache.hadoop.io.LongWritable',
            'org.apache.hadoop.io.Text',
            conf={'textinputformat.record.delimiter': '@'}) #, 
        .map(lambda delim_lines_tup: delim_lines_tup[1])  # keeps just the lines and not the @ delimiter
        .filter(lambda x: x.startswith("SRR")) # gets rid of entries due to '@' appearing in the wrong line
        .map(lambda x: x.split("\n")[:2]) # splits the lines, keeps only the first two
        .filter(lambda x: len(x)==2) # git rid of any cut off records
        .repartition(num_partitions)
        # .cache()
        )

    print("----------------------there are %i reads" % (readLines.count()))

    # get new RDD including lists of kmers (with no Ns), (k+1)mers
    kmers = (
        readLines
        .map(lambda entry: entry[1])
        .flatMap(lambda read: getKmerToNextCharCounts(read, k)) 
    )
   
    print("----------------------there are %i kmers instances" % (kmers.count()))

    kmers_with_exts = (
        kmers.reduceByKey(func=lambda x,y: x+y) 
    )

    print("----------------------there are %i distinct kmers" % (kmers_with_exts.count()))

    junctions = kmers_with_exts.filter(lambda kmer_tup: my_filter(kmer_tup))

    print("----------------------there are %i junctions" % junctions.count())

    # for i in junctions.take(10):
    #     if sum(i[1])>1:
    #         print i
    

    generate_juncs = build_partial_junctions_set()
    junctions_set_rdd = (
        junctions.mapPartitions(generate_juncs)
        .reduceByKey(merge_sets)
        .collect()
    )

    juncs_broadcast = sc.broadcast(junctions_set_rdd[0][1])
    print("----------------------there are %i junctions in broadcast" % len(juncs_broadcast.value))

    # build edge set rdd, filter out edges including a junction at some end

    def read_line_map_function(read_line):
        return utils.map_read_to_anchors_list(read_line[1], k-10, 10, juncs_broadcast.value)

    edges_rdd = (
        readLines.map(lambda read_line: read_line_map_function(read_line))
        .flatMap(lambda anchors: utils.convert_anchors_list_to_seq_edges(anchors), preservesPartitioning=True)
        .filter(lambda (a,b,c): a not in juncs_broadcast.value and b not in juncs_broadcast.value)
        )

    print("----------------------there are %i total edges" % edges_rdd.count())

    # create SQLContext to be able to create dataFrame from rdd
    sqc = SQLContext(sc)
    edges_df = sqc.createDataFrame(edges_rdd, ["src", "dst", "overlap"])
    vertices_df = edges_df.select(concat(col("src"), lit(" "), col("dst")).alias('id')).dropDuplicates()
    g = GraphFrame(vertices_df, edges_df)

    # vertices_df.agg(*[count(c).alias(c) for c in vertices_df.columns]).show()

    print("----------------------there are %i total vertices" % vertices_df.count())

    # get connected components of remaining graph 

    result = g.connectedComponents()
    result.select("id", "component").orderBy("component").show()

    # reads_with_junctions = (
    #     readLines.mapPartitions(filter_reads_by_junctions).collect()
    # )


    # print("----------------------there are %i total reads with junctions" % len(reads_with_junctions))
    # print("----------------------there are %i unique reads with junctions" % len(set(reads_with_junctions)))


    # collapse lists of extensions to counts for each letter
    
    # for i in junctions.take(100):
    #     if sum(i[1])>1:
    #         print i

if __name__ == "__main__":
    launch_spark_job()