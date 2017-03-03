// ⓒ Copyright IBM Corp. 2017

// import org.apache.spark.SparkContext
// import org.apache.spark.SparkContext._
// import org.apache.spark.SparkConf
// import org.apache.hadoop.conf.Configuration

// import org.apache.hadoop.io.LongWritable
// import org.apache.hadoop.io.Text
// import org.apache.hadoop.mapreduce.lib.input.TextInputFormat
import scala.collection.mutable.ListBuffer
import scala.util.control.Breaks._
import scala.collection.mutable.Queue

object SimpleApp {
  val encoding_map = Map('A'-> 0, 'C'-> 1, 'T'-> 2, 'G'-> 3)  // chosen to make rc(x) == XOR(x)

  def rc(dna: String): String = 
    dna.reverse.map(b => b match {
      case 'A' => 'T'
      case 'T' => 'A'
      case 'C' => 'G'
      case 'G' => 'C'
  })

  def getKmers(read: String, k: Int) = {
    for (i <- 0 until read.length - k +1) yield read.slice(i, i+k)     
  }

  def get_pos_kmer_tups(read: String, k: Int) = {
    for (i <- 0 until read.length - k +1) yield (i, read.substring(i, i+k))
  }

  def get_minimizer_w_pos(kmer_tup: (Int, String), l: Int) = {
    get_pos_kmer_tups(kmer_tup._2, l).min(Ordering.by{tup: (Int,String) => tup._2})
  }

  def get_stranded_kmerized_minimizers_list(seq: String, k: Int, l: Int) = {
    /* simple implemenation of getting minimizers while taking into consideration
    both forward and rc orientation. Not efficient because recomputes minimizers from 
    scratch inside each kmer window
    */
    val rc_seq  = rc(seq)
    val kmers = get_pos_kmer_tups(seq,k)
    var rc_kmers = get_pos_kmer_tups(rc_seq,k)
    var res = new ListBuffer[(Int,String)]()
	
//	println("----------------- kmers   " + kmers)
	
	   // AACCTTGCTACGTCCAAG
	   //   012345678901234567

    // change order, indices of rc_kmers to match kmers
    rc_kmers = rc_kmers.reverse.map(pos_kmer_tup => (rc_kmers.length - pos_kmer_tup._1 -1,pos_kmer_tup._2))
	
//	println("----------------- rc_kmers   " + rc_kmers) 
	
	var prev : Int = -1
    for (i <- 0 until kmers.length){
	  breakable {
		  val forward_min = get_minimizer_w_pos(kmers(i), l)
		  val rc_min  = get_minimizer_w_pos(rc_kmers(i), l)
		  // check there are enough bases to append at each end
		  val real_min: (Int,String) = List(forward_min, rc_min).minBy(_._2) //(Ordering.by{tup: (Int,String) => tup._2})
		  val padding_length: Int = (k-l)/2 // NB: choose l and k to be odd --> padding remains whole number
		  
		  val real_start_pos: Int = 
			if (real_min == forward_min) 
				real_min._1 + i
			else
				i + k - l - real_min._1
				
//		 println("-----------------" + real_min + " , " + prev)
		 
		   if (real_start_pos == prev) { break }
		   else { prev = real_start_pos }
		  
//		  println("-----------------" + real_min + " , " + prev)
		  if (real_start_pos - padding_length >= 0 && real_start_pos + padding_length + l -1 <= seq.length-1){
			// return the minimizer after padding it using either the forward of reverse sequence
			val min_tup: (Int, String) =  
			if (real_min == forward_min) 
			(real_start_pos - padding_length, 
			 seq.substring(real_start_pos - padding_length, real_start_pos + padding_length + l))
			else
			// (rc_seq.length - (real_start_pos + l - 1 + padding_length) - 1, 
      // rc_seq.substring(rc_seq.length - (real_start_pos + l - 1 + padding_length) -1, 
      //   rc_seq.length - (real_start_pos - padding_length)))
      (real_start_pos + l + padding_length - k, 
			  seq.substring(real_start_pos + l + padding_length - k, real_start_pos + l + padding_length))
			res += min_tup
		//	println("-----------------" + "," + i + "," + real_min +"," + padding_length +"," + real_start_pos +"," + min_tup)
		  } 
		}
    }
    res
  }

  def get_edge_list(read: String, anchors_array: ListBuffer[(Int,String)], k: Int, l: Int) = {
      // flatMap array to get edge pairs with orientation property; 
      // 1) pair v_i, v_{i+1} as edge s = (u,w)
      val edge_pairs = for (i <- 0 until anchors_array.length-1) yield (anchors_array(i),anchors_array(i+1))

      // 2) determine rep(s) = min(s,s') by extracting substring from start(v_i) to end(v_{i+1}) 
      // using array indices
      val edge_strings = for(i <- 0 until edge_pairs.length) yield read.substring(edge_pairs(i)._1._1, edge_pairs(i)._2._1 + k)
      val edge_rc_strings = for(i <- 0 until edge_pairs.length) yield rc(read.substring(edge_pairs(i)._1._1, edge_pairs(i)._2._1 + k))
      val edge_rep_strings = edge_strings.zip(edge_rc_strings).map{case (a,b) => List(a,b).min}    
      
      // 3) choose polarities of u and w relative to rep(s)
      val edge_list = for {i <- 0 until edge_rep_strings.length
          edge_rep = edge_rep_strings(i)
          left_rep = List(edge_rep.substring(0,k), rc(edge_rep.substring(0,k))).min
          right_rep = List(edge_rep.substring(edge_rep.length - k, edge_rep.length), 
          rc(edge_rep.substring(edge_rep.length - k, edge_rep.length))).min
          edge_tup = 
              (left_rep, right_rep,
               (edge_rep.substring(0,k) == left_rep,
               edge_rep.substring(edge_rep.length - k,edge_rep.length) == right_rep)
              )

      } yield edge_tup
      edge_list
  }

  // def map_read_to_anchors_list(seq: String, k: Int, l: Int, juncs: Map)

  def  getKmerToNextCharCounts(read: String, k: Int) = {
  /*""" given a string x and k val, gets legal kmers in the read, then
        filters kmers containing N chars, translates kmers to (K,V) pairs
        where K is the kmer, and V is a list with only the next character set to one
        e.g., [0,1,0,0] if the (k+1)th character is a C
    """
    */
    var allKmers = getKmers(read, k+1)
    val filtered_kmers = allKmers.filter(x => !x.contains("N") && !x.contains("n"))
    val last_char_value_arrays = filtered_kmers
    .map{kpomer => (kpomer.substring(0,kpomer.length-1), kpomer(kpomer.length-1))} // k.p.o. = "k plus one"
    .map{kmer_parts => 
      val char_counts = new Array[Int](4) // TODO: possibly reduce memory by changing type; effects reduceByKey later
      char_counts(encoding_map(kmer_parts._2)) = 1
      (kmer_parts._1, char_counts)}
    last_char_value_arrays
  }
  def filterTup(x: Array[Int]) : Boolean = {
    var cnt = 0
    // for (element <- x._2) {
    //   if (element > 0) count+=1
    // }
    // if (x._2.count(_ > 0) > 1) return true
    if (x(0) > 0) cnt += 1
    if (x(1) > 0) cnt += 1
    if (x(2) > 0) cnt += 1
    if (x(3) > 0) cnt += 1
    if (cnt > 1) return true 
    else return false 
    // return true
    // else return false
  }

def map_read_to_anchors_list(seq: String, k: Int, l: Int, juncs_set:Set[String]){		
	val mins_lst = get_stranded_kmerized_minimizers_list(seq,k,l)
	val allKmers = get_pos_kmer_tups(seq,k)
	val read_anchors = List[(Int,String)]()
	
	allKmers.filter(kmerT => juncs_set.contains(kmerT._2) || juncs_set.contains(rc(kmerT._2)) || mins_lst.contains(kmerT)).sorted
}


  def main(args: Array[String]) {
    // if (args.length < 3) {
    //   println("Usage: [inputfile] [k] [partitionsNum]")
    //   exit(1)
    // }
    // val inputFile = args(0)
    // val k = args(1).toInt
    // val partitionsNum = args(2).toInt

    // val name = "My Simple App " + partitionsNum  // .concat(num_partitions.toString)
    // val sc = new SparkContext(new SparkConf().setAppName(name))

    // val conf = new Configuration
    // conf.set("textinputformat.record.delimiter", "@")

   //  // starting with fastq format, first separate individual entries, then keep only read lines
   //  val readLines = 
   //  sc.newAPIHadoopFile(inputFile, classOf[TextInputFormat], classOf[LongWritable], classOf[Text], conf)
   //  .map(x => x._2.toString) // get rid of the @
   //  .filter(x => x.startsWith("SRR")) // make sure line isn't due to '@' inside the quality field
   //  .map(x => x.split("\n")(1)) // keep only the DNA sequence (the read), discard the name (header)

   //  // make sure the there is the right count and that lines look right
   //  val xx = readLines.count()   
   //  println(s"---------------- readLines.count: $xx")
   //  readLines.take(10).foreach(println)

   //  // get new RDD including lists of kmers (with no Ns), counts of next base extensions,
   //  // then reduceByKey to get counts per extension
   //  val kmers = readLines
   //  .flatMap(read=> getKmerToNextCharCounts(read, k))
   //  .reduceByKey((a,b) => (a,b).zipped.map(_ + _)) // element-wise sum
    
   //  val yy = kmers.count()
   //  println(s"---------------- kmers.count: $yy")
   //  kmers
   //  .take(100000)
   //  .filter{case (a,b) => b.sum>1}
   //  .foreach(x => println("#### " + x._1 +" "+ x._2.mkString(" "))) // only print if total count > 1

   //  // junctions - keep only kmers with multiple positive extensions
   //  val junctions = kmers.filter(tuple => filterTup(tuple._2)._2)
   // // val junctions = kmers.filter{case (str,arr) => arr.count(_>0) > 1}
   //  val zz = junctions.count()
   //  println(s"---------------- junctions.count: $zz")

   //  // broadcast so that each node gets all junctions instead of a slice
   //  val junctions_broadcast = sc.broadcast(junctions.collect().toSet) 
   //  val jj = junctions_broadcast.value.length
   //  println(s"---------------- junctions in broadcast: $jj")

    println("---------" + get_stranded_kmerized_minimizers_list("AACCTTGCTACGTCCAAG", 5, 3))    

    // sc.stop()
  }
}
