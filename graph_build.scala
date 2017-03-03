// graph building code - following example from:
// https://www.mapr.com/blog/how-get-started-using-apache-spark-graphx-scala

// *** I'm not sure all of these imports are needed
// import org.apache.spark._
// import org.apache.spark.rdd.RDD
// import org.apache.spark.util.IntParam
// // import classes required for using GraphX
// import org.apache.spark.graphx._
// import org.apache.spark.graphx.util.GraphGenerators


object SimpleApp {

    def rc(dna: String): String = {
        dna.reverse.map(b => b match {
        case 'A' => 'T'
        case 'T' => 'A'
        case 'C' => 'G'
        case 'G' => 'C'
        })
    }

    def get_edge_list(read: String, anchors_array: Array[(Int,String)], k: Int, l: Int) = {
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
            val edge_rep = edge_rep_strings(i)
            val left_rep = List(edge_rep.substring(0,k), rc(edge_rep.substring(0,k))).min
            val right_rep = List(edge_rep.substring(edge_rep.length - k, edge_rep.length), 
            rc(edge_rep.substring(edge_rep.length - k, edge_rep.length))).min
            val edge_tup = 
                (left_rep, right_rep,
                 (edge_rep.substring(0,k) == left_rep,
                 edge_rep.substring(edge_rep.length - k,edge_rep.length) == right_rep)
                )

        } yield edge_tup
        edge_list
    }
    def main(args: Array[String]) {

    // spaler scheme: (k+1)mer s is edge, rep(s) = min(s,s') 
    // polarity set relative to F/R of kmers relative to rep(s) [?]
    // e.g., if s = AAAT = (u,v); rep(s) = min(AAAT, ATTT) = AAAT; u = AAA, v = AAT
    // rep(u) = min(AAA,TTT) = AAA = u; pol(u) = false (didn't need to RC)
    // rep(v) = min(AAT, ATT) = AAT = v; pol(v) = false
    // spark edge: (AAA,AAT,(false,false))


    // our scheme:
    // 1) find junctions
    // 2) map each read to list of anchors - junctions + minimizers (extended to k)
    // output of 2) is should be an Array that looks like this:
    // Array((15, "ACGTTCACGGTTA"), (18, "ACGTTTTTTTTTA"), (35, "ATTTTCACAAAAA"))
    // indices are start positions of kmers obtained after extending lmer minimizers 
    // both (k-l)/2 bases in both directions
    val l = 7
    val k = 13
    // 2 edges pointing forward
    val read1 = "S"*15 + "ACGTTCACGGTTATTA" + "ATTTTCACAAAAA"
    val a1 = Array((15, "ACGTTCACGGTTA"), (18, "TTCACGGTTATTA"), (31, "ATTTTCACAAAAA"))
    get_edge_list(read1,a1,k,l)
    
    // 3 edges F,R,F - should see k-mers on middle edge in reverse order (according to direction we read the edge)
    val read2 = "S"*15 + "ACGTTCACGGTTATTA" + "AAATCGGTCGTTT" + "ATTTTCACAAAAA"
    val a2 = Array((15, "ACGTTCACGGTTA"), (18, "TTCACGGTTATTA"), (31, "AAATCGGTCGTTT"), (44, "ATTTTCACAAAAA"))
    get_edge_list(read2,a2,k,l)



    }
}