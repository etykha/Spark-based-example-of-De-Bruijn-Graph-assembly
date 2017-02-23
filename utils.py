# ⓒ Copyright IBM Corp. 2017

# import process_reads as pr
from operator import itemgetter


def get_pos_kmer_tups(seq, k):
    return [(i,seq[i:i+k]) for i in range(len(seq)-k+1)]  

def get_minimizer_w_pos(kmer_tup, l):
    # this gets the minimizer and its position in the window
    assert  l <= len(kmer_tup[1])
    return min(get_pos_kmer_tups(kmer_tup[1], l),key=itemgetter(1))
    
def get_minimizers_list(seq,k,l):
	""" a bit convoluted, but unavoidable - 
		the intuition is to get all minimizers out of a read
		including their start positions with just a single pass on the read
		gets the minimizer for the initial window, then changes minimizer
		only if either a smaller minimizer enters the window, or the current
		minimizer drops out. Returns a list of tuples where item 1 in each tuple
		is the position of the minimizer, and item 2 is the minimizer sequence
	"""

	kmers = get_pos_kmer_tups(seq,k)
	res = []
	curr_min = 'S'* l # chosen b/c DNA max string is 'TT...T'
	window_pos = -1
	for kmer in kmers:
		# print kmer
		if (curr_min > kmer[1][-l:]) or (window_pos < 0):
			window_pos, curr_min = get_minimizer_w_pos(kmer, l)
			# append kmer start position instead of minimizer start position
			res.append((window_pos + kmer[0], curr_min))
			window_pos -= 1
			# print curr_min, "added"
		else:
			# print curr_min, window_pos
			window_pos -= 1
	return res

def get_kmerized_minimizers_list(seq,k,l):
	""" as above, but with the requirement that it be possible
		to extend minimizers to k characters
	"""
	kmers = get_pos_kmer_tups(seq,k)
	res = []
	curr_min = 'S'* l
	window_pos = -1
	for kmer in kmers:
		# print kmer
		if (curr_min > kmer[1][-l:]) or (window_pos < 0):
			window_pos, curr_min = get_minimizer_w_pos(kmer, l)
			# only include minimizer if it can be extended to a kmer
			minimizer_pos = window_pos + kmer[0]
			if (minimizer_pos > len(seq) - k + 1): break

			# append kmer start position instead of minimizer start position
			res.append((minimizer_pos, seq[minimizer_pos: minimizer_pos + k]))
			window_pos -= 1
			# print curr_min, "added"
		else:
			# print curr_min, window_pos
			window_pos -= 1
	return res

def map_read_to_anchors_list(seq, k, l, juncs_set):
	""" given seq, break it up into a list of all of its anchor
	    (junction or minimizer filled up to k characters) points;
	"""
	from collections import deque
	# use deque to make popping front faster
	mins_lst = deque(get_kmerized_minimizers_list(seq,k,l)) 
	allKmers = get_pos_kmer_tups(seq,k)
	read_anchors = []
	for kmer in allKmers:
	    if kmer[1] in juncs_set:
	        read_anchors.append(kmer)
	    elif kmer in mins_lst:
	        read_anchors.append(kmer)
	        mins_lst.popleft()
	return read_anchors

def convert_anchors_list_to_seq_edges(anchors, k=31):
	""" Given an anchors list where entries are (a,b) tuples 
		where a values are seq positions and b values are anchors
		(junctions or minimizers extened to kmers) of the seq,
		converts to list of edges where end-points are consecutive anchors
		and a third entry is a value implying the sequence between the anchors.
		If there is a gap between anchors, the value is the missing sequence as a
		string. If there is an overlap, the value is the length of overlap. If 
		the anchors are directly adjacent (i.e., the last character of the first
		is followed by the first character of the second), the value is 0.	
	"""
	edges_lst = []
	for pos in range(len(anchors)-1):
		start_anc = anchors[pos]
		end_anc = anchors[pos+1]
		overlap_len = -(end_anc[0] - (start_anc[0] + k - 1) - 1)
		# I think there is always either an overlap or direct
		# adjacency, otherwise part of seq isn't covered by a minimizer
		assert(overlap_len >= 0)
		# if overlap_len >= 0:
		edge_val = overlap_len
		
		edges_lst.append((start_anc[1], end_anc[1], edge_val))
	return edges_lst


if __name__ == "__main__":

	read = "231032101233101"
	k=6; l=3
	juncs = set(['32101', '10123'])
	# print read
	# print faster_get_minimizers_list(read,k,l)
	assert(get_minimizers_list(read,5,3) == 
		[(2, '103'), (3, '032'), (6, '101'), (7, '012'), (8, '123'), (9, '233'), (12, '101')])

	assert(get_kmerized_minimizers_list(read,5,3) == 
		[(2, '10321'), (3, '03210'), (6, '10123'), (7, '01233'), (8, '12331'), (9, '23310')])

	print map_read_to_anchors_list(read, k, l, juncs)
	anchors = map_read_to_anchors_list(read, k, l, juncs)
	print convert_anchors_list_to_seq_edges(read, k, anchors)