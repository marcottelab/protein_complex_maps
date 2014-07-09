
class ClusterOutOfBounds(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)


def get_cluster(Y, cluster_id):

	if cluster_id > 2*len(Y):
		raise ClusterOutOfBounds(cluster_id)

	#kdrew: leaf node
	if cluster_id < (len(Y)+1):
		return [cluster_id]

	#kdrew: Y is of size n-1, the first n clusters are leafs (original observations)
	Y_id = cluster_id - (len(Y) + 1)

	original_ids = get_cluster(Y, Y[Y_id][0])
	original_ids += get_cluster(Y, Y[Y_id][1])

	return original_ids


