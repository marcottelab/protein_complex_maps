
import validate_by_interactions as vbi
import MySQLdb

db = MySQLdb.connect("localhost", 'kdrew', 'kdrew_utexas', 'genemania_20131024')


iset = vbi.InteractionSet(db)
iset.get_total_weight(['O75832', 'Q06323', 'Q9UL46'])

iset.get_interactions_from_db(['O75832', 'Q06323', 'Q9UL46'])
iset.get_interactions_from_db(['O75832', 'Q06323', 'Q9UL46'])

total_weight = iset.get_total_weight(['Q06323', 'Q9UL46'])
print total_weight
assert( total_weight == 6.98 )

print iset.interactions
