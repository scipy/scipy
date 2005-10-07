from shelve import Shelf

class DbfilenameShelf(Shelf):
	"""Shelf implementation using the "anydbm" generic dbm interface.

	This is initialized with the filename for the dbm database.
	See the module's __doc__ string for an overview of the interface.
	"""
	
	def __init__(self, filename, flag='c'):
		import dumbdbm_patched
		Shelf.__init__(self, dumbdbm_patched.open(filename, flag))


def open(filename, flag='c'):
	"""Open a persistent dictionary for reading and writing.

	Argument is the filename for the dbm database.
	See the module's __doc__ string for an overview of the interface.
	"""
	
	return DbfilenameShelf(filename, flag)
