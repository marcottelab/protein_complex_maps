# I make an independent file for Struct so that when reloading utils, I don't
# create a different Struct class (re.: bloody Python issues)
# -Cedric

class Struct:
    # Norvig's anonymous class
    # http://norvig.com/python-iaq.html
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        args = ['%s=%s' % (k, repr(v)) for (k,v) in vars(self).items()]
        return 'Struct(%s)' % ', '.join(args)

class sStruct:
    # A silent Struct that does not shout its objects whenever it's mentionned
    # (to prevent REPL nightmares if the objects are very long)
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        if hasattr(self, 'dumb_type'):
            return 'sStruct(%s)' % self.dumb_type
