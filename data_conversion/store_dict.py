import sys

version = float(sys.version[:3])

if version >= 3.6 and version < 3.8:
    import pickle5 as pickle
else:
    import pickle

if version >= 3.6:
    prot = 5
elif version >= 3:
    prot = 4
else:
    prot = 2

def save_obj(obj, name):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, prot)

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)
