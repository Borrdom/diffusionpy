#run in console before running pythoncode the
#git clone https://github.com/simonmb/fragmentation_algorithm
#run python -m pip install rdkit-pypi 

def RunFrag(smiles):
    from fragmentation_algorithm.fragmenter import fragmenter #Better practice to import it this way. It import the class fragmenter from the folder ../fragmentation_algorithm.fragmenter//fragmenter.py
    
    fragmentation_scheme = {
        'CH2' : '[CH2]',
        'OH' : '[OH]',
        'CH3' : '[CH3]',
        'CH2-CH2' : '[CH2][CH2]'
    } # data type with {} signalizes a dictionary which is similar to the matlab struct in matlab 
    #I.E you can try
    print(fragmentation_scheme["CH2"])
    
    # FYI: List are different from matlab vectors in that they are not designed to be enable big calculations in vector or matrix style (for that use numpy instead)
    # Instead they serve as the most very general form of data storage and can even store classes,function etc. 
    fragmentation_scheme_order1 = ['CH2-CH2', 'CH3', 'CH2', 'OH']
    
    # fragmenter is the class. But we need to call an instance of the class
    # Every class has a __init__ function. If you look in the fragmenter.py file you can see the same inputs that are suplied in the example
    # After the __init__ of a class is called an instance of that class is generated here called frg
    frg = fragmenter(fragmentation_scheme, fragmentation_scheme_order=fragmentation_scheme_order1, algorithm='simple')

    # The following line put all available methods of this class instance frg in a list. Good if you want to learn what a class instance can do.
    # a method is like a function. In this regard a class instance is just a big pile of functions
    method_list = [func for func in dir(frg) if callable(getattr(frg, func)) ]
    print(method_list) #Here we can see the fragment method used in the next lines
    

    #Next the method fragment is called. A method from the frg class instance is called by the dot. The rest is like a function call
    fragmentation, success, fragmentation_matches = frg.fragment(smiles)
    
    return smiles, fragmentation, success,fragmentation_matches


smiles = ['CCCCO', 'CCCO', 'CCO', 'CO']
# The example requires one string as input at a time. Therefore the loop here
results=[RunFrag(val) for val in smiles]
# Now you have a list of results. In this list are tuples () which follow the order from the return statement of your function "RunFrag" => smiles, fragmentation, success,fragmentation_matches
result1=results[0]
# A tuple can be easily extrated like this
smiles1, fragmentation1, success1,fragmentation_matches1=result1
print(smiles1)
print(fragmentation1)
print(success1)
print(fragmentation_matches1)
# The result from this is a dictionary e.g fragmentation_matches1. If you dont know the "fields" of this "struct" you can use the "keys" and values methods
#  You will notice that ALMOST EVERYTHING in python is a class as "keys" is also a method from the class instance of a dictionary
print(fragmentation_matches1.keys())
print(fragmentation_matches1.values())



#________________________________

