import numpy as np
from pebble import concurrent #Module needed
from concurrent.futures import TimeoutError
import importlib
from inputmatrix_auto import get_matrix_input

def ImportModel(sname,fname):
    ModelScript = importlib.import_module(sname)
    return getattr(ModelScript,fname)
def ImportInputs(listnames,listarrays):
    return dict(zip(listnames,listarrays))
def ConstructCaseSpace(*args):
    Grid=np.meshgrid(*args)
    return [val.flatten() for val in Grid]
     
################### Needs to be adjusted
NameOfModelScript="cryst_real"
NameOfFunctionInModelScript="cryst"
# ListOfNamesOfInputs=['rho', 'M', 'temp','sigma', 'kt', 'g', 'D','Dw']
ListOfNamesOfInputs=['rho', 'M', 'temp','sigma', 'kt', 'g', 'D']
timeout=10 
######################
# MatrixOfInputs_dict = get_matrix_input()
# for index,val in MatrixOfInputs_dict.items():
#     MatrixOfInputs = val
Model=ImportModel(NameOfModelScript,NameOfFunctionInModelScript)
# input_dic=ImportInputs(ListOfNamesOfInputs,ListOfArraysOfInputs)


@concurrent.process(timeout=timeout) #Function decorator for Pebble 
def Model_TimeOut(*inputs):
    return Model(*inputs)
                
def casestudy(model,inputnames,inputmatrix):
    
    output_dic={}
    for i,val in enumerate(inputmatrix[0,:]):
        inputs=tuple(inputmatrix[:,i])
        resultstring="_".join([inputnames[j]+"="+str(inputs[j]) for j,valj in enumerate(inputnames)])
        try:

            fun_output=model(*inputs).result() #
                    
        except TimeoutError as error:
            fun_output="Function call for inputs "+resultstring+" took longer than %d seconds" % timeout
            print(fun_output)

        except Exception as error:
            fun_output="Function raised %s" % error
            print("Function raised %s" % error)
            
        output_dic[resultstring]=fun_output       
    return output_dic
        
if __name__=="__main__":
    from inputmatrix_auto import get_matrix_input
    MatrixOfInputs_dict = get_matrix_input()
    cs_output_dic = {}
    for key,val in MatrixOfInputs_dict.items():
        MatrixOfInputs = val
        cs_output_dic.update({key:casestudy(Model_TimeOut,ListOfNamesOfInputs,MatrixOfInputs)})
    import pickledic
    import os
    cwd=os.getcwd() # 
    filename=os.path.join(cwd,'cs_ouput_dic.pickle')
    pickledic.pickledicsave(cs_output_dic,filename)

