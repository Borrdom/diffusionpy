import pickledic
from tkinter.filedialog import askopenfilename
import os
from datetime import datetime
import matplotlib.pyplot as plt
import xlsxwriter
import shutil 
plt.tight_layout()
def cs_plt_excel(case_lst,unit_lst,safename_lst,cs_output_dic):
    # To create a chronological folder system for the casestudy, first a subfolder outside the folder of the python scripts needs to be generated.
    # To do so, the parent path has to be found and in the parent path a folder called 'casestudy' + the date and time is created.
    current_path = os.getcwd()
    parent_path = os.path.dirname(current_path)
    now = datetime.now()
    casestudy_name = now.strftime("casestudy#%d-%m-%Y_%H-%M-%S")
    casestudy_folder = os.path.join(parent_path, casestudy_name)
    os.makedirs(casestudy_folder)
    plot_font = {'fontname':'Calibri'}
    
    for index0, (k0,v0) in enumerate(cs_output_dic.items()):
        folder_path = os.path.join(casestudy_folder, k0)
        os.makedirs(folder_path)
        number_datasets = len(v0.items())
        for index1, (k1,v1) in enumerate(v0.items()):
            os.chdir(os.path.join(casestudy_folder, k0))
            row = 0
            col = 0
            if not isinstance(v1,str):
                excelfile = xlsxwriter.Workbook(os.path.join(folder_path,k1[0:]+'.xlsx'))
                exceltable = excelfile.add_worksheet('outputs')
                for index2, (k2,v2) in enumerate(v1.items()):
                    fig,ax = plt.subplots(figsize=(5, 5))
                    #
                    # print(v0.items())
                    for index3, (k3,v3) in enumerate(v0.items()):
                        if not isinstance(v3,str):
                            plt.rc('xtick',labelsize=26)
                            plt.rc('ytick',labelsize=26)
                            plt.ioff()                       
                            print(index2)
                            # ax.set_ylabel('$\\alpha$\n/$-$',r'$\mathregular{ '+case_lst[index2]+r' } \ / \;'+unit_lst[index2]+r'}$',**plot_font,size=26,loc='top',rotation=0)
                            ax.set_ylabel("$\/" + case_lst[index2] + "$\n/$" + unit_lst[index2] + "$",**plot_font,size=24,rotation=0)
                            ax.yaxis.set_label_coords(-0.01,1)
                            ax.set_xlabel(r'$\mathregular{time \ / \; h}$',**plot_font,size=24)
                            ax.plot(v3['t_vec'],v3[k2],'k-', linewidth=4, alpha=float(1/number_datasets*(index3+1)))
                            for axis in ['top','bottom','left','right']:
                                ax.spines[axis].set_linewidth(1.1)
                                
                            ax.tick_params(direction='in', right=True, top=True, labelright=False, labeltop=False, width=1.2, length=6)
                            
                            # plt.rcParams['xtick.major.pad']='12'
                            # ax.set_ylim(bottom=v3[k2][0])
                            # ax.set_xlim(left=0)
                            # ax.grid()
                            #plt.tight_layout()
                    savefile = k0+safename_lst[index2]+'.png'
                    fig.savefig(os.path.join(folder_path, savefile),bbox_inches='tight')
                    plt.close('all')
                    row = 0
                    exceltable.write(row,col,k2)
                    row += 1
                    for i in v2:
                        try:
                            exceltable.write(row,col,i)
                            row += 1
                        except:
                            exceltable.write(row,col,str(i))
                            row += 1
                    col += 1
                excelfile.close()    
            else:
                f= open(k1+".txt","w+")
                f.write(v1)
                f.close() 
    pickledic_path_origin = os.path.join(current_path,'cs_ouput_dic.pickle')
    pickledic_path_target = os.path.join(casestudy_folder,'cs_ouput_dic_'+casestudy_name+'.pickle')
    shutil.copyfile(pickledic_path_origin, pickledic_path_target)
    os.remove(pickledic_path_origin)
    par_excel_path_origin = os.path.join(current_path,'par_bounds.xlsx')
    par_excel_path_target = os.path.join(casestudy_folder,'par_bounds_'+casestudy_name+'.xlsx')
    shutil.copyfile(par_excel_path_origin, par_excel_path_target)
    
                
if __name__=="__main__":
    # filename=askopenfilename()
    current_path = os.getcwd()
    filename = os.path.join(current_path,'cs_ouput_dic.pickle')
    cs_output_dic=pickledic.pickledicload(filename)
    if len(cs_output_dic.items()) == 5:
        case_lst = ['t','N','l',r'\alpha','S',r'\Delta \mu',r'\dfrac{dN}{dt}',r'\dfrac{dr}{dt}','n^*','z','f^*_{dif}','f^*_{inter}']
        unit_lst = ['s','-',r'\mu m','-','-','-','s^{-1}',r'nm \cdot s^{-1}','-','-','s^{-1}','s^{-1}']
        safename_lst = ['t','N','r','alpha','s','dmu','nt','rt','nstar','z','fdif','finter']
    elif len(cs_output_dic.items()) == 6:
        case_lst = ['t','N','l',r'\alpha','S',r'\Delta \mu',r'\dfrac{dN}{dt}',r'\dfrac{dr}{dt}','n^*','z','f^*_{dif}','f^*_{inter}','w_w']
        unit_lst = ['s','-',r'\mu m','-','-','-','s^{-1}',r'nm \cdot s^{-1}','-','-','s^{-1}','s^{-1}','-']
        safename_lst = ['t','N','r','alpha','s','dmu','nt','rt','nstar','z','fdif','finter','ww']
    cs_plt_excel(case_lst,unit_lst,safename_lst,cs_output_dic)

# Three plots:
    # labelsize=28
    # font,size=26
    # linewidth(1.2)
    # width=1.3
    # length=7




