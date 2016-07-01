import numpy as np
import pandas as pd
import sys,os,re,multiprocessing,netCDF4
from netCDF4 import Dataset
import matplotlib.pyplot as plt
#netcdf file path
ncfile = sys.argv[1]
ncores = 16  
global graph, timestep, s_loop
global reactants, products,r_loop,concentrate,getspec, getcoeff 




#switches
onlyCarbon = True
onlyNitrogen = False
force=True
#maxlines = os.system('wc -l %s'%rfile) - 1  #get file lenght (number of timesteps)



######################################################################################## 
###########################################################################  Functions 

def getcoef (x): 
    '''gets specie coefficients from data'''    
    a = re.sub(r'([\.\d]*)\s*\D[\d\D]*', r'\1', x)
    try: return int(a) 
    except: return 1  #assume coeff are Z+      

getspec = lambda x : re.sub(r'([\.\d\s]*)(\D[\d\D]*)', r'\2', x)


def flux_capacitor(item):#  gets flux & coefficients 
    global reactants,products,r_loop,concentrate,getspec,getcoeff
    e = []
    for i in item:
        r_coeff      =  np.vectorize( getcoef )( reactants[i]) #get react coeff
        reactants[i] =  np.vectorize( getspec )( reactants[i])  #rm react coeff
        products [i] =  np.vectorize( getspec )( products [i])  #rm prod coeff    
        r_conc       =  np.vectorize( concentrate.__getitem__)(reactants[i])#get conc @ Timestp
        t_flux       =  np.prod( r_coeff * r_conc ) * r_loop[i]
            
        if t_flux>0: ######## get edge reactions #######
            for r_spec in reactants[i]: 
                for p_spec in products[i]:
                    if (len(p_spec)) > 0 and (len(r_spec)>1): 
                        e.append( [r_spec,p_spec,float(t_flux)] )

            del r_conc, t_flux, r_coeff
    return e




def fluxify (item):
        shead = s_loop.index
        total_flux = pd.Series(np.zeros(len(shead)),index=shead)    
        global reactants,products,r_loop,concentrate,getspec,getcoeff
        e = []
        for i in item:
            r_coeff      =  np.vectorize( getcoef )( reactants[i]) #get react coeff
            reactants[i] =  np.vectorize( getspec )( reactants[i])  #rm react coeff
            products [i] =  np.vectorize( getspec )( products [i])  #rm prod coeff    
            r_conc       =  np.vectorize( concentrate.__getitem__)(reactants[i])#get conc @ Timestp
            t_flux       =  np.prod( r_coeff * r_conc ) * r_loop[i]
            
            if t_flux>0: ######## get edge reactions #######
                for r_spec in reactants[i]: 
                    count = False
                    for p_spec in products[i]:
                        if (len(p_spec)) > 0 and (len(r_spec)>1): 
                            total_flux[p_spec] = total_flux[p_spec] + float(t_flux)
                            if (not count): 
                                total_flux[r_spec] = total_flux[r_spec] - float(t_flux)
                                count = True

            return total_flux      
    
    
def check_negative(x):
    if x<0: x=0 
    return x    
######################################################################################## 
################################################################################## start 


netCDF_data = Dataset(ncfile, mode='r')

for group in netCDF_data.groups:
    print '\n\nReading', group

    #data
    s_file = pd.DataFrame(netCDF_data.groups[group].variables['Spec'][:])
    r_file = pd.DataFrame(netCDF_data.groups[group].variables['Rate'][:])

    #headers
    r_file.columns = str(netCDF_data.groups[group].variables['Rate'].head).split(',')
    s_file.columns = str(netCDF_data.groups[group].variables['Spec'].head).split(',')

    #get m value 
    M = s_file['M'].mean()


    ############# split headers ######################
    rhead = [item.split('-->') for item in r_file.columns]
    reactants = [x[0].replace(' ','').split('+') for x in rhead[6:]] 
    products  = [x[1].replace(' ','').split('+') for x in rhead[6:]] 


    #################################################
    ####multiple timesteps? change runs #############

    graph_list, s_file_list = [],[]
    
    
    
    
    runs = [24*6]# range(2) #### LAST RUN NEEDS TO BE ONE WITH GREATEST NO REACTIONS! ### 144,288,431,

    
    
    
    for timestep in runs:
        
        print '######################## timestep in days', (1+timestep)*10/60/24 
        #reset on each
        s_loop = (s_file.ix[timestep]).map(check_negative) #specie file for loop 
        r_loop = r_file.iloc[timestep] #reaction rate for loop 
        print 'converting to mixing ratio'
        
        
        ###### make dictionaries ######################### 
        #i_spec = {key: value for (key, value) in enumerate(s_loop.index)} #index of species
        #r_i_spec = {value: key for (key, value) in enumerate(s_loop.index)} #index of species
        concentrate = dict(zip( s_loop.index[9:] , s_loop.iloc[9:]))

        ############# get lengths   ######################
        r_len= len(reactants)
        r_range= xrange(r_len)

      
        print 'begin time travel'
        #calculate total flux and retun locations
        edges = multiprocessing.Pool(ncores).map(flux_capacitor, np.array_split(range(r_len),ncores)) 
        total_flux =  multiprocessing.Pool(ncores).map(fluxify , np.array_split(range(r_len),ncores))
        s_loop = s_loop/M # we can now convert concentrations into a normal type
           
        dummy=[]    ;     [dummy.extend(i) for i in edges] #reformat 
        graph = pd.DataFrame(dummy,columns= ['Source', 'Target', 'weight']) #edge result to df
        graph['weight']=graph['weight'].astype(float)
        print 'individual flux calculated as the product of fluxes for each species'     

        del dummy 
        print 'you have arrived'
        
        

        
        
        if force :
            import ropa_to_graph
            ropa_to_graph.graph = graph
            outlist = ropa_to_graph.datamonging_part (s_loop, concentrate,  ncores, total_flux, onlyCarbon , onlyNitrogen)
            del ropa_to_graph
            
        #dictionary for enumeration       
        specdata = pd.read_csv('fullmcmspecs.csv')
        global ids
        ids =  dict(zip(specdata['item'], specdata['id']))
          
 
        #random
        logs= lambda x : x.map(lambda x : np.log(x))
        
      
        spc = outlist[0].fillna(0)#[6:]#.dropna(axis=0,how='all', subset=['conc'])
        nodes = set(outlist[1]['Target']) | set(outlist[1]['Source'])#only reactions in read list
        spc = spc.ix[list(nodes)]
        

        def spc_id(x):
            try: return  ids[x]
            except: 
                print 'failed on', x
                return -1 #error flag..
        
        #convert names into numeric form
        spc['id'] = np.vectorize(spc_id)( np.array(spc.index)) 
        outlist[1][['Source','Target']] = np.vectorize(lambda x: ids[x])( np.array(outlist[1][['Source','Target']])) 
        gh = outlist[1]
        

        
        
        string = r'{"nodes":[ '
        
        for node in spc.iterrows(): 
            string +=  '{"name":"%s","x":%.2f,"y":%.2f,"shape":"circle","s":%.2f,"id":%d},'%(node[0],np.random.random(),np.random.random(),node[1].size,node[1].id)

        string = string[:-1] 
        string += r'],"links":['        
        

        for i in gh.iterrows(): 
            string += '{"source":%d,"target":%d,"value":%.2f},' %(i[1].Source, i[1].Target, i[1].flux)


        string = string[:-1]
        string += r']}'

        filename = './JSON_Files/%s_%d.json'%(group[:3],timestep)
        f = open(filename,'w')
        print 'add date here'
        f.write(string)
        f.close()

        os.system('scp ./%s dp626@research0.york.ac.uk:/usr/userfs/d/dp626/web/force/ics'%filename)

