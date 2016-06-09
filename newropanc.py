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
    global mxlim,mnlim
    mxlim = 1e-8
    mnlim = 1e-13

    ############# split headers ######################
    rhead = [item.split('-->') for item in r_file.columns]
    reactants = [x[0].replace(' ','').split('+') for x in rhead[6:]] 
    products  = [x[1].replace(' ','').split('+') for x in rhead[6:]] 


    #################################################
    ####multiple timesteps? change runs #############

    graph_list, s_file_list = [],[]
    
    
    
    
    runs = [6] #### LAST RUN NEEDS TO BE ONE WITH GREATEST NO REACTIONS! ### 144,288,431,

    
    
    
    for timestep in runs:
        
        print '######################## timestep in days', (1+timestep)*10/60/24 
        #reset on each
        s_loop = (s_file.ix[timestep]/M).map(check_negative) #specie file for loop 
        r_loop = r_file.iloc[timestep] #reaction rate for loop 
        print 'converting to mixing ratio'
        
        
        ###### make dictionaries ######################### 
        #i_spec = {key: value for (key, value) in enumerate(s_loop.index)} #index of species
        #r_i_spec = {value: key for (key, value) in enumerate(s_loop.index)} #index of species
        concentrate = dict(zip( s_loop.index[9:] , s_loop.ix[9:]))

        ############# get lengths   ######################
        r_len= len(reactants)
        r_range= xrange(r_len)

      
        print 'begin time travel'
        #calculate total flux and retun locations
        edges = multiprocessing.Pool(ncores).map(flux_capacitor, np.array_split(range(r_len),ncores)) 
        total_flux =  multiprocessing.Pool(ncores).map(fluxify , np.array_split(range(r_len),ncores))
        
        

        
        

            
        dummy=[]    ;     [dummy.extend(i) for i in edges] #reformat 
        graph = pd.DataFrame(dummy,columns= ['Source', 'Target', 'weight']) #edge result to df
        graph['weight']=graph['weight'].astype(float)
        print 'individual flux calculated as the product of fluxes for each species'     

        del dummy 
        print 'you have arrived'
        

#######################################################
########################################################
        
        
        
        
        
        def adj_spec (x): 
            ''' eliminates outliers from the data'''
            if   ( x <= mnlim and x != 0) : x = mnlim
            elif   x > mxlim: x = mxlim
            return x 
     
            
        def adj_flux (x): 
            ''' eliminates outliers from the data'''
            if   ( x < mn-st and x != 0) : x = mn-st
            elif   x > mn+st: x = mn+st
            return x 
            
            
        def reverse(inp):
            empty,data=[],[]
            for item in inp: 
                    i1= re.sub(r'(.*) (.*)', r'\1', item)
                    i2= re.sub(r'(.*) (.*)', r'\2', item)
                    dummy = graph[ (graph['Target']==i1) & (graph['Source']==i2) ] #gets all duplicates
                    index = graph[ (graph['Source']==i1) & (graph['Target']==i2) ].index  #real value 
                    empty.extend(dummy.index) #reverse duplicates to delete
                    data.append([int(index[0]), float(abs(graph.loc[index,'weight']-dummy.sum()['weight']))])
            return [empty,data]
     
     
        def limit_by_smiles(what, s_cols):
                       
            smiles = pd.read_csv('MCMv3.2inchis.csv')
            inorganics= 'O,O3,O1D,O2,OH,NO,NO2,NO2N2O5,H2O2,HO2,HO2NO2,HONO,HNO3,CO,SO2,SO3,NA'.split(',')
            smiles = smiles[smiles.smiles.map(lambda x : bool(re.match('.*%s.*'%what,x)))]
            others=[]

            for item in inorganics: 
                if bool(re.match('.*'+what+'.*',item)): others.append(item)

            shead=list( set(list(smiles.name)+(others))  &  set(s_cols)  )  
            
                
            print 'taking only', what 
            
            return shead
            
        
       
        
        print 'move this to new module'   
        
        
        
        ### in function 
        ##require edges, s_loop, head, names,  ncores, total flux
        
        
        def datamonging_part ( s_loop, concentrate,  ncores, total_flux, onlyCarbon=True, onlyNitrogen=False):
         
        #if True:        
            ######## calc total flux #######     
            print 'Getting edges...'
            global st,mn,shead,graph

            #brute force and ignorance here. 
      
            #slow
            graph = graph[graph['Target'] != graph['Source']] #no self reactions
            graph = graph.groupby(['Source','Target']).agg(sum).reset_index()#remove duplicates

            '''
            Since on a force graph if you have two springs of different strengths, they apply a net force. 
            We are concerned about the system shape stronger forces tend to take priority. 
            In normal force graphs this ends up added together, this cannot be done when treating the springs as a semi-static edge. 
            Therefore the net calculation has to be applied here. If further details on the reaction are required, seek further diagnostic tools, not just the force graph iteself. 
            '''  

            duplicates=set(graph.Source+' '+graph.Target) & set(graph.Target+' '+graph.Source)#overlap reversed 1140, 2583, 2696, 4341
            
            delete = multiprocessing.Pool(ncores).map(reverse, np.array_split(list(duplicates),ncores))
            
            empty=[]
            for item in delete: 
                empty.extend(item[0])
                for update in item[1]: graph.loc[update[0],'weight'] = update[1]

            graph.drop(graph.index[empty],inplace=True)
            graph = graph[graph['Target'] !='DUMMY']#rm depos reactions

            del duplicates, delete, empty
            print 'duplicates and net reactions netted'

            ###########################################

            if onlyCarbon : shead = limit_by_smiles('C',s_loop.index)
            elif onlyNitrogen : shead = limit_by_smiles('N',s_loop.index)
            
            s_loop[shead]
            graph = graph[graph['Target'].isin(shead) & graph['Source'].isin(shead)]  
            graph = graph.reset_index(drop= True)
            
            ## normalise flux strength
            

            print 'normalising flux'

            #.to_csv('gephi_input.csv',index=False)

            ###split reactions with data
            

            reactions= graph[['Source','Target']]
            # make seperate 
            # graph.drop(['Source','Target'], axis=1, inplace=True)

            
            
            #### filter outliers
            s_loop = s_loop.map(adj_spec)  # remove null results and adjust max / minimums
            s_loop = s_loop[s_loop>0]
            
            flux = graph.weight.map(lambda x : np.log10(x))
            st, mn= 2*flux.std(), flux.mean()
            flux = flux.map(adj_flux) 
            flux = flux+ abs(min(flux)) + 1
            graph['flux'] = (flux-flux.min())/flux.max()
            
            #total flux - kludge
            for i in total_flux[1:]:total_flux[0] = total_flux[0] + i 
            total_flux = total_flux[0]
            total_flux = total_flux / np.array( [0]*9+ list(np.vectorize(concentrate.__getitem__)(np.array(total_flux.index[9:])) ) )
            total_flux.sort()
            pt_flux =  total_flux[total_flux!=0] > 0 # a booian of which total fluxes are positive
            total_flux = total_flux[total_flux!=0].map(lambda x: abs(np.log10(abs(x))))
            total_flux= (total_flux-total_flux.min())/total_flux.max()

            
            s_loop= pd.concat([s_loop,total_flux,pt_flux],axis=1).sort(1)
            s_loop.columns = ['conc','tflux','netdir']
            return [  s_loop ,  graph ]


outlist = datamonging_part ( s_loop, concentrate,  ncores, total_flux, onlyCarbon , onlyNitrogen)




################################################################
###############################################################






print 'string bar'













