import numpy as np
import pandas as pd
import sys,os,re,multiprocessing


global graph, mxlim,mnlim
mxlim = 1e-8
mnlim = 1e-13



def adj_spec (x): 
    ''' eliminates outliers from the data'''
    if   ( x <= mnlim and x != 0) : x = mnlim
    elif   x > mxlim: x = mxlim
    return x 

    
def adj_flux (x): 
    ''' eliminates outliers from the data'''
    print x, mn ,st
    if   ( x < mn-1*st and x != 0) : x = mn-1*st
    if   ( x < mn-(1.5*st) and x != 0) : x = 0
    #elif   x > mn+st: x = mn+st # dont cap upper
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
#limit types using smile strings 
    smiles = pd.read_csv('MCMv3.2inchis.csv')
    inorganics= 'O,O3,O1D,O2,OH,NO,NO2,NO2N2O5,H2O2,HO2,HO2NO2,HONO,HNO3,CO,SO2,SO3,NA'.split(',')
    smiles = smiles[smiles.smiles.map(lambda x : bool(re.match('.*%s.*'%what,x)))]
    others=[]
    for item in inorganics: 
        if bool(re.match('.*'+what+'.*',item)): others.append(item)
    shead=list( set(list(smiles.name)+(others))  &  set(s_cols)  )         
    print 'taking only', what 
    return shead
    
##############################################################################################

def datamonging_part ( s_loop, concentrate,  ncores, total_flux, onlyCarbon=True, onlyNitrogen=False):
 

    ######## calc total flux #######     
    print 'Getting edges...'
    global st,mn,shead,graph
    

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

    print 'normalising flux'

    ###split reactions with data
    reactions= graph[['Source','Target']]
    # make seperate 

    
    #### filter outliers
    s_loop = s_loop.map(adj_spec)  # remove null results and adjust max / minimums
    s_loop = s_loop[s_loop>0]
    
    
    
    
    flux = graph.weight.map(lambda x : np.log10(x))
    st, mn=1*flux.std(), flux.mean()
    print '2 - 1.5 std at 1.5 DELETE OVER 2STD'
    flux = flux.map(adj_flux) 
    flux = flux + abs(min(flux)) + 1
    graph['flux'] = (flux-flux.min())/flux.max()
    
    
    
    #total flux - kludge
    for i in total_flux[1:]:total_flux[0] = total_flux[0] + i 
    total_flux = total_flux[0]
    total_flux = total_flux / np.array( [0]*9+ list(np.vectorize(concentrate.__getitem__)(np.array(total_flux.index[9:])) ) )
    total_flux.sort()
    pt_flux =  total_flux[total_flux!=0] > 0 # a booian of which total fluxes are positive
    total_flux = total_flux[total_flux!=0].map(lambda x: abs(np.log10(abs(x))))
    total_flux= (total_flux-total_flux.min())/total_flux.max()+0.01
    
    size = s_loop.map(lambda x: np.log10(x)).ix[:]
    size = ((size+abs(min(size)))/abs(max(size)))+0.01
    s_loop= pd.concat([s_loop, size, total_flux, pt_flux],axis=1).sort(1)
    s_loop.columns = ['conc', 'size' ,'tflux','netdir']
    return [  s_loop ,  graph ]
    
    
    
