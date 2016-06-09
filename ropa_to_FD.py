    #######################################
        #datamonging part        #clean results 
        #######################################


        graph = graph[graph['Target'] != graph['Source']] #no self reactions
        graph = graph.groupby(['Source','Target']).agg(sum).reset_index()#remove duplicates


        copy=graph

        '''
        Since on a force graph if you have two springs of different strengths, they apply a net force. Since we are concerned about the system shape stronger forces tend to take priority. In normal force graphs this ends up added together, this cannot be done when treating the springs as a semi-static edge. Therefore the net calculation has to be applied here. If further details on the reaction are required, seek further diagnostic tools, not just the force graph iteself. 
        '''    

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

        duplicates=set(graph.Source+' '+graph.Target) & set(graph.Target+' '+graph.Source)#overlap reversed 1140, 2583, 2696, 4341
        delete = multiprocessing.Pool(ncores).map(reverse, np.array_split(list(duplicates),ncores) )
        empty=[]
        for item in delete: 
            empty.extend(item[0])
            for update in item[1]:
                graph.loc[update[0],'weight'] = update[1]

        graph.drop(graph.index[empty],inplace=True)
        graph = graph[graph['Target'] !='DUMMY']#rm depos reactions

        del duplicates, delete, empty
        print 'duplicates and net reactions netted'


        #######################################################################################
        #### To limit species, one must first limit reactions top only those species specified"
        full_results, main_shead = graph, shead  

        def limit_by_smiles(what):
            global graph, inorganics, shead, smiles, timestep,s_loop

            smiles = smiles[smiles.smiles.map(lambda x : bool(re.match('.*C.*',x)))]
            others=[]

            for item in inorganics: 
                if bool(re.match('.*'+what+'.*',item)): others.append(item)

            shead=list(set(list(smiles.name)+(others)) & set(s_loop.columns))  
            s_loop = s_loop.ix[timestep][shead]

            graph = graph[graph['Target'].isin(shead) & graph['Source'].isin(shead)]  
            graph = graph.reset_index(drop= True)    
            print 'taking only', what 


        #######################################################
        if onlyCarbon : limit_by_smiles('C')
        elif onlyNitrogen : limit_by_smiles('N')
        else: s_loop = s_loop.ix[timestep][shead]

        ## normalise flux strength
        graph['flux']= graph.weight
        graph['weight'] = np.log10(graph.weight) + 0.001
        graph['weight'] = graph.weight + abs(graph.weight.min()) + 0.001
        graph['weight'] = 0.1 +0.9*( graph.weight / graph.weight.max() )

        print 'log10(flux) + min / max * .9 + .1'

        #.to_csv('gephi_input.csv',index=False)

        ###split reactions with data
        reactions= graph[['Source','Target']]
        graph.drop(['Source','Target'], axis=1, inplace=True)


        graph.columns = graph.columns+ str(timestep)
        graph = pd.concat([reactions,graph],axis=1)
        graph_list.append(graph)
        s_file_list.append(s_loop)

    #################################################s



      #print make reacts with table


    if force:   

        def conect(inputs):
            empty=[]
            for it in inputs:
                dummy = []
                dummy.extend( list(graph[graph['Source']==it]['Target']) )
                dummy.extend( list(graph[graph['Target']==it]['Source']) )
                empty.append( [dummy, len(dummy)] ) 
            return empty

        connections, graphdata=[],[];

        for item in multiprocessing.Pool(ncores).map(conect, np.array_split(shead,ncores)): connections.extend(item)
        reactswith = dict(zip(shead,connections))

        ##concat new
        graph  = reactions # pick last point as should have most reactions. 

        def concat_g (g_select):
            rowval=[]
            for i in g_select:
                row=[i[0],i[1]]
                for gph in graph_list:
                    reaction = np.array( gph[(gph['Source']==i[0]) & (gph['Target']== i[1])] )
                    try:  row.extend( reaction[0][-2:]) 
                    except:    row.extend([])
                rowval.append(row)     
            return rowval

        for item in multiprocessing.Pool(ncores).map(concat_g, np.array_split(np.array(graph),ncores)): graphdata.extend(item)
        graph = pd.DataFrame(graphdata)


        col = ['Source','Target']
        for run in reversed(runs): col.extend(['weight%s'%run ,'flux%s'%run])    
        graph.columns = col ;  graph = graph.fillna(0)
        s_file = pd.concat(s_file_list, axis=1)


        if onlyCarbon : filename='C_%s.json'%group[:3].upper()
        elif onlyNitrogen : filename= 'N_%s.json'%group[:3].upper()
        else: filename='ALL_%s.json'%group[:3].upper()

        f = open(filename,'w')
        ##input data
        def inform (inputs):
            global smiles, inorganics 
            empty=[]
            for value in inputs:
                if value in inorganics: empty.append('Name: '+value+'\\nInorganic')
                else:

                    info= smiles[smiles['name'] == value]
                    empty.append(re.sub(r"[\[\'\]]", r'', ' Name: %s \\n Smiles: %s \\n %s \\n Mass: %s '%(value, list(info.smiles), list(info.inchi), list(info.mass)) ))
            return empty


        ##carbon oxygen ratio  
        def ocratio(inputs):
            global smiles, inorganics
            empty=[]
            for value in inputs:
                if value in inorganics: empty.append(value.count('O')/value.count('C'))
                else:
                    value= str(smiles[smiles['name'] == value].smiles)
                    empty.append(value.count('O')/value.count('C'))
            return empty
        ##number of nitrogens
        def nN(inputs):
            global smiles, inorganics
            empty=[]
            for value in inputs:
                if value in inorganics: empty.append(value.count('N'))
                else: empty.append(str(smiles[smiles['name'] == value].smiles).count('N'))
            return empty       

        print '\nmaking json file'




        p=multiprocessing.Pool(ncores)

        string, size = [], []
        for item in p.map(inform, np.array_split(list(s_loop.index),ncores)): string.extend(item)

        if onlyCarbon: 
            for item in p.map(ocratio, np.array_split(list(s_loop.index),ncores)): size.extend(item)
        if onlyNitrogen: 
            for item in p.map(nN, np.array_split(list(s_loop.index),ncores)): size.extend(item)

        print 'replacing inf, addingpriary remoe empties from here'

        primary = [i.split(',')[1] for i in netCDF_data.initial_conditions_str.split('\n')[7:-1]]


        for run in runs:

            if total_flux: 
                colnm = 'abflx'+str(run)
                fxnm = 'flux'+str(run)
                s_file[colnm]=0
                for i in graph.iterrows():  
                    s_file.loc[i[1]['Source'],colnm] -= float(i[1][fxnm])
                    s_file.loc[i[1]['Target'],colnm] += float(i[1][fxnm])
                #flux in is +ve out is -ve (so souce is -ve) 

            print 'total flux'



            s_file['conc'+str(run)]=s_file[run]/M
            s_file['size'+str(run)] = np.log10(s_file['conc'+str(run)])
            s_file['size'+str(run)].replace([np.inf, -np.inf], 0,inplace=True)
            s_file['size'+str(run)] = s_file['size'+str(run)] + abs(s_file['size'+str(run)].min()) 
            s_file['size'+str(run)] = 60*(s_file['size'+str(run)] / s_file['size'+str(run)].max())
            #  /s_file['conc'+str(run)].max()
##############        
        s_file.drop(runs,1,inplace=True)

        s_file['Info']=string

        s_file['name']=s_file.index
        #if size == []: size=s_loop['Concentration']
        #s_file['oc']=size    
        s_file['reactions'] = s_file['name'].map(reactswith.__getitem__)
        s_file=s_file.replace([np.inf, -np.inf], 0)
        #s_file=s_file[s_file.sum(axis=1)<1e90]


        print 'writing string'
        nind = dict([[v,k] for k,v in enumerate(shead)]) 
        ####nodes
        string ='{ "runs":'+ str(runs)+', "nodes":['
        #circle params 
        counter = 0 
        xcen = 300
        ycen = 300
        rad = 100;
        numberOfPoints = len(s_file)
        theta = 360/numberOfPoints;


        print ' adding locations -counter'
        for item in xrange(numberOfPoints):
            item = s_file.ix[item]
            counter+= 1
            size=''
            for run in runs:
                size+= '"node_size%s":%e,' %(run, item['size%d'%run]+0.1            )

            if (item['name'] in primary):  shape = 'triangle-up'
            else:  shape = 'circle'


            #get alkene alkane, aromatic 
            dummy = smiles[smiles.name == item['name']].smiles
            if (len(dummy)>0) and re.match('.*(c).*',str(dummy).split('\n')[0]):     
                it =6 # aromatic
            elif (len(dummy)>0) and re.search('^(?![ABD-Z]).*',str(dummy).split('\n')[0]): 
                it = 1 #'alkane' 
            elif (len(dummy)>0) and re.match('.*(C[()]{0,1}\=[()]{0,1}C).*',str(dummy).split('\n')[0]):     
                it = 2 # alkene
            else:
                it =0 #unclassified



            '''
            if (len(smiles[smiles.name == item['name']].smiles)>0) and re.match('.*'+'[Cc]'+'.*',str(smiles[smiles.name == item['name']].smiles).split('\n')[0]): it+='c'        
            if (len(smiles[smiles.name == item['name']].smiles)>0) and re.match('.*'+'[Nn]'+'.*',str(smiles[smiles.name == item['name']].smiles).split('\n')[0]): it+='n' 
            '''    

            xloc = ( rad * np.cos(theta * counter) + xcen )
            yloc = ( rad * np.sin(theta * counter) + ycen )   


            string += '{"name":"%s","x":%s,"y":%s,"it":"%s","shape":"%s","group":%d,%s "alttext":"%s","with":"%s","w":%d,"id":%d },\n' %(item['name'],xloc,yloc,it,shape,0,size,item['Info'],item['reactions'][0],item['reactions'][1],nind[item['name']])


        string=string[:-2]
        string += '],"links":['

        ####links

        for i in graph.index:
            dat= graph.ix[i]

            val=''
            for run in runs:
                val+= ',"value%s":%e'%(run, dat['weight%d'%run])

            try: string += '{"source":%d,"target":%d%s},\n'%( nind[dat['Source']], nind[dat['Target']], val)    
            except Exception as e:          print e
        string=string[:-2]    
        string += ']}'    

     #re.sub(r"\n", r'',
        f.write(string)
        print 'JSON file written', f
        f.close()
        #os.system("perl -i-p -e 's/nan/0/g' %s"%filename)



#read all file permissions     
os.system('chmod 666 *.json')





if barplot: 



    def barchart(item):
            global reactants, products,getspec,r_loop,graph_list,concentrate,getspec ,getcoeff
            e = []
                for i in item:
                    r_coeff =  np.vectorize(getcoeff)(reactants[i]) #get react coeff
                        reactants[i] =  np.vectorize(getspec)(reactants[i]) #rm react coeff
                        products[i]  =  np.vectorize(getspec)(products[i] )

                         #rm prod coeff    
                        r_conc  =  np.vectorize( concentrate.__getitem__)(reactants[i])#get conc @ Timestp
                        t_flux = np.prod(r_coeff*r_conc)*r_loop[i]
                        ######## get edge reactions #######
                        if t_flux>0 : 
                            e.append([reactants[i], products[i],t_flux])

                        del r_conc, t_flux, r_coeff
                return e

    dummy=[]
    edges = multiprocessing.Pool(ncores).map(barchart, np.array_split(range(r_len),ncores)) 
    [dummy.extend(i) for i in edges]




    edges=[]

    for sp in shead:

        for i in dummy: 
            if (sp in i[0]): 
                if not (sp in i[1]):

                    keep=True
                    for j in [0,1]: 
                        for k in i[j]:
                            keep *= k in shead
                    if keep != 0: 
                        edges.append( ['-->'.join(['+'.join(i[0]), '+'.join(i[1])]), np.log10(i[2]),-1,sp])#react

            if (sp in i[1]):

                keep=True
                    for j in [0,1]: 
                        for k in i[j]:
                            keep *= k in shead
                    if keep != 0: 
                        edges.append( ['-->'.join(['+'.join(i[0]), '+'.join(i[1])]), np.log10(i[2]),+1,sp])#prod

    edges=np.array(edges)
    nmin = abs(min(edges[:,1].astype(float)))
    edges[:,1] = edges[:,1].astype(float) + nmin

    edges[:,1]= edges[:,1].astype(float) * edges[:,2].astype(float)
    ##edges.columns = ['reaction', 'flux', 'inout', 'inquestion']  

df=pd.DataFrame(edges)       


print 'string bar'

rs ='{ "name": "reactions",  "children": [ '
ps ='{ "name": "reactions",  "children": [ '

total ='{ "name": "reactions",  "children": [ '


for sp in shead:

    df1 = df[df[3]==sp]
    print df1[1].sum()





    if len(df1)>0: 
        rs+=' {  "name": "%s", "children": [  '%sp 

        for i in df1.iterrows():

            rs += '{"name": "%s", "size": %s},'%(i[1][0], i[1][1])

        rs=rs[:-1]    
        rs+=']},'
    else:
        print 'no reactions: ', sp 

rs=rs[:-1]                         
rs+=']}'  


f = open('bar1.json','w')
f.write(rs)
f.close()

os.system('chmod 666 bar1.json')









