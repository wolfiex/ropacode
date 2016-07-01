from ropa_lib import *
    

spec_list = pd.read_csv('mcm_subset_mass.txt', skiprows = 19 , delimiter = '~' , names = ['item','smiles','inchi','mass'])
inorganics= 'O,O3,O1D,O2,OH,NO,NO2,NO2N2O5,H2O2,HO2,HO2NO2,HONO,HNO3,CO,SO2,SO3,NA,NO3,HSO3,N2O5,RO2,H2,CL,A,H2O'.split(',')

inorg = pd.DataFrame(np.zeros((len(inorganics),4)), columns =  ['item','smiles','inchi','mass'])
inorg.item = inorganics
inorg.smiles = inorganics

spec_list = pd.concat([spec_list, inorg])
del inorganics, inorg 



spec_list = spec_list.sort('mass')

n,r=len(spec_list),700.
theta = 360. / n
locs = np.array([ [r*np.cos(theta*i), r*np.sin(theta*i)] for i in xrange(n)]) #get points around a circle
indexes = np.array([i for i in [range(0+j,n,int(n/8)) for j in xrange(0,int(n/8))]])
idx=[]
for i in indexes: idx.extend(i)
    
spec_list['id'] = idx 

spec_list = spec_list.sort('id') #seperate by mass - spaced

spec_list['x']= locs[:,0]
spec_list['y']= locs[:,1]

spec_list = spec_list.sort('item')

spec_list.index = spec_list['item']

spec_list.to_csv('fullmcmspecs.csv')



     
