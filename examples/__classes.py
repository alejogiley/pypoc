class dPCAplus:
    '''Principal Component Analysis of dihedral angles
    
    Implementation of the dPCA+ method from references:
     
    J Comput Chem. 2009; 30(3):479-92. doi: 10.1002/jcc.21076.
    J Chem Phys. 2017; 147(24):244101. doi: 10.1063/1.4998259.
        
    with an Optimization step added to maximize shifting 
    1) the circular mean of the dihedral angles distribution 
    to the center of the space and 2) the maximal gap of the 
    sampling to the periodic boundary. 
    
    '''
    
    def __init__(self, angles, frames):
        '''Initialize class
        
        Parameters
        ----------
        angles : float
            numpy array of dihedral angles
            
        frames : int
            number of trajectory frames
        
        '''
        self.ang = angles
        self.fms = frames
        self.shp = angles.shape
    
    def lowdb(self, width=5):
        '''Estimate lowest density bin
        
        Parameters
        ----------
        width : int
            binwidth of the density 
            distribution (default is 5)
            
        Returns
        -------
        angle : float
            angle bin of lowest density
        
        '''
        # define pi
        pi = np.pi
        # binwidth to radians
        bw = np.deg2rad(width)
        # bin edges
        ed = np.arange(-pi,pi,bw)
        # calculate density
        dens, bx, by = np.histogram2d(self.ang[:,0],
                                      self.ang[:,1],
                                      bins=(ed, ed), 
                                      density=True)
        # get projections
        a, b = np.sum(dens, axis=0), np.sum(dens, axis=1)
        # get lowest density value
        c, d = np.sort(a)[0], np.sort(b)[0]
        # index of lowest density
        i, j = np.where(a == c), np.where(b == d)
        # return first bin of lowest density
        return np.array([ed[j][0], ed[i][0]])
    
    def swrap(self, rho):
        '''Shift and wrap angles
        
        shift angles by rho
        and wrap them around [-pi, pi]
        
        Parameters
        ----------
        rho : float
            shifts in radians
            rho.shape = 2
        
        '''
        # shift
        self.ang -= rho.T
        # wrap & rescale
        self.ang %= np.array([2*np.pi, 2*np.pi])
        self.ang -= np.array([np.pi, np.pi])
        
    def shape(self, rows, cols):
        angles = self.ang
        #print(type(angles))
        return angles.reshape(rows, cols)
    
    def optimize(self):
        '''Optimize the maximal smapling gap 
           at the Ramachandran map borders
        
        '''
        # data dimensions
        size  = np.prod(self.shp)
        size /= 2
        
        # format data
        self.ang = self.shape(int(size), 2)
        
        # shift low density to boundary 
        self.swrap(self.lowdb(5))
    
    def dopca(self):
        '''sklearn standard PCA
        '''
        # optimize
        self.optimize()
	
      	# data dimensions
        size  = np.prod(self.shp)
        size /= self.fms
        
        # format data
        data = self.shape(int(self.fms), 
                          int(size))

        # standardize data
        #scaler = StandardScaler()
        #scdata = scaler.fit_transform(data)
        
        # get PCA projection 
        pca = PCA(n_components=100)
        vec = pca.fit_transform(data)
        
        #return results
        return pca, vec
    
    @classmethod
    def cirmean(self):
        '''Calculate circular mean
            
        Returns
        -------
        mean : float
            circular mean of angles
        
        '''
        # transformation
        c = np.average(np.cos(self.ang))
        s = np.average(np.sin(self.ang))
        # circular mean
        return atan2(s, c)
    
    @classmethod
    def stdmean(self):
        '''Compute standard mean
        
        Returns
        -------
        mean : float
            standard mean of angles
        
        '''
        return np.average(self.ang)
    
    @staticmethod
    def dangle(x, y):
        '''Delta angle
         
        Parameters
        ----------
         
        x: float
           target torsion angle
        y: float 
           source torsion angle
           
        Returns
        -------
        angle : float
            minumum angle 
            between x and y
        
        '''
        s = sin(x-y) 
        c = cos(x-y)
        return atan2(s,c)
#END
