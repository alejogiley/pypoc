# load libraries
import json
import requests

from urllib.error import URLError
from urllib.request import Request, urlopen

# add test to function
def find_json(json, name):
    """Read json objects and
       return value from item.
       
       Keyword arguments:
       name -- name of item
       json -- JSON object
    """
    try:
        mylist = []
        for item in json:
            if item.get(name) is not None:
                mylist += item.get(name)
        return mylist
    
    except TypeError:
        return [item.get(name) for item in json]

# add test to function 
def find_urls(json, keyy):
    """Read json objects and
       retrieve raw experimental 
       data *url* from PDB entry
       
       Keyword arguments:
       keyy -- name of PDB entry 
       json -- JSON object
    """
    # basic urls
    emdb_url = "ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/"
    xray_url = "http://www.ebi.ac.uk/pdbe/coordinates/files/"
    bnmr_url = "http://www.ebi.ac.uk/pdbe/entry-files/download/"
    
    # read json object
    for item in json:
        
        # Electron microscopy entries
        if item.get("experimental_method_class") == "em":
            # retrieve EMDB ID
            emdb_id = [ sub.get("emdb_id") 
                        for sub in item.get("additional_experimental_details") ]
            # convert e.g. "EMD-5776" into "emd_5776"
            return [ emdb_url + emdb_id[0] + "/map/" + 
                     emdb_id[0].lower().replace("-","_") + ".map.gz" ]
        
        # X-ray entries
        elif item.get("experimental_method_class") == "x-ray":
            return [ xray_url + keyy + ".ccp4" ]
        
        # NMR entries
        elif item.get("experimental_method_class") == "nmr":
            return [ bnmr_url + keyy + ".mr" ]
        
        else:
            return "None"

# add test to function 
def check_urls(url):
    """Handle ERROR response code 
       from server using urllib
       
       Keyword arguments:
       url -- web page link
    """
    req = Request(url)
    
    try:
        response = urlopen(req)
    except URLError as e:
        if hasattr(e, 'reason'):
            print('We failed to reach server %s' % url)
            print('Reason: ', e.reason)
        elif hasattr(e, 'code'):
            print('The server %s couldn\'t fulfill the request' % url)
            print('Error code: ', e.code)
    else:
        pass

# download link
def wget_urls(url):
    """Download http/ftp links
       save files on current path
       
       Keyword arguments:
       urls -- web link http/ftp
    """
    name = url.split("/")[-1]
    # this should work well with ftp links
    try:
        urlretrieve(url, name)
        print('File %s downloaded' % name)
    except:
        print('Error during download of %s' % url)

# reparameterization trick
# instead of sampling from Q(z|X), sample eps = N(0,I)
# z = z_mean + sqrt(var)*eps
from keras import backend as K

def sampling(args):
    """Reparameterization trick by sampling fr an isotropic unit Gaussian.

    # Arguments
        args (tensor): mean and log of variance of Q(z|X)

    # Returns
        z (tensor): sampled latent vector
    """

    z_mean, z_log_var = args
    batch = K.shape(z_mean)[0]
    dim = K.int_shape(z_mean)[1]
    # by default, random_normal has mean=0 and std=1.0
    epsilon = K.random_normal(shape=(batch, dim))
    return z_mean + K.exp(0.5 * z_log_var) * epsilon
