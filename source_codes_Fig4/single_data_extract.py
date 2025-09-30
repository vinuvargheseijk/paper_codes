import pandas as pd
import numpy as np
import dendShape
import writeXML
import matplotlib.pyplot as plt
import totalEnergy as TE
import tag

phi_entire = 9000
num_spines = 1
scalF = 0.65
initial_phitot = 5000
spacing = 5e-6
p_code = 2
perturb_num = 200
tag_string = tag.tag_string


class data_extract():
  """
  This class extracts the features of the spine specified
  """
  def __init__(self, ns):
    """
    Takes in the number of spine as the parameter
    """
    self.ns = ns    
    self.phi_list = []
    self.theta_list = []
    self.Hmean_list = []
    self.rm_list = []
    self.rp_list = []
    self.mConc_list = []
    self.cConc_list = []
    self.data_save_lists()

  def data_save_lists(self):
    """
    Gives out the values of phi, rm, theta, hmean, mConc, cConc in that order of the specified spine
    """
    phi = pd.read_csv(tag_string + "phi.csv")
    theta = pd.read_csv(tag_string + "theta.csv")
    hmean = pd.read_csv(tag_string + "hmean.csv")
    mConc = pd.read_csv(tag_string + "mConc.csv")
    cConc = pd.read_csv(tag_string + "cConc.csv")
    Rp = pd.read_csv(tag_string + "rp.csv")




    y_save = []
    x_save = []
    for i in range(len(phi)):
      if np.isnan(theta.iloc[i][self.ns]) == False:
          self.theta_list.append( theta.iloc[i][self.ns] )
          self.Hmean_list.append( hmean.iloc[i][self.ns] )
          self.rm_list.append( 1 / hmean.iloc[i][self.ns] )
          self.rp_list.append( Rp.iloc[i][self.ns] )
          self.mConc_list.append( mConc.iloc[i][self.ns] )
          self.cConc_list.append( cConc.iloc[i][self.ns] )
          self.phi_list.append( phi.iloc[i][self.ns] )
    plot = False      
    if plot == True:
      fig, ( (ax1, ax2), (ax3, ax4) ) = plt.subplots(2, 2)
      ax1.set_title("$\phi$")
      ax2.set_title("$R^{d}$")
      ax3.set_title("$\Theta$")

      ax1.plot( phi_list )
      ax2.plot( rm_list )
      ax3.plot( theta_list )
      ax4.plot( Hmean_list )
      plt.show()
    #return phi_list, rm_list, theta_list, Hmean_list, mConc_list, cConc_list


#ds = data_extract(1)
#phi_list, rm_list, theta_list, Hmean_list = ds.data_save_lists()

